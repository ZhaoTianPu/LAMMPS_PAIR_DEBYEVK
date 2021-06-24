/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Tianpu Zhao (TPZ)
   tz1416@ic.ac.uk, pacosynthesis@gmail.com
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Coulombic potential with various Debye screening and cutoff between 
   species pairs. Modified based on pair_lj_cut.cpp and pair_coul_debye.cpp.
   The essential idea is to make the kappa and cutoff in the original pair
   coul/debye a pair-specific parameter rather than a global one.
   All modifications from pair_coul_debye.cpp are commented.

   2021.05.07 Created                                                 TPZ
              (progress: finished compute, allocate
              to be finished: settings, coeff, init_style, init_one,
              write_restart, read_restart, write_restart_settings,
              read_restart_settings, write_data, write_data_all, single
              extract)
   2021.05.11 Modified                                                TPZ
              (progress: finished compute, allocate, settings, coeff,
              init_style, init_one, write_restart, read_restart, 
              write_restart_settings, read_restart_settings, 
              write_data, write_data_all, single, extract
              to be finished: test and debug)
   2021.05.31 Modified                                                TPZ
              (debug for single())
   2021.06.01 Modified                                                TPZ
              (changed line 64 from Pair to PairCoulCut)
   2021.06.10 Modified                                                TPZ
              (remove scale)
------------------------------------------------------------------------- */

#include "pair_coul_debye_vk.h"

#include <cmath>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "math_const.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairCoulDebyeVK::PairCoulDebyeVK(LAMMPS *lmp) : Pair(lmp) {} 

/* ---------------------------------------------------------------------- */

// add destroy (and allocate) memory

PairCoulDebyeVK::~PairCoulDebyeVK()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(kappa); 
  }
}

/* ---------------------------------------------------------------------- */

void PairCoulDebyeVK::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair;
  double rsq,r2inv,r,rinv,forcecoul,factor_coul,screening;
  int *ilist,*jlist,*numneigh,**firstneigh;

  ecoul = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        r = sqrt(rsq);
        rinv = sqrt(r2inv);
        // calculate screening for each pair
        screening = exp(-kappa[itype][jtype]*r);
        // calculate force with each individual kappa
        forcecoul = qqrd2e * 
          qtmp*q[j] * screening * (kappa[itype][jtype] + rinv);
        fpair = factor_coul*forcecoul * r2inv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag)
          // calculate energy
          ecoul = factor_coul * qqrd2e *
            qtmp*q[j] * rinv * screening;

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             0.0,ecoul,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

// add allocate memory here

void PairCoulDebyeVK::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(kappa,n+1,n+1,"pair:kappa");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

/* one argument, global cutoff */

void PairCoulDebyeVK::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = utils::numeric(FLERR,arg[0],false,lmp);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

/* 3 parameters, i j kappa */

void PairCoulDebyeVK::coeff(int narg, char **arg)
{
  if (narg < 3 || narg > 4)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();
  
  // processing species numbers
  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);
  
  // kappa
  double kappa_one = utils::numeric(FLERR,arg[2],false,lmp);

  // pair cutoff setting, if arg[3] is there, then specify it for pair
  double cut_one = cut_global;
  if (narg == 4) cut_one = utils::numeric(FLERR,arg[3],false,lmp);
  
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      kappa[i][j] = kappa_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

/* just to check if the species has charge or not */

void PairCoulDebyeVK::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style coul/cut requires atom attribute q");

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairCoulDebyeVK::init_one(int i, int j)
{
  if (setflag[i][j] == 0){
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
    // add kappa[i][j]
    kappa[i][j] = mix_distance(kappa[i][i],kappa[j][j]);
  }

  

  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulDebyeVK::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        // set cut and kappa
        fwrite(&cut[i][j],sizeof(double),1,fp);
        fwrite(&kappa[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulDebyeVK::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
          // add kappa
          utils::sfread(FLERR,&kappa[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
        // add kappa
        MPI_Bcast(&kappa[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

/* no kappa or tail, just mix_flag, cut_global and offset_flag */

void PairCoulDebyeVK::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulDebyeVK::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

/* write i kappa */

void PairCoulDebyeVK::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g\n",i,kappa[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

/* write i j kappa cut */

void PairCoulDebyeVK::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g\n",i,j,kappa[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

/* calculate single pair potential */

double PairCoulDebyeVK::single(int i, int j, int /*itype*/, int /*jtype*/,
                           double rsq, double factor_coul, double /*factor_lj*/,
                           double &fforce)
{
  double r2inv,rinv,forcecoul,phicoul,r,screening;

  r2inv = 1.0/rsq;
  r = sqrt(rsq);
  rinv = 1.0/r;
  screening = exp(-kappa[i][j]*r);
  forcecoul = force->qqrd2e * atom->q[i]*atom->q[j] *
    screening * (kappa[i][j] + rinv);
  fforce = factor_coul*forcecoul * r2inv;

  phicoul = force->qqrd2e * atom->q[i]*atom->q[j] * rinv * screening;
  return factor_coul*phicoul;
}

/* ---------------------------------------------------------------------- */

void *PairCoulDebyeVK::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"cut_coul") == 0) return (void *) &cut;
  if (strcmp(str,"kappa") == 0) return (void *) kappa;
  return nullptr;
}
