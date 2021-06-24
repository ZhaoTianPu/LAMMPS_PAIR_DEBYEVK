/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------
   Coulombic potential with various Debye screening and cutoff between 
   species pairs. Modified based on pair_lj_cut.cpp and pair_coul_debye.cpp.
   The essential idea is to make the kappa and cutoff in the original pair
   coul/debye a pair-specific parameter rather than a global one.
   All modifications from pair_coul_debye.cpp are commented.

   2021.05.07 Created                                                 TPZ
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(coul/debye/vk,PairCoulDebyeVK)

#else

#ifndef LMP_PAIR_COUL_DEBYE_VK_H
#define LMP_PAIR_COUL_DEBYE_VK_H

// #include "pair.h"
#include "pair.h"

namespace LAMMPS_NS {

class PairCoulDebyeVK : public Pair {
 public:
  PairCoulDebyeVK(class LAMMPS *);
  virtual ~PairCoulDebyeVK();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);

 protected:
  double cut_global;
  double **cut;
  double **kappa;

  virtual void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style coul/cut requires atom attribute q

The atom style defined does not have these attributes.

*/
