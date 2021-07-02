### pair_style coul/debye/vk command

#### Syntax

```
pair_style coul/debye/vk global_cutoff
pair_coeff i j kappa cutoff
```

* global_cutoff = global cutoff for Debye-shielded Coulombic interactions
* i, j = atom types
* kappa = (inverse) Debye length
* cutoff = (optional) cutoff between individual pairs

#### Examples

```
pair_style coul/debye/vk 5.0
pair_coeff * * 1.0 

pair_style coul/debye/vk 5.0
pair_coeff 1 1 1.0
pair_coeff 2 2 2.0 2.5
```

#### Description

The style *coul/debye/vk* adds an additional exp() damping factor to the Coulombic term, given by
$$
E = \frac{Cq_iq_j}{\epsilon r}\exp(-\kappa_{ij}r)\quad r < r_{c,ij}\,,
$$
where $C$ is an energy-conversion constant, $q_i$ and $q_j$ are the charges on the 2 atoms, and $\epsilon$ is the dielectric constant which can be set by the [dielectric](https://docs.lammps.org/dielectric.html) command. The cutoff $r_{c,ij}$ truncates the interaction distance.  

Notice that this style differs from *coul/debye* in the Debye screening length: the inverse Debye screening length $\kappa$ is a global parameter in *coul/debye*, i.e. it applies for all pairs, whereas in *coul/debye/vk*, each pair is assigned with an individual $\kappa_{ij}$ value.

At least all the $\kappa_{ij}$ values for i = j need to be set; for i != j, the inverse screening length can be automatically obtained by applying mixing rules. The default mixing rule is *geometric*, i.e. $\kappa_{ij} = \sqrt{\kappa_{ii}\kappa_{jj}}$.

The cutoff parameter in pair_coeff command is optional. If it is not used (as in some of the examples above), the default global value specified in the pair_style command is used.

#### Mixing, shift, table, tail correction, restart, rRESPA info

For atom type pairs i,j and i != j, the cutoff distance and kappa can be mixed. The default mix value for both are *geometric*. See the “pair_modify” command for details.

The [pair_modify](https://docs.lammps.org/pair_modify.html) shift option is not relevant for these pair styles.

This pair style do not support the [pair_modify](https://docs.lammps.org/pair_modify.html) tail option for adding long-range tail corrections to energy and pressure.

These pair styles write their information to [binary restart files](https://docs.lammps.org/restart.html), so pair_style and pair_coeff commands do not need to be specified in an input script that reads a restart file.

These pair styles can only be used via the *pair* keyword of the [run_style respa](https://docs.lammps.org/run_style.html) command. They do not support the *inner*, *middle*, *outer* keywords.

