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

#ifdef FIX_CLASS

FixStyle(wall/hydrophobic,FixWallHydrophobic)

#else

#ifndef LMP_FIX_WALL_HYDROPHOBIC_H
#define LMP_FIX_WALL_HYDROPHOBIC_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWallHydrophobic : public Fix {
 public:
  int nwall;
  int wallwhich[6];
  double coord0[6];
  int xflag;           // 1 if any wall position is a variable
  int xstyle[6];
  int xindex[6];
  char *xstr[6];

  FixWallHydrophobic(class LAMMPS *, int, char **);
  virtual ~FixWallHydrophobic();
  int setmask();
  void init();
  void precompute(int);
  void wall_particle(int, int, double);
  void setup(int);
  void min_setup(int);
  void pre_force(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);

 protected:
  int wallstyle[6];
  char *varstr[6];
  int varindex[6];
  double coeff1[6],coeff2[6],coeff3[6],coeff4[6],offset[6];
  double lambda[6], sigmaS[6];
  double epsilon[6],sigma[6],cutoff[6];
  double ewall[7],ewall_all[7];
  double xscale,yscale,zscale;
  int estyle[6],sstyle[6],wstyle[6];
  int eindex[6],sindex[6];
  char *estr[6],*sstr[6];
  int varflag;                // 1 if any wall position,epsilon,sigma is a var
  int eflag;                  // per-wall flag for energy summation
  int ilevel_respa;
  int fldflag;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Wall defined twice in fix wall/hydrophobic command

Self-explanatory.

E: Cannot use fix wall/hydrophobic in periodic dimension

Self-explanatory.

E: Cannot use fix wall/hydrophobic zlo/zhi for a 2d simulation

Self-explanatory.

E: Variable name for fix wall/hydrophobic does not exist

Self-explanatory.

E: Variable for fix wall/hydrophobic is invalid style

Only equal-style variables can be used.

W: Should not allow rigid bodies to bounce off relecting walls

LAMMPS allows this, but their dynamics are not computed correctly.

*/
