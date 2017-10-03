/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <stdlib.h>
#include <string.h>
#include "fix_wall_hydrophobic.h"
#include "atom.h"
#include "comm.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "lattice.h"
#include "input.h"
#include "variable.h"
#include "error.h"
#include "force.h"
#include <math.h>
#include "respa.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{XLO=0,XHI=1,YLO=2,YHI=3,ZLO=4,ZHI=5};
enum{NONE=0,EDGE,CONSTANT,VARIABLE};


/* ---------------------------------------------------------------------- */

FixWallHydrophobic::FixWallHydrophobic(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  nwall(0)
{
  // The inputs after the xlo keyword are position epsilon sigma lambda cutoff
  scalar_flag = 1;
  vector_flag = 1;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  respa_level_support = 1;
  ilevel_respa = 0;
  virial_flag = 1;

  // parse args

  int scaleflag = 1;
  fldflag = 0;
  int pbcflag = 0;

  for (int i = 0; i < 6; i++) xstr[i] = estr[i] = sstr[i] = NULL;

  int iarg = 3;
  while (iarg < narg) {
    if ((strcmp(arg[iarg],"xlo") == 0) || (strcmp(arg[iarg],"xhi") == 0) ||
        (strcmp(arg[iarg],"ylo") == 0) || (strcmp(arg[iarg],"yhi") == 0) ||
        (strcmp(arg[iarg],"zlo") == 0) || (strcmp(arg[iarg],"zhi") == 0)) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal fix wall command");

      int newwall;
      if (strcmp(arg[iarg],"xlo") == 0) newwall = XLO;
      else if (strcmp(arg[iarg],"xhi") == 0) newwall = XHI;
      else if (strcmp(arg[iarg],"ylo") == 0) newwall = YLO;
      else if (strcmp(arg[iarg],"yhi") == 0) newwall = YHI;
      else if (strcmp(arg[iarg],"zlo") == 0) newwall = ZLO;
      else if (strcmp(arg[iarg],"zhi") == 0) newwall = ZHI;

      for (int m = 0; (m < nwall) && (m < 6); m++)
        if (newwall == wallwhich[m])
          error->all(FLERR,"Wall defined twice in fix wall command");

      wallwhich[nwall] = newwall;
      if (strcmp(arg[iarg+1],"EDGE") == 0) {
        xstyle[nwall] = EDGE;
        int dim = wallwhich[nwall] / 2;
        int side = wallwhich[nwall] % 2;
        if (side == 0) coord0[nwall] = domain->boxlo[dim];
        else coord0[nwall] = domain->boxhi[dim];
      } else if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        xstyle[nwall] = VARIABLE;
        int n = strlen(&arg[iarg+1][2]) + 1;
        xstr[nwall] = new char[n];
        strcpy(xstr[nwall],&arg[iarg+1][2]);
      } else {
        xstyle[nwall] = CONSTANT;
        coord0[nwall] = force->numeric(FLERR,arg[iarg+1]);
      }

      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) {
        int n = strlen(&arg[iarg+2][2]) + 1;
        estr[nwall] = new char[n];
        strcpy(estr[nwall],&arg[iarg+2][2]);
        estyle[nwall] = VARIABLE;
      } else {
        epsilon[nwall] = force->numeric(FLERR,arg[iarg+2]);
        estyle[nwall] = CONSTANT;
      }

      if (strstr(arg[iarg+3],"v_") == arg[iarg+3]) {
        int n = strlen(&arg[iarg+3][2]) + 1;
        sstr[nwall] = new char[n];
        strcpy(sstr[nwall],&arg[iarg+3][2]);
        sstyle[nwall] = VARIABLE;
      } else {
        sigma[nwall] = force->numeric(FLERR,arg[iarg+3]);
        sstyle[nwall] = CONSTANT;
      }

      lambda[nwall] = force->numeric(FLERR,arg[iarg+4]);
      cutoff[nwall] = force->numeric(FLERR,arg[iarg+5]);
      nwall++;
      iarg += 6;

    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix wall command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal fix wall command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"fld") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix wall command");
      if (strcmp(arg[iarg+1],"no") == 0) fldflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) fldflag = 1;
      else error->all(FLERR,"Illegal fix wall command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"pbc") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix wall command");
      if (strcmp(arg[iarg+1],"yes") == 0) pbcflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) pbcflag = 0;
      else error->all(FLERR,"Illegal fix wall command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix wall command");
  }

  size_vector = nwall;

}

/* ---------------------------------------------------------------------- */

void FixWallHydrophobic::precompute(int m)
{
  coeff1[m] = 48.0 * epsilon[m] * pow(sigma[m],12.0);
  coeff2[m] = 24.0 * epsilon[m] * pow(sigma[m],6.0);
  coeff3[m] = 4.0 * epsilon[m] * pow(sigma[m],12.0);
  coeff4[m] = 4.0 * epsilon[m] * pow(sigma[m],6.0);

  double r2inv = 1.0/(cutoff[m]*cutoff[m]);
  double r6inv = r2inv*r2inv*r2inv;
  offset[m] = r6inv*(coeff3[m]*r6inv - coeff4[m]);
  z0[m] = pow( 2. , 1./6.)*sigma[m];
}

/* ----------------------------------------------------------------------
   interaction of all particles in group with a wall
   m = index of wall coeffs
   which = xlo,xhi,ylo,yhi,zlo,zhi
   error if any particle is on or behind wall
------------------------------------------------------------------------- */

void FixWallHydrophobic::wall_particle(int m, int which, double coord)
{
  double delta,rinv,r2inv,r6inv,fwall,W0,W1,fW0,fW1;
  double vn;

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int dim = which / 2;
  int side = which % 2;
  if (side == 0) side = -1;

  int onflag = 0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (side < 0) delta = x[i][dim] - coord;
      else delta = coord - x[i][dim];
      if (delta >= cutoff[m]) continue;
      if (delta <= 0.0) {
        onflag = 1;
        continue;
      }
      rinv = 1.0/delta;
      r2inv = rinv*rinv;
      r6inv = r2inv*r2inv*r2inv;
      if (delta<z0[m]) {
         fW0 = side * r6inv*(coeff1[m]*r6inv - coeff2[m]) * rinv;
         fW1 = 0.;
         W0 = r6inv*(coeff3[m]*r6inv - coeff4[m]) + epsilon[m];
         W1 = - epsilon[m];
      } else {
         fW0 = 0.;
         fW1 = side * r6inv*(coeff1[m]*r6inv - coeff2[m]) * rinv;
         W0 = 0;
         W1 = r6inv*(coeff3[m]*r6inv - coeff4[m]);
      }
      fwall = fW0 + lambda[m]*fW1;
      f[i][dim] -= fwall;
      ewall[0] += W0 + lambda[m]*W1;
      ewall[m+1] += fwall;

      if (evflag) {
        if (side < 0) vn = -fwall*delta;
        else vn = fwall*delta;
        v_tally(dim, i, vn);
      }
    }

  if (onflag) error->one(FLERR,"Particle on or inside fix wall/hydrophobic surface");
}

/* ---------------------------------------------------------------------- */

FixWallHydrophobic::~FixWallHydrophobic()
{
  if (copymode) return;

  for (int m = 0; m < nwall; m++) {
    delete [] xstr[m];
    delete [] estr[m];
    delete [] sstr[m];
  }
}

/* ---------------------------------------------------------------------- */

int FixWallHydrophobic::setmask()
{
  int mask = 0;

  // FLD implicit needs to invoke wall forces before pair style

  if (fldflag) mask |= PRE_FORCE;
  else mask |= POST_FORCE;

  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallHydrophobic::init()
{
  for (int m = 0; m < nwall; m++) {
    if (xstyle[m] == VARIABLE) {
      xindex[m] = input->variable->find(xstr[m]);
      if (xindex[m] < 0)
        error->all(FLERR,"Variable name for fix wall/hydrophobic does not exist");
      if (!input->variable->equalstyle(xindex[m]))
        error->all(FLERR,"Variable for fix wall/hydrophobic is invalid style");
    }
    if (estyle[m] == VARIABLE) {
      eindex[m] = input->variable->find(estr[m]);
      if (eindex[m] < 0)
        error->all(FLERR,"Variable name for fix wall/hydrophobic does not exist");
      if (!input->variable->equalstyle(eindex[m]))
        error->all(FLERR,"Variable for fix wall/hydrophobic is invalid style");
    }
    if (sstyle[m] == VARIABLE) {
      sindex[m] = input->variable->find(sstr[m]);
      if (sindex[m] < 0)
        error->all(FLERR,"Variable name for fix wall/hydrophobic does not exist");
      if (!input->variable->equalstyle(sindex[m]))
        error->all(FLERR,"Variable for fix wall/hydrophobic is invalid style");
    }
  }

  // setup coefficients

  for (int m = 0; m < nwall; m++) precompute(m);

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixWallHydrophobic::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet")) {
    if (!fldflag) post_force(vflag);
  } else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixWallHydrophobic::min_setup(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   only called if fldflag set, in place of post_force
------------------------------------------------------------------------- */

void FixWallHydrophobic::pre_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWallHydrophobic::post_force(int vflag)
{

  // energy and virial setup

  eflag = 0;
  if (vflag) v_setup(vflag);
  else evflag = 0;
  for (int m = 0; m <= nwall; m++) ewall[m] = 0.0;

  // coord = current position of wall
  // evaluate variables if necessary, wrap with clear/add
  // for epsilon/sigma variables need to re-invoke precompute()

  if (varflag) modify->clearstep_compute();

  double coord;
  for (int m = 0; m < nwall; m++) {
    if (xstyle[m] == VARIABLE) {
      coord = input->variable->compute_equal(xindex[m]);
      if (wallwhich[m] < YLO) coord *= xscale;
      else if (wallwhich[m] < ZLO) coord *= yscale;
      else coord *= zscale;
    } else coord = coord0[m];
    if (wstyle[m] == VARIABLE) {
      if (estyle[m] == VARIABLE) {
        epsilon[m] = input->variable->compute_equal(eindex[m]);
        if (epsilon[m] < 0.0)
          error->all(FLERR,"Variable evaluation in fix wall/hydrophobic gave bad value");
      }
      if (sstyle[m] == VARIABLE) {
        sigma[m] = input->variable->compute_equal(sindex[m]);
        if (sigma[m] < 0.0)
          error->all(FLERR,"Variable evaluation in fix wall/hydrophobic gave bad value");
      }
      precompute(m);
    }

    wall_particle(m,wallwhich[m],coord);
  }

  if (varflag) modify->addstep_compute(update->ntimestep + 1);
}

/* ---------------------------------------------------------------------- */

void FixWallHydrophobic::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWallHydrophobic::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy of wall interaction
------------------------------------------------------------------------- */

double FixWallHydrophobic::compute_scalar()
{
  // only sum across procs one time

  if (eflag == 0) {
    MPI_Allreduce(ewall,ewall_all,nwall+1,MPI_DOUBLE,MPI_SUM,world);
    eflag = 1;
  }
  return ewall_all[0];
}

/* ----------------------------------------------------------------------
   components of force on wall
------------------------------------------------------------------------- */

double FixWallHydrophobic::compute_vector(int n)
{
  // only sum across procs one time

  if (eflag == 0) {
    MPI_Allreduce(ewall,ewall_all,nwall+1,MPI_DOUBLE,MPI_SUM,world);
    eflag = 1;
  }
  return ewall_all[n+1];
}
