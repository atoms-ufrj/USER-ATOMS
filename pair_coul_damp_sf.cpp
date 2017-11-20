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

#include "string.h"
#include "pair_coul_damp_sf.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "math_const.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairCoulDampSF::PairCoulDampSF(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 1;
  self_flag = 0;
}

/* ---------------------------------------------------------------------- */

PairCoulDampSF::~PairCoulDampSF()
{
  if (!copymode)
    if (allocated)
      memory->destroy(setflag);
}

/* ---------------------------------------------------------------------- */

void PairCoulDampSF::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype,intra;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,vr,fr,evdwl,ecoul,fpair;
  double r,rsq,r2inv,r6inv,forcelj,prefactor,forcecoul,factor_lj,factor_coul;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  double *special_coul = force->special_coul;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // Compute self energy:
  if (eflag && self_flag)
    for (i = 0; i < nlocal; i++)
      ev_tally(i,i,nlocal,0,0.0,e_self*q[i]*q[i],0.0,0.0,0.0,0.0);

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = qqrd2e*q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      intra = sbmask(j);
      factor_lj = special_lj[intra];
      factor_coul = special_coul[intra];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cut_coulsq) {
        r2inv = 1.0/rsq;

        prefactor = factor_coul*qtmp*q[j];
        if (intra)
          forcecoul = prefactor*sqrt(r2inv);
        else {
          r = sqrt(rsq);
          unshifted( r, vr, fr );
          forcecoul = prefactor*(fr - f_shift)*r;
        }

        fpair = forcecoul*r2inv;
        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag)
          ecoul = intra ? forcecoul : prefactor*(vr + r*f_shift - e_shift);

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,ecoul,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairCoulDampSF::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairCoulDampSF::settings(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal pair_style command");

  alpha = force->numeric(FLERR,arg[0]);
  cut_coul = force->numeric(FLERR,arg[1]);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairCoulDampSF::coeff(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairCoulDampSF::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style coul/damp/sf requires atom charges");

  neighbor->request(this,instance_me);

  cut_coulsq = cut_coul * cut_coul;
  unshifted( cut_coul, e_shift, f_shift );
  e_shift += f_shift*cut_coul;
  e_self = -(e_shift/2.0 + alpha/sqrt(MY_PI))*force->qqrd2e;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairCoulDampSF::init_one(int i, int j)
{
  return cut_coul;
}

/* ---------------------------------------------------------------------- */

void PairCoulDampSF::modify_params(int narg, char **arg)
{
  if (narg == 0)
    error->all(FLERR,"Illegal pair_modify command");

  int iarg, ns, skip[narg];
  iarg = ns = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"self") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0)
        self_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0)
        self_flag = 0;
      else
        error->all(FLERR,"Illegal pair_modify command");
      single_enable = !self_flag;
      iarg += 2;
    }
    else // no keyword found - skip argument:
      skip[ns++] = iarg++;
  }

  // Call parent-class routine with skipped arguments:
  if (ns > 0) {
    for (int i = 0; i < ns; i++)
      arg[i] = arg[skip[i]];
    Pair::modify_params(ns, arg);
  }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulDampSF::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++)
      fwrite(&setflag[i][j],sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulDampSF::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulDampSF::write_restart_settings(FILE *fp)
{
  fwrite(&alpha,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&self_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulDampSF::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&alpha,sizeof(double),1,fp);
    fread(&cut_coul,sizeof(double),1,fp);
    fread(&self_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&alpha,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&self_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

double PairCoulDampSF::single(int i, int j, int itype, int jtype, double rsq,
                                double factor_coul, double factor_lj,
                                double &fforce)
{
  double r2inv,r6inv,r,vr,fr,prefactor;

  r2inv = 1.0/rsq;
  fforce = 0.0;
  if (rsq < cut_coulsq) {
    r = sqrt(rsq);
    unshifted( r, vr, fr );
    prefactor = factor_coul * force->qqrd2e * atom->q[i] * atom->q[j];
    fforce += prefactor*(fr-f_shift)*r;
  }
  fforce *= r2inv;

  double eng = 0.0;
  if (rsq < cut_coulsq)
    eng += prefactor * (vr + r*f_shift - e_shift);
 
  return eng;
}

/* ---------------------------------------------------------------------- */

void *PairCoulDampSF::extract(const char *str, int &dim)
{
  if (strcmp(str,"cut_coul") == 0) {
    dim = 0;
    return (void *) &cut_coul;
  }
  return NULL;
}

