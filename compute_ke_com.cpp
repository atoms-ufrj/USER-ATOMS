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

#include <mpi.h>
#include "compute_ke_com.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

#define BIG 2000000000

/* ---------------------------------------------------------------------- */

ComputeKECOM::ComputeKECOM(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute ke/com command");

  scalar_flag = 1;
  extscalar = 1;

  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;

  tagint mmax = -1;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) mmax = MAX(mmax,molecule[i]);
  MPI_Allreduce( &mmax, &maxmol, 1, MPI_LMP_TAGINT, MPI_MAX, world );
}

/* ---------------------------------------------------------------------- */

void ComputeKECOM::init()
{
  pfactor = 0.5 * force->mvv2e;
}

/* ---------------------------------------------------------------------- */

double ComputeKECOM::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  int nlocal = atom->nlocal;
  double **v = atom->v;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *mask = atom->mask;
  int *type = atom->type;
  tagint *molecule = atom->molecule;

  double data[maxmol][4];
  double all[maxmol][4];

  for (int i = 0; i < maxmol; i++)
    data[i][0] = data[i][1] = data[i][2] = data[i][3] = 0.0;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        int imol = molecule[i] - 1;
        double m = rmass[i];
        data[imol][0] += m*v[i][0];
        data[imol][1] += m*v[i][1];
        data[imol][2] += m*v[i][2];
        data[imol][3] += m;
      }
  }
  else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        int imol = molecule[i] - 1;
        double m = mass[type[i]];
        data[imol][0] += m*v[i][0];
        data[imol][1] += m*v[i][1];
        data[imol][2] += m*v[i][2];
        data[imol][3] += m;
      }
  }

  MPI_Allreduce( &data, &all, 4*maxmol, MPI_DOUBLE, MPI_SUM, world );

  double ke = 0.0;
  for (int i = 0; i < maxmol; i++)
    if (all[i][3] > 0.0)
      ke += (all[i][0]*all[i][0] + 
             all[i][1]*all[i][1] +
             all[i][2]*all[i][2])/all[i][3];

  scalar = pfactor*ke;
  return scalar;
}

/* ---------------------------------------------------------------------- */

