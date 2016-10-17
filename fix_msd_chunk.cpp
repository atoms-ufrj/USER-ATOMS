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

/* ----------------------------------------------------------------------
   Contributing authors: Charlles Abreu (abreu@eq.ufrj.br)
                         Ana Silveira (asilveira@plapiqui.edu.ar)
------------------------------------------------------------------------- */

#include <stdlib.h>
#include <string.h>
#include "fix_msd_chunk.h"
#include "update.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "atom.h"
#include "group.h"
#include "compute_chunk_atom.h"
#include "domain.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMSDChunk::FixMSDChunk(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 8) error->all(FLERR,"Illegal fix msd/chunk command");

  // Syntax: fix ID group msd/chunk chunkID nevery blocksize file Nfreq

  // Set a pointer to the corresponding compute chunk/atom:
  int icompute = modify->find_compute(arg[3]);
  if (icompute < 0 || strcmp(modify->compute[icompute]->style,"chunk/atom") != 0)
    error->all(FLERR,"ChunkID passed to fix msd/chunk is not a compute chunk/atom ID");
  cchunk = (ComputeChunkAtom *) modify->compute[icompute];

  nevery = force->inumeric(FLERR,arg[4]);
  blocksize = force->inumeric(FLERR,arg[5]);

  int n = strlen(arg[6]) + 1;
  filename = new char[n];
  strcpy(filename,arg[6]);

  nfreq = force->inumeric(FLERR,arg[7]);

  // Allocate arrays:
  nchunks = cchunk->setup_chunks();
  memory->create(com,nchunks,3,"fix_msd_chunk:com");
  allocate(1);

  // Set fix-related flags:
  array_flag = 1;                   // Array is computed (mean square displacements)
  size_array_cols = 5;              // Columns are [time msdx msdy msdz msd]
  size_array_rows_variable = 1;     // Number of rows is not knwon in advance
  size_array_rows = blocksize;
  extarray = 0;                     // Array elements are intensive variables
}

/* ---------------------------------------------------------------------- */

FixMSDChunk::~FixMSDChunk()
{
  memory->destroy(com);
  deallocate();
}

/* ---------------------------------------------------------------------- */

int FixMSDChunk::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMSDChunk::init()
{
  int ntimestep = update->nsteps;

  // Determine maximum number of levels in this run:
  nlevels = 1;
  int level = ntimestep/blocksize;
  while (level != 0) {
    nlevels++;
    level /= blocksize;
  }

  // If necessary, reallocate arrays:
  if (nlevels > maxlevel) {
    deallocate();
    allocate(nlevels);
  }

  // Initialize variables:
  step = 0;
  compute_centers_of_mass();
  int m = 0;
  msd[0][0] = msd[0][1] = msd[0][2] = msd[0][3] = msd[0][4] = 0.0;
  for (int level = 0; level < maxlevel; level++) {
    samples[level] = 0;
    for (int i = 0; i < nchunks; i++) {
      comtab[level][0][i][0] = com[i][0];
      comtab[level][0][i][1] = com[i][1];
      comtab[level][0][i][2] = com[i][2];
    }
    for (int jump = 1; jump < blocksize; jump++) {
      counter[level][jump] = 0;
      disptab[level][jump][0] = disptab[level][jump][1] = disptab[level][jump][2] = 0.0;
      msd[++m][0] = (double) (nevery*jump*ipow(blocksize,level));
    }
  }

  uptodate = 0;
  nlevels = 0;
}

/* ---------------------------------------------------------------------- */

void FixMSDChunk::end_of_step()
{
  int i, j, level, current, origin, jump, end, ini;
  double dx, dy, dz;

  step++;

  compute_centers_of_mass();

  nlevels = 1;
  level = step/blocksize;
  while (level != 0) {
    nlevels++;
    level /= blocksize;
  }

  for (level = 0; level < nlevels; level++)
    if (step % ipow(blocksize,level) == 0) {
      current = samples[level] + 1;
      end = current % blocksize;
      for (i = 0; i < nchunks; i++) {
        comtab[level][end][i][0] = com[i][0];
        comtab[level][end][i][1] = com[i][1];
        comtab[level][end][i][2] = com[i][2];
        for (origin = MAX(0,current-blocksize); origin < current; origin++) {
          jump = current - origin;
          counter[level][jump]++;
          ini = origin % blocksize;
          dx = com[i][0] - comtab[level][ini][i][0];
          dy = com[i][1] - comtab[level][ini][i][1];
          dz = com[i][2] - comtab[level][ini][i][2];
          disptab[level][jump][0] += dx*dx;
          disptab[level][jump][1] += dy*dy;
          disptab[level][jump][2] += dz*dz;
  	}
      }
      samples[level] = current;
    }

  uptodate = 0;

  if (update->ntimestep % nfreq == 0) {
    compute_array(0,0);
    FILE *fp = fopen(filename,"a");
    fprintf(fp,"timestep = %d:\n",update->ntimestep);
    fprintf(fp,"time msdx msdy msdz msd\n");
    int m = 0;
    for (int level = 0; level < nlevels; level++) {
      int levelsize = MIN(samples[level] + 1, blocksize);
      for (int jump = 1; jump < levelsize; jump++) {
        m++;
        fprintf(fp,"%f %f %f %f %f\n",msd[m][0],msd[m][1],msd[m][2],msd[m][3],msd[m][4]);
      }
    }
    fprintf(fp,"\n");
    fclose(fp);
  }
}

/* ---------------------------------------------------------------------- */

double FixMSDChunk::compute_array(int i, int j)
{
  if (!uptodate) {
    int m = 0;
    for (int level = 0; level < nlevels; level++) {
      int levelsize = MIN(samples[level] + 1, blocksize);
      for (int jump = 1; jump < levelsize; jump++) {
        m++;
        msd[m][1] = disptab[level][jump][0]/counter[level][jump];
        msd[m][2] = disptab[level][jump][1]/counter[level][jump];
        msd[m][3] = disptab[level][jump][2]/counter[level][jump];
        msd[m][4] = msd[m][1] + msd[m][2] + msd[m][3];
      }
    }
    size_array_rows = m + 1;
    uptodate = 1;
  }
  if (i < size_array_rows && j < size_array_cols)
    return msd[i][j];
  else
    return 0.0;
}

/* ---------------------------------------------------------------------- */

int FixMSDChunk::ipow(int base, int exp)
{
  int result = 1;
  while (exp) {
    if (exp & 1) result *= base;
    exp >>= 1;
    base *= base;
  }
  return result;
}

/* ---------------------------------------------------------------------- */

void FixMSDChunk::allocate(int nlevels)
{
  maxlevel = nlevels;
  memory->create(counter,maxlevel,blocksize,"fix_msd_chunk:counter");
  memory->create(samples,maxlevel,"fix_msd_chunk:samples");
  memory->create(comtab,maxlevel,blocksize,nchunks,3,"fix_msd_chunk:comtab");
  memory->create(disptab,maxlevel,blocksize,3,"fix_msd_chunk:disptab");
  memory->create(msd,maxlevel*(blocksize-1)+1,5,"fix_msd_chunk:msd");
}

/* ---------------------------------------------------------------------- */

void FixMSDChunk::deallocate()
{
  memory->destroy(counter);
  memory->destroy(samples);
  memory->destroy(comtab);
  memory->destroy(disptab);
  memory->destroy(msd);
}

/* ---------------------------------------------------------------------- */

void FixMSDChunk::compute_centers_of_mass()
{
  int index;
  double massone;
  double unwrap[3];

  int n = cchunk->setup_chunks();
  if (n != nchunks) error->all(FLERR,"Fix msd/chunk nchunks is not static");
  cchunk->compute_ichunk();
  int *ichunk = cchunk->ichunk;

  double data[nchunks][4];

  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nchunks; i++)
    data[i][0] = data[i][1] = data[i][2] = data[i][3] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      index = ichunk[i]-1;
      if (index < 0) continue;
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      domain->unmap(x[i],image[i],unwrap);
      data[index][0] += unwrap[0] * massone;
      data[index][1] += unwrap[1] * massone;
      data[index][2] += unwrap[2] * massone;
      data[index][3] += massone;
    }

  double alldata[nchunks][4];
  MPI_Allreduce(&data[0][0],&alldata[0][0],4*nchunks,MPI_DOUBLE,MPI_SUM,world);

  for (int i = 0; i < nchunks; i++) {
    massone = alldata[i][3];
    if (massone > 0.0) {
      com[i][0] = alldata[i][0]/massone;
      com[i][1] = alldata[i][1]/massone;
      com[i][2] = alldata[i][2]/massone;
    }
  }
}

/* ---------------------------------------------------------------------- */
