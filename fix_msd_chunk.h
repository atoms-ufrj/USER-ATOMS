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

FixStyle(msd/chunk,FixMSDChunk)

#else

#ifndef LMP_FIX_MSD_CHUNK_H
#define LMP_FIX_MSD_CHUNK_H

#include <stdio.h>
#include "fix.h"
#include "compute.h"

namespace LAMMPS_NS {

class FixMSDChunk : public Fix {
 public:
  FixMSDChunk(class LAMMPS *, int, char **);
  ~FixMSDChunk();
  int setmask();
  void init();
  void end_of_step();
  double compute_array(int, int);

 private:
  class ComputeChunkAtom *cchunk;

  int nfreq;
  char *filename;

  int maxlevel;
  int blocksize;
  int nchunks;
  int step;

  double **com;
  int nlevels;
  int *samples;
  int **counter;
  double ****comtab;
  double ***disptab;

  double **msd;
  int uptodate;

  int ipow(int,int);
  void allocate(int);
  void deallocate();
  void compute_centers_of_mass();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot open fix print file %s

The output file generated by the fix print command cannot be opened

*/
