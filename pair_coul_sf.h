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

#ifdef PAIR_CLASS

PairStyle(coul/sf,PairCoulSF)

#else

#ifndef LMP_PAIR_COUL_SF_H
#define LMP_PAIR_COUL_SF_H

#include "pair.h"

namespace LAMMPS_NS {

class PairCoulSF : public Pair {
 public:
  PairCoulSF(class LAMMPS *);
  ~PairCoulSF();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  virtual void init_style();
  virtual double init_one(int, int);
//  void modify_params(int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);

 protected:
  double cut_coul,cut_coulsq;
  double f_shift,e_shift;
  double e_self;
  int self_flag;

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

E: Pair style lj/cut/coul/dsf requires atom attribute q

The atom style defined does not have these attributes.

*/
