#ifndef POWER_SPECTRUM_H_INCLUDED
#define POWER_SPECTRUM_H_INCLUDED

#include <complex>
#include <boost/numeric/ublas/vector.hpp>

#include "molecule.h"

#define L_MAX 5
#define N_MAX 5

typedef boost::numeric::ublas::vector<std::complex<double> > Power_spectrum;
Power_spectrum coords2power_spectrum(Position *coords, int coords_no);
int cart2sph(Position *coords, int coords_no, double *phi, double *theta, double *r);
int get_ps_idx(int l, int n1, int n2);
double radial_basis_function(double r,double cutoff,int n,int n_max);

#endif // POWER_SPECTRUM_H_INCLUDED
