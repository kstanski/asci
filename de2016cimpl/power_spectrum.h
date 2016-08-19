#ifndef POWER_SPECTRUM_H_INCLUDED
#define POWER_SPECTRUM_H_INCLUDED

#include <complex>
#include <boost/numeric/ublas/vector.hpp>

#include "molecule.h"

#define L_MAX 12
#define N_MAX 10

typedef float ps_element_type;
typedef boost::numeric::ublas::vector<ps_element_type> Power_spectrum;
Power_spectrum coords2power_spectrum(Position *coords, int coords_no);
int cart2sph(Position *coords, int coords_no, double *phi, double *theta, double *r);
ps_element_type sh_real_form(int l, int m, double theta, double phi);
int get_ps_idx(int l, int n1, int n2);
ps_element_type radial_basis_function(double r,double cutoff,int n,int n_max);
double dot_prod(Power_spectrum *A, Power_spectrum *B);

#endif // POWER_SPECTRUM_H_INCLUDED
