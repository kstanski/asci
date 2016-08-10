#ifndef POWER_SPECTRUM_H_INCLUDED
#define POWER_SPECTRUM_H_INCLUDED

#include <complex>
#include <boost/numeric/ublas/vector.hpp>

#define L_MAX 5
#define N_MAX 5

typedef boost::numeric::ublas::vector<std::complex<double> > Power_spectrum;
Power_spectrum coords2power_spectrum(double **coords, int coords_no);
int cart2sph(double **coords, int coords_no, double *phi, double *theta, double *r);
double get_element(Power_spectrum ps, int l, int n1, int n2);
int set_element(Power_spectrum ps,int l, int n1, int n2, double element);
int spherical_harmonics(std::complex<double> **sh, int coords_no, int l, double *theta, double *phi);
double radial_basis_function(double r,double cutoff,int n,int n_max);
double rbf_poly(double r, double cutoff, int alpha);
int free_power_spectrum(Power_spectrum ps);

#endif // POWER_SPECTRUM_H_INCLUDED
