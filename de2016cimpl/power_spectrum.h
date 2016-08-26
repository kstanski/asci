#ifndef POWER_SPECTRUM_H_INCLUDED
#define POWER_SPECTRUM_H_INCLUDED

#include "molecule.h"

#define L_MAX 12
#define N_MAX 10
#define PS_LEN (L_MAX+1)*N_MAX*N_MAX

typedef float ps_element_type;
typedef ps_element_type *Power_spectrum;
Power_spectrum coords2power_spectrum(Position *coords, int coords_no);
int cart2sph(Position *coords, int coords_no, double *phi, double *theta, double *r);
ps_element_type sh_real_form(int l, int m, double theta, double phi);
int get_ps_idx(int l, int n1, int n2);
ps_element_type radial_basis_function(double r,double cutoff,int n,int n_max);
double dot_prod(Power_spectrum A, Power_spectrum B);
int normalise(Power_spectrum ps);

#endif // POWER_SPECTRUM_H_INCLUDED
