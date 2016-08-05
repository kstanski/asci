#ifndef POWER_SPECTRUM_H_INCLUDED
#define POWER_SPECTRUM_H_INCLUDED

#define L_MAX 12
#define N_MAX 10

typedef double *Power_spectrum;
Power_spectrum coords2power_spectrum(double **coords);
int free_power_spectrum(Power_spectrum ps);

#endif // POWER_SPECTRUM_H_INCLUDED
