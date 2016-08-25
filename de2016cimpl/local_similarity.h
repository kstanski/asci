#ifndef LOCAL_SIMILARITY_H_INCLUDED
#define LOCAL_SIMILARITY_H_INCLUDED

#include "power_spectrum.h"
#include "descriptor.h"

#define LOCAL_ZETA 1

double local_similarity(Power_spectrum *S1, Power_spectrum *S2, double *diag);
double **create_local_similarity_array(Descriptor *desc_arr,int desc_no,double *diag);
int free_ls_arr(double **ls_arr, int desc_no);

#endif // LOCAL_SIMILARITY_H_INCLUDED
