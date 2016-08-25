#ifndef STRUCTURAL_SIMILARITY_H_INCLUDED
#define STRUCTURAL_SIMILARITY_H_INCLUDED

#define REG_PARAM 0.5
#define ITERATIONS_NO 20

#include "power_spectrum.h"
#include "descriptor.h"

double structural_similarity(Power_spectrum **A, Power_spectrum **B, double *LSA, double *LSB, double *diag);
double *create_structural_similarity_array(Descriptor *desc_arr, double **ls_arr, int desc_no, double *diag);
int free_ss_arr(double *ss_arr);

#endif // STRUCTURAL_SIMILARITY_H_INCLUDED
