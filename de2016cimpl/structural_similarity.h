#ifndef STRUCTURAL_SIMILARITY_H_INCLUDED
#define STRUCTURAL_SIMILARITY_H_INCLUDED

#define REG_PARAM 0.5
#define ITERATIONS_NO 20

#include "power_spectrum.h"

double structural_similarity(Power_spectrum **A, Power_spectrum **B, double *LSA, double *LSB);

#endif // STRUCTURAL_SIMILARITY_H_INCLUDED
