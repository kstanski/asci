#ifndef RUN_H_INCLUDED
#define RUN_H_INCLUDED

#include "setup.h"
#include "molecule.h"

typedef struct params
{
    double lamdba;
    double zeta;
    double diag[ATOM_TYPES];
} Params;

Stats run(Dataset *dset, Params params);

#endif // RUN_H_INCLUDED
