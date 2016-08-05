#ifndef NEIGHBOURHOOD_H_INCLUDED
#define NEIGHBOURHOOD_H_INCLUDED

#include "molecule.h"

#define MAX_H 16
#define MAX_C 7
#define MAX_N 3
#define MAX_O 3
#define MAX_S 1
#define CUTOFF 3

typedef struct neighbourhood
{
    double coords[ATOM_TYPES][MAX_ATOMS][DIMENSIONS];    // 5 atom species
    int last_atom_idx[ATOM_TYPES];
} Neighbourhood;

Neighbourhood *molecule2neighbourhoods(Molecule *mol_ptr);
int free_nhoods(Neighbourhood *nhoods);

#endif // NEIGHBOURHOOD_H_INCLUDED
