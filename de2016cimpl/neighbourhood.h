#ifndef NEIGHBOURHOOD_H_INCLUDED
#define NEIGHBOURHOOD_H_INCLUDED

#include "molecule.h"

#define MAX_H 16
#define MAX_C 7
#define MAX_N 3
#define MAX_O 3
#define MAX_S 1
#define MAX_TOTAL (MAX_H+MAX_C+MAX_N+MAX_O+MAX_S)
#define CUTOFF 3

typedef struct neighbourhood
{
    Position *coords[ATOM_TYPES];    // 5 atom species
    int last_atom_idx[ATOM_TYPES];
} Neighbourhood;

Neighbourhood *molecule2neighbourhoods(Molecule *mol_ptr);
Position pos_diff(Position a, Position b);
Neighbourhood *create_nhoods(int nhood_no);
int free_nhoods(Neighbourhood *nhoods, int nhood_no);

#endif // NEIGHBOURHOOD_H_INCLUDED
