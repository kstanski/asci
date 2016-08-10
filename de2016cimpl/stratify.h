#ifndef STRATIFY_H_INCLUDED
#define STRATIFY_H_INCLUDED

#include "molecule.h"

int stratify(Molecule *sorted_arr, Molecule **train_arr, int train_no, Molecule **validate_arr, int validate_no);

#endif // STRATIFY_H_INCLUDED
