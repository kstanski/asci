#ifndef DESCRIPTOR_H_INCLUDED
#define DESCRIPTOR_H_INCLUDED

#include "molecule.h"
#include "power_spectrum.h"

typedef Power_spectrum **Descriptor;

int molecule2descriptor(Molecule *mol_ptr, Descriptor desc);
Descriptor *create_descriptor_arr(int desc_no);
int free_desc_arr(Descriptor *desc_arr, int desc_no);

#endif // DESCRIPTOR_H_INCLUDED
