#ifndef DESCRIPTOR_H_INCLUDED
#define DESCRIPTOR_H_INCLUDED

#include "molecule.h"
#include "power_spectrum.h"

typedef Power_spectrum **Descriptor;

int molecule2descriptor(Molecule *mol_ptr, Descriptor *desc_ptr);

#endif // DESCRIPTOR_H_INCLUDED
