
#include <stdlib.h>

#include "descriptor.h"
#include "neighbourhood.h"
#include "power_spectrum.h"

int molecule2descriptor(Molecule *mol_ptr, Descriptor desc)
{
    Neighbourhood *nhoods = molecule2neighbourhoods(mol_ptr);
    for (int atom_idx=0; atom_idx<MAX_TOTAL; atom_idx++)
    {
        for (int type=0; type<ATOM_TYPES; type++)
        {
            int coords_no = nhoods[atom_idx].last_atom_idx[type]+1;
            Position *coords = nhoods[atom_idx].coords[type];
            desc[atom_idx][type] = coords2power_spectrum(coords,coords_no);
        }
    }
    free_nhoods(nhoods,MAX_TOTAL);
    return 0;
}

Descriptor *create_descriptor_arr(int desc_no)
{
    Descriptor *desc_arr = (Descriptor *) malloc(desc_no*sizeof(Descriptor));
    for (int desc_idx=0; desc_idx<desc_no; desc_idx++)
    {
        desc_arr[desc_idx] = (Descriptor) malloc(MAX_TOTAL*sizeof(Power_spectrum *));
        for (int atom_idx=0; atom_idx<MAX_TOTAL; atom_idx++)
        {
            desc_arr[desc_idx][atom_idx] = (Power_spectrum *) malloc(ATOM_TYPES*sizeof(Power_spectrum));
        }
    }
    return desc_arr;
}

int free_desc_arr(Descriptor *desc_arr, int desc_no)
{
    for (int desc_idx=0; desc_idx<desc_no; desc_idx++)
    {
        for (int atom_idx=0; atom_idx<MAX_TOTAL; atom_idx++)
        {
            for (int type=0; type<ATOM_TYPES; type++)
                free(desc_arr[desc_idx][atom_idx][type]);
            free(desc_arr[desc_idx][atom_idx]);
        }
        free(desc_arr[desc_idx]);
    }
    free(desc_arr);
    return 0;
}
