
#include <stdlib.h>

#include "descriptor.h"
#include "neighbourhood.h"
#include "power_spectrum.h"

int molecule2descriptor(Molecule *mol_ptr, Descriptor *desc_ptr)
{
    Neighbourhood *nhoods = molecule2neighbourhoods(mol_ptr);
    //Descriptor *desc = (Descriptor *) malloc(MAX_TOTAL*ATOM_TYPES*sizeof(Power_spectrum));
    int coords_no;
    double **coords;
    coords = (double **) malloc(MAX_ATOMS * sizeof(double *));
    for (int a=0; a<MAX_ATOMS; a++)
    {
        coords[a] = (double *) malloc(DIMENSIONS * sizeof(double));
    }

    for (int atom_idx=0; atom_idx<MAX_TOTAL; atom_idx++)
    {
        for (int type=0; type<ATOM_TYPES; type++)
        {
            coords_no = nhoods[atom_idx].last_atom_idx[type]+1;
            for (int i=0; i<coords_no; i++)
            {
                for (int d=0; d<DIMENSIONS; d++)
                {
                    coords[i][d] = nhoods[atom_idx].coords[type][i][d];
                }
            }
            Power_spectrum temp = coords2power_spectrum(coords,coords_no);
            (*desc_ptr)[atom_idx][type] = temp;
        }
    }

    return 0;//desc;
}
