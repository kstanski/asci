#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "neighbourhood.h"

Neighbourhood *molecule2neighbourhoods(Molecule *mol_ptr)
{
    /* create neighbourhood array */
    Neighbourhood *nhood_arr;
    int atoms_no = mol_ptr->atoms_no;
    int max_total = MAX_H+MAX_C+MAX_N+MAX_O+MAX_S;
    int max_types[] = {MAX_H,MAX_C,MAX_N,MAX_O,MAX_S};
    nhood_arr = (Neighbourhood *) malloc(max_total*sizeof(Neighbourhood));
    if (nhood_arr == NULL)
    {
        fprintf(stderr,"Not enough memory: neighbourhoods");
        exit (1);
    }

    /* set indices of empty arrays to -1 */
    for (int idx=0; idx<max_total; idx++)
    {
        for (int type=0 ; type<ATOM_TYPES; type++)
        {
            nhood_arr[idx].last_atom_idx[type] = -1;
        }
    }

    /* search for close neighbours */
    double c_i[DIMENSIONS];
    for (int i=0; i<atoms_no; i++)  // for center atom
    {
        for (int idx=0;idx<DIMENSIONS;idx++) c_i[idx] = mol_ptr->ff_coords[i][idx];
        for (int j=i; j<atoms_no; j++)  // for neighbouring atom
        {
            double norm = 0;
            double c_j[DIMENSIONS];
            for (int idx=0;idx<DIMENSIONS;idx++)
            {
                c_j[idx] = mol_ptr->ff_coords[j][idx];
                norm += pow((c_i[idx] - c_j[idx]),2);
            }
            norm = sqrt(norm);

            if (norm < CUTOFF)
            {
                /* add neighbour of i */
                int neighbour_type = mol_ptr->atom_types[j];
                nhood_arr[i].last_atom_idx[neighbour_type]++;
                int last = nhood_arr[i].last_atom_idx[neighbour_type];
                for (int idx=0;idx<DIMENSIONS;idx++)
                {
                    nhood_arr[i].coords[neighbour_type][last][idx] = c_j[idx] - c_i[idx];
                }

                if (i != j)
                {
                    /* add neighbour of j */
                    neighbour_type = mol_ptr->atom_types[i];
                    nhood_arr[j].last_atom_idx[neighbour_type]++;
                    last = nhood_arr[j].last_atom_idx[neighbour_type];
                    for (int idx=0;idx<DIMENSIONS;idx++)
                    {
                        nhood_arr[j].coords[neighbour_type][last][idx] = c_i[idx] - c_j[idx];
                    }
                }
            }
        }
    }

    /* top up with dummy atoms */
    int last = atoms_no - 1;
    for (int type=0; type<ATOM_TYPES; type++)
    {
        int lack = max_types[type] - (mol_ptr->types_total)[type];
        while (lack > 0)
        {
            last++;
            nhood_arr[last].last_atom_idx[type] = 0;
            for (int idx=0;idx<DIMENSIONS;idx++)
                nhood_arr[last].coords[type][0][idx] = 0;
            lack--;
        }
    }

    return nhood_arr;
}

int free_nhoods(Neighbourhood *nhoods)
{
    int len = sizeof(nhoods)/sizeof(Neighbourhood);
    for (int idx=0; idx<len; idx++)
    {
        free(&nhoods[idx]);
    }
    free(nhoods);
    return 0;
}
