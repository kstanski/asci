#include "neighbourhood.h"
#include "molecule.h"

Neighbourhood *molecule2neighbourhoods(Molecule *mol_ptr)
{
    /* create neighbourhood array */
    Neighbourhood *nhood_arr = create_nhoods(MAX_TOTAL);
    int atoms_no = mol_ptr->atoms_no;
    int max_types[] = {MAX_H,MAX_C,MAX_N,MAX_O,MAX_S};

    /* set indices of empty arrays to -1 */
    for (int idx=0; idx<MAX_TOTAL; idx++)
    {
        for (int type=0 ; type<ATOM_TYPES; type++)
        {
            nhood_arr[idx].last_atom_idx[type] = -1;
        }
    }

    /* search for close neighbours */
    Position c_i;
    for (int i=0; i<atoms_no; i++)  // for center atom
    {
        c_i = mol_ptr->ff_coords[i];
        for (int j=i; j<atoms_no; j++)  // for neighbouring atom
        {
            Position c_j;
            c_j = mol_ptr->ff_coords[j];

            if (bg::distance(c_i,c_j) < CUTOFF)
            {
                /* add neighbour of i */
                int neighbour_type = mol_ptr->atom_types[j];
                nhood_arr[i].last_atom_idx[neighbour_type]++;
                int last = nhood_arr[i].last_atom_idx[neighbour_type];
                Position diff = pos_diff(c_j,c_i);
                nhood_arr[i].coords[neighbour_type][last] = diff;

                if (i != j)
                {
                    /* add neighbour of j */
                    neighbour_type = mol_ptr->atom_types[i];
                    nhood_arr[j].last_atom_idx[neighbour_type]++;
                    last = nhood_arr[j].last_atom_idx[neighbour_type];
                    nhood_arr[j].coords[neighbour_type][last] = diff;
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
            Position dummy(0,0,0);
            nhood_arr[last].coords[type][0] = dummy;
            lack--;
        }
    }

    return nhood_arr;
}

Position pos_diff(Position a, Position b)
{
    Position diff;
    diff.set<0>(a.get<0>() - b.get<0>());
    diff.set<1>(a.get<1>() - b.get<1>());
    diff.set<2>(a.get<2>() - b.get<2>());
    return diff;
}

Neighbourhood *create_nhoods(int nhood_no)
{
    Neighbourhood *nhood_arr;
    nhood_arr = (Neighbourhood *) malloc(nhood_no*sizeof(Neighbourhood));
    for (int nhood_idx=0; nhood_idx<nhood_no; nhood_idx++)
    {
        for (int type=0; type<ATOM_TYPES; type++)
        {
            nhood_arr[nhood_idx].coords[type] = (Position *) malloc(MAX_ATOMS*sizeof(Position));
        }
    }
    if (nhood_arr == NULL)
    {
        fprintf(stderr,"Not enough memory: neighbourhoods");
        exit (1);
    }
    return nhood_arr;
}

int free_nhoods(Neighbourhood *nhoods, int nhood_no)
{
    for (int nhood_idx=0; nhood_idx<nhood_no; nhood_idx++)
    {
        for (int type=0; type<ATOM_TYPES; type++)
        {
            free(nhoods[nhood_idx].coords[type]);
        }
    }
    free(nhoods);
    return 0;
}
