#include <fstream>
#include <string>
#include <iostream>

#include "molecule.h"
#include "neighbourhood.h"

int main()
{
    //char filename[] = "test_data.xyz";
    char filename[] = "dsgdb7ae2.xyz";
    int molecules_no = 16;
    Molecule *mol_arr = read_molecules(filename,molecules_no);
    Neighbourhood *nhoods_arr[molecules_no];

    for (int i=0; i<molecules_no; i++)
    {
        nhoods_arr[i] = molecule2neighbourhoods(&mol_arr[i]);
    }
    free_mol_array(mol_arr);



    for (int i=0; i<molecules_no; i++)
    {
        free_nhoods(nhoods_arr[i]);
    }

    return 0 ;
}
