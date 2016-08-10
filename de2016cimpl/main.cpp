#include <fstream>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <math.h>

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "molecule.h"
#include "neighbourhood.h"
#include "power_spectrum.h"
#include "invert_matrix.h"
#include "stratify.h"
#include "descriptor.h"
#include "local_similarity.h"
#include "structural_similarity.h"

int main()
{
    namespace bnu = boost::numeric::ublas;
    /* read molecules from file */
    //char filename[] = "test_data.xyz";
    char filename[] = "dsgdb7ae2.xyz";
    int molecules_no = 100;
    Molecule *mol_arr = read_molecules(filename,molecules_no);

    /* stratify and divide into training and validation */
    std::sort(mol_arr, mol_arr+molecules_no, compare_molecules);
    int train_no = 20;
    Molecule *train_mol[train_no];
    int validate_no = molecules_no - train_no;
    Molecule *validate_mol[validate_no];
    stratify(mol_arr, train_mol, train_no, validate_mol, validate_no);


    /* TRAINING */


    /* compute power spectra descriptors */
    Descriptor train_desc[train_no];
    for (int desc_idx=0; desc_idx<train_no; desc_idx++)
    {
        train_desc[desc_idx] = (Descriptor) malloc(MAX_TOTAL*sizeof(Power_spectrum *));
        for (int atom_idx=0; atom_idx<MAX_TOTAL; atom_idx++)
        {
            train_desc[desc_idx][atom_idx] = (Power_spectrum *) malloc(ATOM_TYPES*sizeof(Power_spectrum));
        }
    }

    bnu::vector<double> energy(train_no);
    for (int mol_idx=0; mol_idx<train_no; mol_idx++)
    {
        molecule2descriptor(train_mol[mol_idx],&(train_desc[mol_idx]));
        energy(mol_idx) = train_mol[mol_idx]->energy;
    }

    /* compute local similarity */
    double *LS[train_no];
    for (int mol_idx=0; mol_idx<train_no; mol_idx++)
    {
        LS[mol_idx] = (double *) malloc(MAX_TOTAL*sizeof(double));
        for (int atom_idx=0; atom_idx<MAX_TOTAL; atom_idx++)
        {
            double temp = local_similarity(train_desc[mol_idx][atom_idx],train_desc[mol_idx][atom_idx]);
            //std::cout << temp;
            LS[mol_idx][atom_idx] = temp;
        }
        //std::cout << std::endl;
    }

    for (int idx=0; idx<train_no; idx++)
    {
        for (int atom2=0; atom2<MAX_TOTAL; atom2++)
        {
            //std::cout << LS[idx][atom2];
        }
        //std::cout << std::endl;
    }

    /* self structural similarty */
    double SS[train_no];
    for (int mol_idx=0; mol_idx<train_no; mol_idx++)
    {
        Descriptor mol_desc = train_desc[mol_idx];
        double *LSA = LS[mol_idx];
        SS[train_no] = structural_similarity(mol_desc,mol_desc,LSA,LSA);
    }

    for (int i=0; i<train_no; i++) std::cout << SS[i] << std::endl;

    /* cross structural similarity */
    bnu::matrix<double> K(train_no,train_no), invK(train_no,train_no);
    for (int i=0; i<train_no; i++)
    {
        K(i,i) = 1;
        for (int j=i+1; j<train_no; j++)
        {
            double k_ij = structural_similarity(train_desc[i],train_desc[j],LS[i],LS[j]);
            k_ij /= sqrt(SS[i]*SS[j]);
            K(i,j) = k_ij;
            K(j,i) = k_ij;
        }
    }

    for (int idx=0; idx<train_no; idx++)
    {
        for (int idx2=0; idx2<train_no; idx2++)
        {
            std::cout << K(idx,idx2);
        }
        std::cout << std::endl;
    }

    invert_matrix(K,invK);
    bnu::vector<double> alpha = prod(invK,energy);





    /* VALIDATION */




    return 0;
}
