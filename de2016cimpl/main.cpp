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

#define LAMBDA pow(10,-12)
#define ZETA 1

int main()
{
    namespace bnu = boost::numeric::ublas;
    /* read molecules from file */
    //char filename[] = "test_data.xyz";
    char filename[] = "dsgdb7ae2.xyz";
    int molecules_no = 40;
    Molecule *mol_arr = read_molecules(filename,molecules_no);

    /* stratify and divide into training and validation arrays */
    std::sort(mol_arr, mol_arr+molecules_no, compare_molecules);
    int train_no = 20;
    Molecule *train_mol[train_no];
    int validate_no = molecules_no - train_no;
    Molecule *validate_mol[validate_no];
    stratify(mol_arr, train_mol, train_no, validate_mol, validate_no);


    /* TRAINING */


    /* compute power spectra descriptors */
    Descriptor *train_desc = create_descriptor_arr(train_no);
    bnu::vector<double> energy(train_no);
    for (int mol_idx=0; mol_idx<train_no; mol_idx++)
    {
        molecule2descriptor(train_mol[mol_idx],train_desc[mol_idx]);
        energy(mol_idx) = train_mol[mol_idx]->energy;
    }

    /* compute self local similarity array */
    double **LS = create_local_similarity_array(train_desc,train_no);
    for (int i=0; i<train_no; i++)
    {
        for (int j=0; j<MAX_TOTAL; j++)
            std::cout << LS[i][j];
        std::cout << std::endl;
    }

    /* self structural similarty */
    double *SS = create_structural_similarity_array(train_desc,LS,train_no);
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

    std::cout << std::endl;
    for (int idx=0; idx<train_no; idx++)
    {
        for (int idx2=0; idx2<train_no; idx2++)
            std::cout << K(idx,idx2) << " ";
        std::cout << std::endl;
    }

    bnu::identity_matrix<double> eye(train_no);
    K += LAMBDA*eye;
    invert_matrix(K,invK);
    bnu::vector<double> alpha = prod(invK,energy);


    /* VALIDATION */


    /* compute power spectra descriptors */
    Descriptor *validate_desc = create_descriptor_arr(validate_no);
    bnu::vector<double> energy_p(validate_no);
    for (int mol_idx=0; mol_idx<validate_no; mol_idx++)
    {
        molecule2descriptor(validate_mol[mol_idx],validate_desc[mol_idx]);
        energy_p(mol_idx) = validate_mol[mol_idx]->energy;
    }

    /* compute self local similarity array */
    double **LS_p = create_local_similarity_array(validate_desc,validate_no);
    for (int i=0; i<validate_no; i++)
    {
        for (int j=0; j<MAX_TOTAL; j++)
            std::cout << LS_p[i][j];
        std::cout << std::endl;
    }

    /* self structural similarty */
    double *SS_p = create_structural_similarity_array(validate_desc,LS_p,validate_no);
    for (int i=0; i<validate_no; i++) std::cout << SS_p[i] << std::endl;

    /* cross structural similarity */
    bnu::matrix<double> L(train_no,validate_no);
    for (int i=0; i<train_no; i++)
    {
        for (int j=0; j<validate_no; j++)
        {
            double l_ij = structural_similarity(train_desc[i],validate_desc[j],LS[i],LS_p[j]);
            L(i,j) = l_ij / sqrt(SS[i]*SS_p[j]);
        }
    }

    bnu::vector<double> f = prod(trans(L),alpha);

    double mae = 0;
    for (int p_idx=0; p_idx<validate_no; p_idx++)
    {
        mae += std::abs(energy_p(p_idx) - f(p_idx));
    }
    mae /= validate_no;

    std::cout << "MAE: " << mae << std::endl;

    return 0;
}
