#include <fstream>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <time.h>

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
    int molecules_no = 20;
    Molecule *mol_arr = read_molecules(filename,molecules_no);

    /* stratify and divide into training and validation arrays */
    std::sort(mol_arr, mol_arr+molecules_no, compare_molecules);
    int train_no = 10;
    Molecule *train_mol[train_no];
    int validate_no = molecules_no - train_no;
    Molecule *validate_mol[validate_no];
    stratify(mol_arr, train_mol, train_no, validate_mol, validate_no);

    /* compute power spectra descriptors */
    std::cout << "computing power spectra descriptors" << std::endl;
    Descriptor *train_desc = create_descriptor_arr(train_no);
    Descriptor *validate_desc = create_descriptor_arr(validate_no);
    bnu::vector<double> energy(train_no);
    for (int mol_idx=0; mol_idx<train_no; mol_idx++)
    {
        molecule2descriptor(train_mol[mol_idx],train_desc[mol_idx]);
        energy(mol_idx) = train_mol[mol_idx]->energy;
    }
    bnu::vector<double> energy_p(validate_no);
    for (int mol_idx=0; mol_idx<validate_no; mol_idx++)
    {
        molecule2descriptor(validate_mol[mol_idx],validate_desc[mol_idx]);
        energy_p(mol_idx) = validate_mol[mol_idx]->energy;
    }
    free_mol_array(mol_arr);

    std::cout << "computing self similarity" << std::endl;   //for efficient normalisation
    /* compute self local similarity array */
    double **LS = create_local_similarity_array(train_desc,train_no);
    double **LS_p = create_local_similarity_array(validate_desc,validate_no);

    /* self structural similarty */
    double *SS = create_structural_similarity_array(train_desc,LS,train_no);
    double *SS_p = create_structural_similarity_array(validate_desc,LS_p,validate_no);

    std::cout << "cross structural similarity:" << std::endl;
    /* cross structural similarity */
    std::cout << "training matrix..." << std::endl;

    clock_t start, end;
    start = clock();
    structural_similarity(train_desc[train_no-1],train_desc[train_no-2],LS[train_no-1],LS[train_no-2]);
    end = clock();
    std::cout << ((double) (end - start)) / CLOCKS_PER_SEC << std::endl;

    bnu::matrix<double> K(train_no,train_no);
    //#pragma omp parallel for schedule(dynamic)
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
    std::cout << "done" << std::endl;
    std::cout << "validation matrix..." << std::endl;
    bnu::matrix<double> L(train_no,validate_no);
    #pragma omp parallel for schedule(static)
    for (int i=0; i<train_no; i++)
    {
        for (int j=0; j<validate_no; j++)
        {
            double l_ij = structural_similarity(train_desc[i],validate_desc[j],LS[i],LS_p[j]);
            if (isnan(l_ij)) std::cout << l_ij << std::endl;
            L(i,j) = l_ij / sqrt(SS[i]*SS_p[j]);
        }
    }
    std::cout << "done" << std::endl;
    /* free self similarity arrays */
    free_ls_arr(LS,train_no);
    free_ls_arr(LS_p,validate_no);
    free_ss_arr(SS);
    free_ss_arr(SS_p);

    /* TRAINING */
    std::cout << "inverting matrix" << std::endl;
    for (int i=0; i<train_no; i++)
        for (int j=0; j<train_no; j++)
            K(i,j) = pow(K(i,j),ZETA);
    bnu::identity_matrix<double> eye(train_no);
    K += LAMBDA*eye;    //apply ridge parameter
    bnu::matrix<double> invK(train_no,train_no);
    invert_matrix(K,invK);
    bnu::vector<double> alpha = prod(invK,energy);

    /* VALIDATION */
    std::cout << "producing predictions" << std::endl;
    bnu::vector<double> f = prod(trans(L),alpha);
    std::cout << std::endl;

    /* stats */
    std::cout << "stats:" << std::endl;
    double mre = 0;
    double mae = 0;
    double rmse = 0;
    for (int p_idx=0; p_idx<validate_no; p_idx++)
    {
        double err = energy_p(p_idx) - f(p_idx);
        mre += std::abs(err/energy_p(p_idx));
        mae += std::abs(err);
        rmse += pow(err,2);
    }
    mre /= validate_no;
    mae /= validate_no;
    rmse /= validate_no;
    rmse = sqrt(rmse);

    std::cout << "MRE: " << mre*100 << "%" << std::endl;
    std::cout << "MAE: " << mae << std::endl;
    std::cout << "RMSE: " << rmse << std::endl;

    return 0;
}
