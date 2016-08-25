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
#include "solver.h"
#include "stratify.h"
#include "descriptor.h"
#include "local_similarity.h"
#include "structural_similarity.h"
#include "stats.h"
#include "setup.h"

#include "run.h"

Stats run(Dataset *dset, Params params)
{
    /* read dataset */
    int train_no = dset->train_no;
    int validate_no = dset->validate_no;
    Descriptor *train_desc = dset->train_desc;
    Descriptor *validate_desc = dset->validate_desc;
    bnu::vector<double> energy(train_no);
    for (int i=0; i<train_no; i++) energy(i) = dset->train_val[i];
    bnu::vector<double> energy_p(validate_no);
    for (int i=0; i<validate_no; i++) energy_p(i) = dset->validate_val[i];
    double *diag = params.diag;

    std::cout << "computing self similarity" << std::endl;   //for efficient normalisation
    /* compute self local similarity array */
    double **LS = create_local_similarity_array(train_desc,train_no,diag);
    double **LS_p = create_local_similarity_array(validate_desc,validate_no,diag);

    /* self structural similarty */
    double *SS = create_structural_similarity_array(train_desc,LS,train_no,diag);
    double *SS_p = create_structural_similarity_array(validate_desc,LS_p,validate_no,diag);

    std::cout << "cross structural similarity:" << std::endl;
    /* cross structural similarity */
    std::cout << "training matrix..." << std::endl;

    bnu::matrix<double> K(train_no,train_no);
    #pragma omp parallel for schedule(dynamic)
    for (int i=0; i<train_no; i++)
    {
        K(i,i) = 1.0;
        for (int j=i+1; j<train_no; j++)
        {
            double k_ij = structural_similarity(train_desc[i],train_desc[j],LS[i],LS[j],diag);
            k_ij /= sqrt(SS[i]*SS[j]);
            K(i,j) = k_ij;
            K(j,i) = k_ij;
        }
    }
    std::cout << "done" << std::endl;
    std::cout << "validation matrix..." << std::endl;
    bnu::matrix<double> L(train_no,validate_no);
    #pragma omp parallel for schedule(dynamic)
    for (int i=0; i<train_no; i++)
    {
        for (int j=0; j<validate_no; j++)
        {
            double l_ij = structural_similarity(train_desc[i],validate_desc[j],LS[i],LS_p[j],diag);
            L(i,j) = l_ij / sqrt(SS[i]*SS_p[j]);
        }
    }
    std::cout << "done" << std::endl;
    /* free unused arrays */
    free_desc_arr(train_desc,train_no);
    free_desc_arr(validate_desc,validate_no);
    free_ls_arr(LS,train_no);
    free_ls_arr(LS_p,validate_no);
    free_ss_arr(SS);
    free_ss_arr(SS_p);

    /* TRAINING */
    std::cout << "solving linear system" << std::endl;
    for (int i=0; i<train_no; i++)
        for (int j=0; j<train_no; j++)
            K(i,j) = pow(K(i,j),params.zeta);
    bnu::identity_matrix<double> eye(train_no);
    K += params.lamdba*eye;    //apply ridge parameter

    bnu::vector<double> alpha (train_no);
    int res = solve_linear_system(K,alpha,energy);
    if (res != 0) std::cout << "cannot solve the linear system" << std::endl;;

    /* VALIDATION */
    std::cout << "producing predictions" << std::endl;
    bnu::vector<double> f = prod(trans(L),alpha);
    std::cout << std::endl;

    /* stats */
    std::cout << "stats:" << std::endl;

    output_plot_data(energy_p,f);

    return produce_stats(energy_p,f);
}
