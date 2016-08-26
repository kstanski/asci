#include <iostream>
#include <math.h>
#include <boost/numeric/ublas/matrix.hpp>

#include "solver.h"
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
#if VERBOSE
    std::cout << "computing self similarity" << std::endl;
#endif // VERBOSE
    /* compute self local similarity array */   //for efficient normalisation
    double **LS = create_local_similarity_array(train_desc,train_no,diag);
    double **LS_p = create_local_similarity_array(validate_desc,validate_no,diag);

    /* self structural similarty */
    double *SS = create_structural_similarity_array(train_desc,LS,train_no,diag);
    double *SS_p = create_structural_similarity_array(validate_desc,LS_p,validate_no,diag);
#if VERBOSE
    std::cout << "cross structural similarity:" << std::endl;
    std::cout << "training matrix..." << std::endl;
#endif // VERBOSE
    /* cross structural similarity */
    bnu::matrix<double> K(train_no,train_no);
    #pragma omp parallel for schedule(dynamic)
    for (int i=0; i<train_no; i++)
    {
        K(i,i) = 1.0;
        for (int j=i+1; j<train_no; j++)
        {
            double k_ij = structural_similarity(train_desc[i],train_desc[j],LS[i],LS[j],diag);
            k_ij /= sqrt(SS[i]*SS[j]);
            k_ij = pow(k_ij,params.zeta);
            K(i,j) = k_ij;
            K(j,i) = k_ij;
        }
    }
#if VERBOSE
    std::cout << "done" << std::endl;
    std::cout << "validation matrix..." << std::endl;
#endif // VERBOSE
    bnu::matrix<double> L(train_no,validate_no);
    #pragma omp parallel for schedule(dynamic)
    for (int i=0; i<train_no; i++)
    {
        for (int j=0; j<validate_no; j++)
        {
            double l_ij = structural_similarity(train_desc[i],validate_desc[j],LS[i],LS_p[j],diag);
            l_ij /= sqrt(SS[i]*SS_p[j]);
            L(i,j) = pow(l_ij,params.zeta);
        }
    }
#if VERBOSE
    std::cout << "done" << std::endl;
#endif // VERBOSE
    /* free unused arrays */
    free_ls_arr(LS,train_no);
    free_ls_arr(LS_p,validate_no);
    free_ss_arr(SS);
    free_ss_arr(SS_p);

    /* TRAINING */
#if VERBOSE
    std::cout << "solving linear system" << std::endl;
#endif // VERBOSE
    bnu::identity_matrix<double> eye(train_no);
    K += params.lamdba*eye;    //apply ridge parameter

    bnu::vector<double> alpha (train_no);
    int res = solve_linear_system(K,alpha,energy);
    if (res != 0) std::cout << "cannot solve the linear system" << std::endl;;

    /* VALIDATION */
#if VERBOSE
    std::cout << "producing predictions" << std::endl;
#endif // VERBOSE
    bnu::vector<double> f = prod(trans(L),alpha);

    /* stats */
#if VERBOSE
    output_plot_data(energy_p,f);
#endif // VERBOSE
    return produce_stats(energy_p,f);
}
