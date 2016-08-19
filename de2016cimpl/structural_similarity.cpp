
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <math.h>

#include "structural_similarity.h"
#include "descriptor.h"
#include "neighbourhood.h"
#include "local_similarity.h"

double structural_similarity(Descriptor A, Descriptor B, double *LSA, double *LSB)
{
    namespace bnu = boost::numeric::ublas;
    bnu::matrix<double> C(MAX_TOTAL, MAX_TOTAL);
    bnu::matrix<double> Sink(MAX_TOTAL, MAX_TOTAL);

    for (int i=0; i<MAX_TOTAL; i++)
    {
        for (int j=0; j<MAX_TOTAL; j++)
        {
            double ls = local_similarity(A[i],B[j])/sqrt(LSA[i]*LSB[j]);    //normalised
            C(i,j) = ls;
            Sink(i,j) = exp((ls-1)/REG_PARAM);
        }
    }

    /* Sinkhorn algorithm */
    bnu::vector<double> v(MAX_TOTAL), u(MAX_TOTAL), en(MAX_TOTAL), temp(MAX_TOTAL);
    for (int i=0; i<MAX_TOTAL; i++) en(i) = 1.0/MAX_TOTAL;
    v = en;
    for (int counter=0; counter<ITERATIONS_NO; counter++)
    {
        temp = prod(Sink,v);
        for (int i=0; i<MAX_TOTAL; i++) u(i) = en(i)/temp(i);
        temp = prod(trans(Sink),u);
        for (int i=0; i<MAX_TOTAL; i++) v(i) = en(i)/temp(i);
    }

    for (int i=0; i<MAX_TOTAL; i++)
    {
        for (int j=0; j<MAX_TOTAL; j++)
        {
            Sink(i,j) *= u(i)*v(j);   //not sure about the order of i,j and u,v in this line (there is a typo in the original paper)
        }
    }

    Sink = prod(trans(Sink),C);
    double trace = 0;
    for (int i=0; i<MAX_TOTAL; i++) trace += Sink(i,i);

    return trace;
}

double *create_structural_similarity_array(Descriptor *desc_arr, double **ls_arr, int desc_no)
{
    double *ss_arr = (double *) malloc(desc_no*sizeof(double));
    #pragma omp parallel for schedule(dynamic)
    for (int mol_idx=0; mol_idx<desc_no; mol_idx++)
    {
        Descriptor mol_desc = desc_arr[mol_idx];
        double *LSA = ls_arr[mol_idx];
        ss_arr[mol_idx] = structural_similarity(mol_desc,mol_desc,LSA,LSA);
    }
    return ss_arr;
}

int free_ss_arr(double *ss_arr)
{
    free(ss_arr);
    return 0;
}
