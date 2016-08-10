
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
            double ls = local_similarity(A[i],B[j])/sqrt(LSA[i]*LSB[j]);
            C(i,j) = ls;
            Sink(i,j) = exp((ls-1)/REG_PARAM);
        }
    }

    /* Sinkhorn algorithm */
    bnu::vector<double> v(MAX_TOTAL), u(MAX_TOTAL), en(MAX_TOTAL), temp(MAX_TOTAL);
    for (int i=0; i<MAX_TOTAL; i++) en(i) = 1/MAX_TOTAL;

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
            Sink(i,j) *= u(i)*v(j);   //not sure about the order of i,j and u,v in this line
        }
    }

    Sink = prod(trans(Sink),C);
    double trace = 0;
    for (int i=0; i<MAX_TOTAL; i++) trace += Sink(i,i);

    return trace;
}
