#include <math.h>

#include "local_similarity.h"
#include "neighbourhood.h"

#define LOCAL_ZETA 1

double local_similarity(Power_spectrum *S1, Power_spectrum *S2, double *diag)
{
    int ntypes = ATOM_TYPES;
    double similarity_kernel[ntypes][ntypes];

    /* diagonal matrix */
    for (int i=0; i<ntypes; i++)
    {
        for (int j=0; j<ntypes; j++)
        {
            if (i == j) similarity_kernel[i][i] = diag[i];
            else similarity_kernel[i][j] = 0;
        }
    }

    double sum = 0;
    for (int i=0; i<ntypes; i++)
    {
        for (int j=0; j<ntypes; j++)
        {
            if (similarity_kernel[i][j] != 0)
            {
                double dot = dot_prod(S1[i],S2[j]);
                sum += similarity_kernel[i][j] * pow(dot,LOCAL_ZETA);
            }
        }
    }
    return sum;
}

double **create_local_similarity_array(Descriptor *desc_arr,int desc_no,double *diag)
{
    double **ls_arr = (double **) malloc(desc_no*sizeof(double *));
    #pragma omp parallel for schedule(dynamic)
    for (int mol_idx=0; mol_idx<desc_no; mol_idx++)
    {
        ls_arr[mol_idx] = (double *) malloc(MAX_TOTAL*sizeof(double));
        for (int atom_idx=0; atom_idx<MAX_TOTAL; atom_idx++)
        {
            Power_spectrum *ps_arr = desc_arr[mol_idx][atom_idx];
            ls_arr[mol_idx][atom_idx] = local_similarity(ps_arr,ps_arr,diag);
        }
    }
    return ls_arr;
}

int free_ls_arr(double **ls_arr, int desc_no)
{
    for (int mol_idx=0; mol_idx<desc_no; mol_idx++)
        free(ls_arr[mol_idx]);

    free(ls_arr);
    return 0;
}

