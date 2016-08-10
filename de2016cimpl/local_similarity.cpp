#include <math.h>
#include <stdlib.h>

#include "local_similarity.h"
#include "molecule.h"

double local_similarity(Power_spectrum *S1, Power_spectrum *S2)
{
    int ntypes = ATOM_TYPES;
    double similarity_kernel[ntypes][ntypes];

    /* unity matrix */
    for (int i=0; i<ntypes; i++)
    {
        for (int j=0; j<ntypes; j++)
        {
            if (i == j) similarity_kernel[i][j] = 1;
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
                double dot = std::real(inner_prod(S1[i],S2[j]));
                sum += similarity_kernel[i][j] * pow(dot,LOCAL_ZETA);
            }
        }
    }
    //std::cout << sum << std::endl;
    return sum;
}
