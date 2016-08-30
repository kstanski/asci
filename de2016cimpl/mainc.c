#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "soap_c_wrap.h"

int main()
{
    char filename[] = "dsgdb7ae2.xyz";
    int molecules_no = 7102;
    int train_no = 50;
    int validate_no = 20;
    void *dset = setup_soap(filename,molecules_no,train_no,validate_no);

    double *params = (double *) malloc(7*sizeof(double));
    double diag[] = {1,1,1,1,1};    //H,C,N,O,S
    int i;
    for (i=0; i<5; i++) params[i] = diag[i];
    params[5] = pow(10,-3); //lambda
    params[6] = 1;          //zeta
    double mae = run_soap(dset,params);
    free(params);
    printf("MAE: %f\n", mae);
    free_dset(dset);

    return 0;
}
