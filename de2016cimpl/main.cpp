#include <stdlib.h>
#include <iostream>
#include <math.h>

#include "stats.h"
#include "setup.h"
#include "run.h"

int main()
{
    char filename[] = "dsgdb7ae2.xyz";
    int molecules_no = 7102;
    int train_no = 100;
    int validate_no = 70;
    Dataset *dset = setup(filename,molecules_no,train_no,validate_no);

    Params params;
    params.lamdba = pow(10,-3);
    params.zeta = 1;
    double diag[] = {1,2,2,3,1};    //H,C,N,O,S
    for (int type=0; type<ATOM_TYPES; type++)
        params.diag[type] = diag[type];
    Stats s = run(dset,params);

    return 0;
}
