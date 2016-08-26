#include <iostream>
#include <math.h>

#include "setup.h"
#include "run.h"
#include "stats.h"

int main()
{
    char filename[] = "dsgdb7ae2.xyz";
    int molecules_no = 7102;
    int train_no = 70;
    int validate_no = 30;
    Dataset *dset = setup(filename,molecules_no,train_no,validate_no);

    Params params;
    params.lamdba = pow(10,-3);
    params.zeta = 1;
    double diag[] = {1,1,1,1,1};    //H,C,N,O,S
    params.diag = diag;
    Stats s = run(dset,params);
    std::cout << s.mae << std::endl;
    free_dataset(dset);

    return 0;
}
