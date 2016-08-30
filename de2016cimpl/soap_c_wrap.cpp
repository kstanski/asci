#include "setup.h"
#include "run.h"
#include "stats.h"
#include "molecule.h"
#include "soap_c_wrap.h"

extern "C" {
    void *setup_soap(const char *filename, int molecules_no, int train_no, int validate_no)
    {
        return setup(filename,molecules_no,train_no,validate_no);
    }

    double run_soap(void *dset, double *params)
    {
        Params p;
        double diag[ATOM_TYPES];    //H,C,N,O,S
        for (int i=0; i<ATOM_TYPES; i++) diag[i] = params[i];
        p.diag = diag;
        p.lamdba = params[ATOM_TYPES];
        p.zeta = params[ATOM_TYPES+1];
        Stats s = run((Dataset *)dset,p);
        return s.mae;
    }

    int free_dset(void *dset)
    {
        free_dataset((Dataset *)dset);
        return 0;
    }
}
