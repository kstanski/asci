#ifndef SOAP_C_WRAP_H_INCLUDED
#define SOAP_C_WRAP_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

typedef struct dset Dset;

Dset *setup_soap(const char *filename, int molecules_no, int train_no, int validate_no);
double run_soap(Dset *dset, double *params);
int free_dset(Dset *dset);

#ifdef __cplusplus
}
#endif

#endif // SOAP_C_WRAP_H_INCLUDED
