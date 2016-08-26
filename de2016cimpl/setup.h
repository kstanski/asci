#ifndef SETUP_H_INCLUDED
#define SETUP_H_INCLUDED

#define VERBOSE 1

#include "descriptor.h"

typedef struct dataset
{
    int train_no;
    Descriptor *train_desc;
    double *train_val;
    int validate_no;
    Descriptor *validate_desc;
    double *validate_val;
} Dataset;

Dataset *setup(const char *filename, int molecules_no, int train_no, int validate_no);
Dataset *create_dataset(int train_no, int validate_no);
int free_dataset(Dataset *dset);

#endif // SETUP_H_INCLUDED
