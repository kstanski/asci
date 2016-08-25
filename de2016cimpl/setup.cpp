#include <boost/numeric/ublas/vector.hpp>

#include "setup.h"
#include "stratify.h"



#include <fstream>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <time.h>

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "molecule.h"
#include "neighbourhood.h"
#include "power_spectrum.h"
#include "solver.h"
#include "stratify.h"
#include "descriptor.h"
#include "local_similarity.h"
#include "structural_similarity.h"
#include "stats.h"

#define LAMBDA pow(10,-3)
#define ZETA 1



Dataset *setup(const char *filename, int molecules_no, int train_no, int validate_no)
{
    Dataset *dset = create_dataset(train_no,validate_no);

    Molecule **mol_arr = read_molecules(filename,molecules_no);

    /* stratify and divide into training and validation arrays */
    std::sort(mol_arr, mol_arr+molecules_no, compare_molecules);
    int sample_no = train_no + validate_no;
    int hold_out_no = molecules_no - sample_no;
    Molecule *sample_mol[sample_no];
    Molecule *hold_out[hold_out_no];
    stratify(mol_arr, sample_mol, sample_no, hold_out, hold_out_no);

    Molecule *train_mol[train_no];
    Molecule *validate_mol[validate_no];
    stratify(sample_mol, train_mol, train_no, validate_mol, validate_no);

    /* compute power spectra descriptors */
    std::cout << "computing power spectra descriptors" << std::endl;

    #pragma omp parallel for schedule(dynamic)
    for (int mol_idx=0; mol_idx<train_no; mol_idx++)
    {
        molecule2descriptor(train_mol[mol_idx],dset->train_desc[mol_idx]);
        dset->train_val[mol_idx] = train_mol[mol_idx]->energy;
    }
    #pragma omp parallel for schedule(dynamic)
    for (int mol_idx=0; mol_idx<validate_no; mol_idx++)
    {
        molecule2descriptor(validate_mol[mol_idx],dset->validate_desc[mol_idx]);
        dset->validate_val[mol_idx] = validate_mol[mol_idx]->energy;
    }
    free_mol_array(mol_arr, molecules_no);
    return dset;
}

Dataset *create_dataset(int train_no, int validate_no)
{
    Dataset *dset = (Dataset *) malloc(sizeof(Dataset));
    dset->train_no = train_no;
    dset->validate_no = validate_no;
    dset->train_desc = create_descriptor_arr(train_no);
    dset->validate_desc = create_descriptor_arr(validate_no);
    dset->train_val = (double *) malloc(train_no*sizeof(double));
    dset->validate_val = (double *) malloc(validate_no*sizeof(double));
    return dset;
}

