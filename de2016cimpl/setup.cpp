#include <iostream>
#include <algorithm>

#include "setup.h"
#include "stratify.h"
#include "molecule.h"
#include "descriptor.h"

#define BELOW_5_NON_H 43

Dataset *setup(const char *filename, int molecules_no, int train_no, int validate_no)
{
    Dataset *dset = create_dataset(train_no,validate_no);
    Molecule **mol_arr = read_molecules(filename,molecules_no);

    /* stratify and divide into training and validation arrays */
    std::sort(mol_arr, mol_arr+molecules_no, compare_molecules);
    /* molecules below 5 non-H atoms (lowest energy) all go into training set */
    Molecule *train_mol[train_no];
    Molecule *validate_mol[validate_no];
    if (train_no <= BELOW_5_NON_H)
        stratify(mol_arr, train_mol, train_no, validate_mol, validate_no);
    else
    {
        int large_mol_no = molecules_no - BELOW_5_NON_H;
        Molecule *large_mol_arr[large_mol_no];
        for (int mol_idx=0; mol_idx<molecules_no; mol_idx++)
        {
            if (mol_idx<BELOW_5_NON_H)
                train_mol[mol_idx] = mol_arr[mol_idx];
            else large_mol_arr[mol_idx-BELOW_5_NON_H] = mol_arr[mol_idx];
        }

        int large_train_no = train_no - BELOW_5_NON_H;
        int large_mol_sample_no = large_train_no + validate_no;
        Molecule *take_mol[large_mol_sample_no];
        Molecule **sample_large_mol;
        if (0<large_mol_sample_no && large_mol_sample_no<large_mol_no)
        {
            int leave_no = large_mol_no - large_mol_sample_no;
            Molecule *leave_mol[leave_no];
            stratify(large_mol_arr, take_mol, large_mol_sample_no, leave_mol, leave_no);
            sample_large_mol = take_mol;
        } else if (large_mol_sample_no == large_mol_no)
            sample_large_mol = large_mol_arr;
        else exit(1);

        Molecule *large_train_mol[large_train_no];
        stratify(sample_large_mol, large_train_mol, large_train_no, validate_mol, validate_no);
        for (int mol_idx=0; mol_idx<large_train_no; mol_idx++)
            train_mol[mol_idx+BELOW_5_NON_H] = large_train_mol[mol_idx];
    }
#if VERBOSE
    std::cout << "computing power spectra descriptors" << std::endl;
#endif // VERBOSE
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

int free_dataset(Dataset *dset)
{
    free_desc_arr(dset->train_desc,dset->train_no);
    free_desc_arr(dset->validate_desc,dset->validate_no);
    free(dset->train_val);
    free(dset->validate_val);
    free(dset);
    return 0;
}

