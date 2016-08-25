#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "molecule.h"

Molecule **read_molecules(const char *filename, int molecules_no)
{
    /* create molecule array */
    Molecule **mol_arr = (Molecule **) malloc(molecules_no*sizeof(Molecule *));
    for (int i=0; i<molecules_no; i++)
    {
        mol_arr[i] = (Molecule *) malloc(sizeof(Molecule));
    }

    /* open file for reading */
    FILE *fp;
    fp = fopen(filename, "r");
    if(fp == NULL)
    {
        fprintf(stderr,"Error opening file\n");
        exit (1);
    }

    /* read line by line */
    char buf[256];
    int count_molecules = 0;
    int count_atoms = 0;
    Molecule *molecule = *mol_arr;   // points to the first Molecule.

    while (fgets(buf, sizeof(buf), fp) && count_molecules < molecules_no)
    {
        /* read words in a line */
        char *words[7]; // max 7 words in a line
        int word_count = 0;
        char *word;
        word = strtok(buf," ");
        while (word != NULL)
        {
            words[word_count] = word;
            word = strtok(NULL, " ");
            word_count++;
        }

        /* extract data from each line */
        switch (word_count)
        {
        case 1:
            count_atoms = 0;
            memset(molecule->types_total,0,ATOM_TYPES);
            molecule->atoms_no = atoi(words[0]);
            break;
        case 2:
            strncpy(molecule->id,words[0],ID_LEN);
            molecule->energy = atof(words[1]);
            break;
        case 7:
            molecule->atom_types[count_atoms] = type2index(words[0]);
            double ff_pos[DIMENSIONS], dft_pos[DIMENSIONS];
            for (int i=0; i<DIMENSIONS; i++)
            {
                ff_pos[i] = atof(words[i+1]);
                dft_pos[i] = atof(words[i+1+DIMENSIONS]);
            }
            molecule->ff_coords[count_atoms] = make_position(ff_pos);
            molecule->dft_coords[count_atoms] = make_position(dft_pos);
            molecule->types_total[type2index(words[0])]++;

            if (++count_atoms == molecule->atoms_no) // all atoms processed
            {
                count_molecules++;
                molecule = mol_arr[count_molecules]; // move pointer to the next Molecule.
            }
            break;
        default :
            fprintf(stderr,"Invalid entry\n");
        }
    }

    if (ferror(fp))
    {
        fprintf(stderr,"Error reading the file\n");
        exit (1);
    }
    return mol_arr;
}

int type2index(char *type)
{
    if (strcmp("H",type) == 0) return 0;
    if (strcmp("C",type) == 0) return 1;
    if (strcmp("N",type) == 0) return 2;
    if (strcmp("O",type) == 0) return 3;
    if (strcmp("S",type) == 0) return 4;
    fprintf(stderr,"Unknown atom type\n");
    return -1;
}

Position make_position(double *val_arr)
{
    Position pos(val_arr[0],val_arr[1],val_arr[2]);
    return pos;
}

bool compare_molecules(Molecule *mol1, Molecule *mol2)
{
    return (*mol1).energy > (*mol2).energy;   //since energy is negative
}

int free_mol_array(Molecule **mol_arr, int molecules_no)
{
    for (int i=0; i<molecules_no; i++)
    {
        free(mol_arr[i]);
    }
    free(mol_arr);
    return 0;
}
