#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "molecule.h"

Molecule *read_molecules(const char *filename, int molecules_no)
{
    /* create molecule array */
    Molecule *mol_arr;
    mol_arr = (Molecule *) malloc(molecules_no*sizeof(Molecule));
    if (mol_arr == NULL)
    {
        fprintf(stderr,"Not enough memory: molecules");
        exit (1);
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
    Molecule *molecule = mol_arr;   // points to the first Molecule.

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
            molecule->energy = atof(words[1]);
            break;
        case 7:
            molecule->atom_types[count_atoms] = type2index(words[0]);
            for (int i=0; i<3; i++)
            {
                molecule->ff_coords[count_atoms][i] = atof(words[i+1]);
                molecule->dft_coords[count_atoms][i] = atof(words[i+4]);
            }

            molecule->types_total[type2index(words[0])]++;

            if (++count_atoms == molecule->atoms_no) // all atoms processed
            {
                count_molecules++;
                molecule++; // move pointer to the next Molecule.
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

int free_mol_array(Molecule *mol_arr)
{
    int len = sizeof(mol_arr)/sizeof(Molecule);
    for (int idx=0; idx<len; idx++)
    {
        free(&mol_arr[idx]);
    }
    free(mol_arr);
    return 0;
}
