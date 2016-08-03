#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_ATOMS 23
#define ATOM_TYPES 5
#define DIMENSIONS 3
#define MAX_H 16
#define MAX_C 7
#define MAX_N 3
#define MAX_O 3
#define MAX_S 1
#define CUTOFF 3

typedef struct molecule
{
    float energy;
    int atoms_no;
    int atom_types[MAX_ATOMS];  // index of [H,C,N,O,S]
    int types_total[ATOM_TYPES];
    float ff_coords[MAX_ATOMS][DIMENSIONS];  // 3D
    float dft_coords[MAX_ATOMS][DIMENSIONS];
} Molecule;

typedef struct neighbourhood
{
    float coords[ATOM_TYPES][MAX_ATOMS][DIMENSIONS];    // 5 atom species
    int last_atom_idx[ATOM_TYPES] = {-1};
} Neighbourhood;

int type2index(char *type);
Neighbourhood *molecule2neighbourhoods(Molecule *mol_ptr);
int free_nhood_array(Neighbourhood *nhood_arr);

int read_data(const char *filename, int molecules_no)
{
    /* open file for reading */
    FILE *fp;
    fp = fopen(filename, "r");
    if(fp == NULL)
    {
        fprintf(stderr,"Error opening file");
        return(-1);
    }

    /* read line by line */
    char buf[256];
    int count_molecules = 0;
    int count_atoms = 0;
    Molecule molecule;

    while (fgets(buf, sizeof(buf), fp))
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
            count_molecules++;
            count_atoms = 0;
            memset(molecule.types_total,0,ATOM_TYPES);
            molecule.atoms_no = atoi(words[0]);
            break;
        case 2:
            molecule.energy = atof(words[1]);
            break;
        case 7:
            molecule.atom_types[count_atoms] = type2index(words[0]);
            for (int i=0; i<3; i++)
            {
                molecule.ff_coords[count_atoms][i] = atof(words[i+1]);
                molecule.dft_coords[count_atoms][i] = atof(words[i+4]);
            }

            molecule.types_total[type2index(words[0])]++;

            if (++count_atoms == molecule.atoms_no) // all atoms processed
            {
                Neighbourhood *nhoods = molecule2neighbourhoods(&molecule);
                for (int idx=0; idx<molecule.atoms_no; idx++)
                {
                    printf("%f ",nhoods[0].coords[0][idx][0]);
                    printf("%f ",nhoods[0].coords[0][idx][1]);
                    printf("%f\n",nhoods[0].coords[0][idx][2]);
                }
                printf("%d\n",nhoods[0].last_atom_idx[1]);
                printf("\n");
                free_nhood_array(nhoods);
            }
            break;
        default :
            fprintf(stderr,"Invalid entry\n");
        }
    }

    if (ferror(fp))
    {
        fprintf(stderr,"Error reading the file\n");
        abort();
    }
    return 0;
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

Neighbourhood *molecule2neighbourhoods(Molecule *mol_ptr)
{
    /* create neighbourhoods array */
    Neighbourhood *nhood_arr;
    int atoms_no = mol_ptr->atoms_no;
    int max_total = MAX_H+MAX_C+MAX_N+MAX_O+MAX_S;
    int max_types[] = {MAX_H,MAX_C,MAX_N,MAX_O,MAX_S};
    nhood_arr = (Neighbourhood *) malloc(max_total*sizeof(Neighbourhood));
    if (nhood_arr == NULL) exit (1);

    for (int idx=0; idx<max_total; idx++)
    {
        //memset(nhood_arr[idx].last_atom_idx,-1,ATOM_TYPES);
    }

    /* search for close neighbours */
    float c_i[DIMENSIONS];
    for (int i=0; i<atoms_no; i++)  // for center atom
    {
        for (int idx;idx<DIMENSIONS;idx++) c_i[idx] = mol_ptr->ff_coords[i][idx];
        for (int j=i; j<atoms_no; j++)  // for neighbouring atom
        {

            float norm = 0;
            float c_j[DIMENSIONS];
            for (int idx;idx<DIMENSIONS;idx++)
            {
                c_j[idx] = mol_ptr->ff_coords[j][idx];
                norm += pow((c_i[idx] - c_j[idx]),2);
            }
            norm = sqrt(norm);

            if (norm < CUTOFF)
            {
                /* add neighbour of i */
                int neighbour_type = mol_ptr->atom_types[j];
                nhood_arr[i].last_atom_idx[neighbour_type]++;
                int last = nhood_arr[i].last_atom_idx[neighbour_type];
                printf("%d\n",last);
                for (int idx;idx<DIMENSIONS;idx++)
                {
                    nhood_arr[i].coords[neighbour_type][last][idx] = c_j[idx] - c_i[idx];
                }

                if (i != j)
                {
                    /* add neighbour of j */
                    neighbour_type = mol_ptr->atom_types[i];
                    nhood_arr[j].last_atom_idx[neighbour_type]++;
                    last = nhood_arr[j].last_atom_idx[neighbour_type];
                    for (int idx;idx<DIMENSIONS;idx++)
                    {
                        nhood_arr[j].coords[neighbour_type][last][idx] = c_i[idx] - c_j[idx];
                    }
                }
            }
        }
    }

    /* top up with dummy atoms */
    int last = atoms_no - 1;
    for (int type=0; type<ATOM_TYPES; type++)
    {
        int lack = max_types[type] - (mol_ptr->types_total)[type];
        while (lack > 0)
        {
            last++;
            nhood_arr[last].last_atom_idx[type] = 0;
            for (int idx;idx<DIMENSIONS;idx++)
                nhood_arr[last].coords[type][0][idx] = 0;
            lack--;
        }
    }

    return nhood_arr;
}

int free_nhood_array(Neighbourhood *nhood_arr)
{
    int len = sizeof(nhood_arr)/sizeof(Neighbourhood);
    for (int idx=0; idx<len; idx++)
    {
        free(&nhood_arr[idx]);
    }
    free(nhood_arr);
    return 0;
}
