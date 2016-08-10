#ifndef MOLECULE_H_INCLUDED
#define MOLECULE_H_INCLUDED

#define MAX_ATOMS 23
#define ATOM_TYPES 5
#define DIMENSIONS 3

typedef struct molecule
{
    double energy;
    int atoms_no;
    int atom_types[MAX_ATOMS];  // index of [H,C,N,O,S]
    int types_total[ATOM_TYPES];
    double ff_coords[MAX_ATOMS][DIMENSIONS];  // 3D
    double dft_coords[MAX_ATOMS][DIMENSIONS];
} Molecule;

Molecule *read_molecules(const char *filename, int molecules_no);
int free_mol_array(Molecule *mol_arr);
int type2index(char *type);
bool compare_molecules(Molecule mol1, Molecule mol2);

#endif // MOLECULE_H_INCLUDED
