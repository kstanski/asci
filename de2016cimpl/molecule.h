#ifndef MOLECULE_H_INCLUDED
#define MOLECULE_H_INCLUDED

#define MAX_ATOMS 23
#define ATOM_TYPES 5
#define DIMENSIONS 3

#include <boost/geometry.hpp>
namespace bg = boost::geometry;
typedef bg::model::point<double, DIMENSIONS, bg::cs::cartesian> Position;

typedef struct molecule
{
    double energy;
    int atoms_no;
    int atom_types[MAX_ATOMS];  // index of [H,C,N,O,S]
    int types_total[ATOM_TYPES];
    Position ff_coords[MAX_ATOMS];
    Position dft_coords[MAX_ATOMS];
} Molecule;

Molecule **read_molecules(const char *filename, int molecules_no);
int free_mol_array(Molecule **mol_arr, int molecules_no);
int type2index(char *type);
Position make_position(double *val_arr);
bool compare_molecules(Molecule *mol1, Molecule *mol2);

#endif // MOLECULE_H_INCLUDED
