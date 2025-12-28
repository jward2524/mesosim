#ifndef UTILS_H
#define UTILS_H

#include "State.h"

#define DEFAULT_SEED 1

extern double normal_cart[6][3];
// extern double normal_lat[6][3];

double drand();

int get_env_index(int neighbor_idx, int atom_type, int neighbor_type, struct SimulationEnv *se);
int get_num_bond_types(int num_elements);
int nn_bondidx_2_envidx(int nn, int bond_idx, int num_bond_types);
int get_bond_index(int a, int b, int num_elements);

void get_system_rw_radius(struct SimulationEnv *se); // XXX: unused

void findzone(int *zone_u, int *zone_v, int *zone_w, int u, int v, int w, struct SimulationEnv *se);
void adjust_pbc(int *u, int *v, int *w, struct SimulationEnv *se);

void check_pbc(int *u, int *v, int *w, double basis[3][3]);

void pbc_translate(int coords_lat[3], int translation_vector[3]);

void cartesian2lattice_site(double ccart[3], double invert_primitive_basis[3][3], int clattice[3]);
void cartesian2lattice(double ccart[3], double invert_primitive_basis[3][3], double clattice[3]);
void lattice2cartesian(int clattice[3], double primitive_basis[3][3], double ccart[3]);

int round_towards(double val, int target);

// void cell_to_latmat(double c[6], double ltmt[3][3]);
void primitive_basis2ucell_params(double ltmt[3][3], double c[6]);

int get_type_from_name(char *atom_name, char **atom_names, int atom_names_cnt, unsigned char *atom_type);

#endif // UTILS_H
