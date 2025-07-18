#ifndef SIMULATION_AUX_H
#define SIMULATION_AUX_H

#include "Common.h"

void get_system_rw_radius(struct SimulationEnv* se);
void get_shifts(struct SimulationEnv* se);
void initialize_neighbor_offsets(int lattice_type, int* max_neighbors, CrystalOffset *jump_offset, int *opposite_offset);

void initialize_jump_offsets(int lattice_type, CrystalOffset *jump_offset, int *opposite_offset);
void initialize_zones(Zone zone_arr[ZONES_IN_X][ZONES_IN_Y][ZONES_IN_Z], struct SimulationEnv* se);
void get_system_normal(double basis[3][3]);
void set_default_orientation(Atom** atom_arr, int atom_cnt, int lattice_type, double rmat[3][3], double basis[3][3]); //double check if this is used for anything other than graphicsw
int get_initial_configuration(int atom_index, int max_neighbors, Atom** atom_arr, int initial_config[]);
int get_final_configuration(int at, int offset_idx, struct SimulationState *ss, struct SimulationEnv *se, int final_config[]);
void findzone(int *zone_u, int *zone_v, int *zone_w, int u, int v, int w, struct SimulationEnv* se);
void adjust_pbc(int *u, int *v, int *w, struct SimulationEnv* se);
void set_primitive_basis(struct SimulationEnv* se);

void calculate_internal_energy(Atom** atom_arr, int atom_cnt, double* total_internal_energy, struct SimulationEnv* se);

void initialize_initial_structure(struct SimulationState* ss, struct SimulationEnv* se);

void initialize_simulation_variables(struct SimulationState* ss, struct SimulationEnv* se);

void initialize_flat_sheet(struct SimulationState *ss, struct SimulationEnv* se);
void initialize_spherical_cluster(struct SimulationState *ss, struct SimulationEnv* se);
void initialize_from_file(char* filename);

void initialize_simulation_box(struct SimulationEnv* se);
void check_pbc(int* u, int* v, int* w, double basis[3][3]);

void pbc_translate(int coords_lat[3], int translation_vector[3]);
void corners2limits(double corners_cart[8][3], int limits_lat[3][2], double inv_basis[3][3]);

#endif // SIMULATION_AUX_H
