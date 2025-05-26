// struct SimulationEnv
// {
//     int simulation_type;
//     bool evaporation_flag;
//     // double total_internal_energy;

//     int zixshift, ziyshift, zizshift;
//     int ssxshift, ssyshift, sszshift;
//     int zsh, ysh, xsh;							// shifts

//     // [ ]: what are the units for this? how does it relate to atomic spacing?
//     int ssx, ssy, ssz;					// system size x, y, z
//     double ssr;
//     int zix, ziy, ziz;

//     int lattice_type;
//     int max_neighbors;

//     int sheet_thickness;
//     int cluster_radius;
//     char atoms_filename[256];

//     double initial_overpotential;
//     double overpotential_ramp_rate;
//     double max_overpotential;

//     double substrate_percent_a;
//     double substrate_percent_b;

//     // double log_interval;
//     //double vacancy_density = 0.01; //can always add back in
//     //int ncsk = 0;

//     double normal_x, normal_y, normal_z; // XXX: likely vistigal

//     int simbox_limits_lat[3][2];

//     double nnE[6];
//     double nnnE[6];
    
//     bool solubility[3]; //all elements cannot dissolve by default
    
//     double temperature;
//     double overpotential;

//     int dissolution; // flag for whether dissolution events can occur
//     char atom_names[3][3];

// };

#ifndef SIMULATION_AUX_H
#define SIMULATION_AUX_H

#include "Common.h"

void get_system_rw_radius(struct SimulationEnv* se);
void get_shifts(struct SimulationEnv* se);
void initialize_neighbor_offsets(int lattice_type, int* max_neighbors);

void initialize_jump_offsets(int lattice_type);
void initialize_zones(Zone zone_arr[ZONES_IN_X][ZONES_IN_Y][ZONES_IN_Z], struct SimulationEnv* se);
void get_system_normal(void);
void set_default_orientation(Atom *atom_arr[], int atom_cnt, int lattice_type); //double check if this is used for anything other than graphicsw
int get_initial_configuration2(int atom_index, int offset_index, int max_neighbors, Atom *atom_arr[], int initial_config[]);
int get_final_configuration2(int at, int offset_idx, struct SimulationState *ss, struct SimulationEnv *se, int final_config[]);
void findzone(int *xz, int *yz, int *zz, int xxx, int yyy, int zzz, struct SimulationEnv* se);
void adjust_pbc(int *x, int *y, int *z);
void set_primitive_basis(int lt);

void calculate_internal_energy(Atom *atom_arr[], int atom_cnt, double nnE[6], int max_neighbors, double* total_internal_energy);

void do_initialize_simulation(struct SimulationState* ss, struct SimulationEnv* se);

void general_simulation_initialization(struct SimulationState* ss, struct SimulationEnv* se, struct LoggingState* ls); //will probably need to change arguments

void initialize_flat_sheet_1(struct SimulationState *ss, struct SimulationEnv* se);
void initialize_spherical_cluster(struct SimulationState *ss, struct SimulationEnv* se);
void initialize_from_file(char* filename);

void initialize_simulation_box(struct SimulationEnv* se);
void check_pbc(int* x, int* y, int* z);

void pbc_translate(int coords_lat[3], int translation_vector[3]);
void corners2limits(double corners_cart[8][3], int limits_lat[3][2]);

#endif // SIMULATION_AUX_H
