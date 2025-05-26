#include "Common.h"

void get_system_rw_radius(void);
void get_shifts(void);
void initialize_neighbor_offsets(void);
// void write_initial_information(FILE *datfil); //defined in Simulation Aux for some reason?
void initialize_jump_offsets(int l_t);
void initialize_zones(struct SimulationState *ss);
void get_system_normal(void);
void set_default_orientation(void); //double check if this is used for anything other than graphicsw
int get_initial_configuration2(int at, int vc, int initial_config[]);
int get_final_configuration2(int at, int vc, int final_config[], struct SimulationState *ss);
void findzone(int *xz, int *yz, int *zz, int xxx, int yyy, int zzz);
void adjust_pbc(int *x, int *y, int *z);
void set_primitive_basis(int lt);

void calculate_internal_energy(int atom_cnt, struct SimulationState *ss);

//bool do_initialize_simulation_parameters(int simulation_index); //modified AND needs additional arguments
void do_initialize_simulation(int simulation_index, struct SimulationState *ss);

//void initialize_simulation_parameters(void); //this might only be a file thing?
void general_simulation_initialization(struct SimulationState *ss); //will probably need to change arguments

void initialize_flat_sheet_1(int z, struct SimulationState *ss);
void initialize_spherical_cluster(int radius_of_sphere, struct SimulationState *ss);
void initialize_from_file(char* filename, struct SimulationState *ss);

void initialize_simulation_box(double system_size_x, double system_size_y, double system_size_z);
void check_pbc(int* x, int* y, int* z);

void pbc_translate(int coords_lat[3], int translation_vector[3]);
void corners2limits(double corners_cart[8][3], int limits_lat[3][2]);

//void initialize_vacancies_1(void); //might be good to have
