#include "Common.h"

void get_system_rw_radius(struct SimulationEnv *se);
void get_shifts(struct SimulationEnv *se);
void initialize_neighbor_offsets(struct SimulationEnv *se);
// void write_initial_information(FILE *datfil); //defined in Simulation Aux for some reason?
void initialize_jump_offsets(int l_t);
void initialize_zones(struct SimulationState *ss, struct SimulationEnv *se);
void get_system_normal(void);
void set_default_orientation(struct SimulationState *ss, struct SimulationEnv *se); //double check if this is used for anything other than graphicsw
int get_initial_configuration2(int at, int vc, int initial_config[], struct SimulationState *ss, struct SimulationEnv *se);
int get_final_configuration2(int at, int vc, int final_config[], struct SimulationState *ss, struct SimulationEnv *se);
void findzone(int *xz, int *yz, int *zz, int xxx, int yyy, int zzz, struct SimulationEnv *se);
void adjust_pbc(int *x, int *y, int *z);
void set_primitive_basis(int lt, struct SimulationState *ss);

void calculate_internal_energy(int atom_cnt, struct SimulationState *ss, struct SimulationEnv *se);

//bool do_initialize_simulation_parameters(int simulation_index); //modified AND needs additional arguments
void do_initialize_simulation(int simulation_index, struct SimulationState *ss, struct SimulationEnv *se);

//void initialize_simulation_parameters(void); //this might only be a file thing?
void general_simulation_initialization(struct SimulationState *ss, struct SimulationEnv *se); //will probably need to change arguments

void initialize_flat_sheet_1(int z, struct SimulationState *ss, struct SimulationEnv *se);
void initialize_spherical_cluster(int radius_of_sphere, struct SimulationState *ss, struct SimulationEnv *se);
void initialize_from_file(char* filename, struct SimulationState *ss, struct SimulationEnv *se);

void initialize_simulation_box(double system_size_x, double system_size_y, double system_size_z);
void check_pbc(int* x, int* y, int* z);

void pbc_translate(int coords_lat[3], int translation_vector[3]);
void corners2limits(double corners_cart[8][3], int limits_lat[3][2]);

//void initialize_vacancies_1(void); //might be good to have
