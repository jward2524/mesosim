#ifndef INITIALIZATION_H
#define INITIALIZATION_H
#include "State.h"

void initialize_neighbor_offsets(struct SimulationEnv *se);
void initialize_zones(Zone ****zone_arr, struct SimulationEnv *se);
void set_primitive_basis(struct SimulationEnv *se);
void initialize_lattice_geometry(struct SimulationEnv *sim_env);
void initialize_initial_structure(struct SimulationState *ss, struct SimulationEnv *se,
                                  struct LoggingState *ls);
void initialize_simulation_variables(struct SimulationState *ss, struct SimulationEnv *se);
void initialize_flat_sheet(struct SimulationState *ss, struct SimulationEnv *se);
void initialize_spherical_cluster(struct SimulationState *ss, struct SimulationEnv *se);
void initialize_from_file(struct SimulationState *ss, struct SimulationEnv *se,
                          struct LoggingState *ls);
void initialize_simulation_box(struct SimulationEnv *se);
void get_shifts(struct SimulationEnv *se);
void get_system_normal(double basis[3][3]);

// double check if this is used for anything other than graphics
void set_default_orientation(int lattice_type, double rmat[3][3]);

void corners2limits(double corners_cart[8][3], int limits_lat[3][2], double inv_basis[3][3]);
void initialize_simulation(struct SimulationState *sim_state, struct SimulationEnv *sim_env,
                           struct LoggingState *log_state);
#endif // INITIALIZATION_H
