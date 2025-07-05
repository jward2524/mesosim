#ifndef SIMULATION_H
#define SIMULATION_H

#include "Common.h"
#include <stdbool.h>

unsigned long perform_simulation(struct SimulationState* ss, struct SimulationEnv* se, struct LoggingState* ls);
void compute_transition_array(struct SimulationState*);
int refresh_transitions(int atom_idx, struct SimulationState* ss, struct SimulationEnv* se);
int is_on_rate_list(struct SimulationState* ss, double rate) ;
int create_new_rate(double rate, struct SimulationState* ss, struct SimulationEnv* se);
void add_to_transition_list(int rate_idx, int atom_idx, int offset_idx, struct SimulationState* ss, struct SimulationEnv* se);
void take_off_transition_list(int atom_idx, int offset_idx, struct SimulationState* ss);
void check_system(struct SimulationState* ss, struct SimulationEnv* se);
int calculate_surf_diffusion_rate(
    int initial_configuration[],
    int final_configuration[],
    int atom_type,
    double temperature,
    double overpotential,
    double *rate,
    struct SimulationEnv* se
);

int	calculate_evaporation_rate(
    int initial_configuration[],
    int atom_type,
    double temperature,
    double overpotential,
    double *rate,
    struct SimulationEnv* se
);

int get_bond_index(int a, int b, struct SimulationEnv* se);

#endif // SIMULATION_H
