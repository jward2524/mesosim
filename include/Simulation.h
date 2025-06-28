#ifndef SIMULATION_H
#define SIMULATION_H

#include "Common.h"
#include <stdbool.h>

unsigned long perform_simulation(struct SimulationState* ss, struct SimulationEnv* se, struct LoggingState* ls);
void compute_transition_array(struct SimulationState*);
int refresh_transitions(int atom_idx, struct SimulationState* ss, struct SimulationEnv* se);
int is_on_transition_list(struct SimulationState* ss, double rate) ;
int create_new_transition(struct SimulationState* ss, double rate);
void add_to_transition_list(int rate_idx, int atom_idx, int offset_idx, struct SimulationState* ss, struct SimulationEnv* se);
void take_off_transition_list(int atom_idx, int offset_idx, struct SimulationState* ss);
void check_system(struct SimulationState* ss, struct SimulationEnv* se);
int calculate_surf_diffusion_rate(
    int initial_configuration[],
    int final_configuration[],
    int number_of_neighbors,
    int atom_type,
    double nnE[6],
    double temperature,
    double overpotential,
    double *rate
);

int	calculate_evaporation_rate(
    int initial_configuration[],
    int number_of_neighbors,
    int atom_type,
    double nnE[6],
    bool solubility[3],
    double temperature,
    double overpotential,
    double *rate
);

#endif // SIMULATION_H
