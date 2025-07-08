#ifndef SIMULATION_H
#define SIMULATION_H

#include "Common.h"
#include <stdbool.h>

unsigned long perform_simulation(struct SimulationState* ss, struct SimulationEnv* se, struct LoggingState* ls);
void compute_transition_array(struct SimulationState* ss, struct SimulationEnv* se);
int refresh_transitions(int atom_idx, struct SimulationState* ss, struct SimulationEnv* se);
int get_rate(unsigned char* atom_env, unsigned char is_evaporation, struct SimulationState* ss, struct SimulationEnv* se);
int create_new_rate(unsigned char *atom_env, unsigned char is_evaporation, struct SimulationState* ss, struct SimulationEnv* se);
void add_to_transition_list(int rate_idx, int atom_idx, int offset_idx, struct SimulationState* ss, struct SimulationEnv* se);
void take_off_transition_list(int atom_idx, int offset_idx, struct SimulationState* ss);
void check_system(struct SimulationState* ss, struct SimulationEnv* se);
int get_bond_index(int a, int b, struct SimulationEnv* se);

double calculate_evaporation_rate(
	unsigned char *atom_env,
	double temperature,						// system temperature
	double overpotential,					// system overpotential
	struct SimulationEnv* se
);

double calculate_surf_diffusion_rate(
	unsigned char* atom_env,
	double temperature,						// system temperature
	// double overpotential,					// system overpotential
	struct SimulationEnv* se
);

#endif // SIMULATION_H
