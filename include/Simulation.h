#include "Common.h"

unsigned long perform_simulation(struct SimulationState *ss, struct SimulationEnv *se);
void check_system(struct SimulationState *ss, struct SimulationEnv *se);
void take_off_transition_list(int at, int vc, struct SimulationState *ss);
int is_on_transition_list(double r, struct SimulationState *ss);
int create_new_transition(double r, struct SimulationState *ss);
void add_to_transition_list(int tl, int at, int vc, struct SimulationState *ss);
int refresh_transitions(int at, struct SimulationState *ss, struct SimulationEnv *se);
void compute_transition_array(struct SimulationState *ss);
int calculate_surf_diffusion_rate(	int initial_configuration[],
    int final_configuration[],
    int number_of_neighbors,
    int atom_type,
    double nnE[6],
    double temperature,
    double overpotential,
    double *rate
);

int	calculate_evaporation_rate( int initial_configuration[],
    int number_of_neighbors,
    int atom_type,
    double nnE[6],
    double temperature,
    double overpotential,
    double *rate,
    bool solubility[3]
);
