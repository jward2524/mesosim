#ifndef SIMULATION_H
#define SIMULATION_H

#include "State.h"
#include <stdbool.h>

int refresh_transitions(long atom_idx, struct SimulationState *ss, struct SimulationEnv *se);
int get_rate(unsigned char *atom_env, int final_config_neighbor_cnt, unsigned char is_evaporation,
             struct SimulationState *ss, struct SimulationEnv *se);
long create_new_rate(unsigned char *atom_env, int final_config_neighbor_cnt,
                     unsigned char is_evaporation, struct SimulationState *ss,
                     struct SimulationEnv *se);
void add_to_transition_list(long rate_idx, long atom_idx, unsigned char offset_idx,
                            struct SimulationState *ss, struct SimulationEnv *se);
void remove_from_transition_list(long atom_idx, int offset_idx, struct SimulationState *ss);

void update_outdated_transitions(int old_x, int old_y, int old_z, long transitioned_atom_idx,
                                 struct SimulationState *ss, struct SimulationEnv *se);

#endif // SIMULATION_H
