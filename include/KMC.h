#ifndef KMC_H
#define KMC_H

#include "State.h"

unsigned long perform_kmc(struct SimulationState *ss, struct SimulationEnv *se,
                          struct LoggingState *ls);

void compute_transition_array(struct SimulationState *ss, struct SimulationEnv *se);

#endif // KMC_H
