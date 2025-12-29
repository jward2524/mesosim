#ifndef MC_H
#define MC_H

#include "State.h"

unsigned long perform_metropolis_mc(struct SimulationState *ss, struct SimulationEnv *se,
                                    struct LoggingState *ls);

#endif // MC_H
