#ifndef ERROR_H
#define ERROR_H
#include "State.h"

void set_state(struct SimulationState *ss, struct SimulationEnv *se, struct LoggingState *ls);
void initialize_states(struct SimulationState **ss, struct SimulationEnv **se,
                        struct LoggingState **ls);
void clean_and_exit(int error);

#endif // ERROR_H
