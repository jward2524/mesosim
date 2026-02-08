#ifndef ERROR_H
#define ERROR_H
#include "State.h"

#ifdef TEST
// set when exit is called during tests
extern int exit_flag;
extern int exit_errno;
// expose to check whether they are correctly set in tests
extern struct SimulationState *sim_state;
extern struct SimulationEnv *sim_env;
extern struct LoggingState *log_state;
#endif

void set_state(struct SimulationState *ss, struct SimulationEnv *se, struct LoggingState *ls);
void initialize_states(struct SimulationState **ss, struct SimulationEnv **se,
                        struct LoggingState **ls);
void clean_and_error(int error);
void call_exit(int errnum);

#ifdef TEST
// expose for testing
void free_if_exists(void **pointer);
#endif

#ifdef HAVE_FORK
void create_coredump(void);
#endif

#endif // ERROR_H
