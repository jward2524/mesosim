#ifndef ERROR_H
#define ERROR_H
#include "State.h"
#include "setjmp.h"

#ifdef TEST
// set when exit is called during tests
extern volatile int exit_flag;
extern volatile int exit_errno;
extern volatile int jmp_set;
extern int expected_exit_errno;
extern jmp_buf test_exit_jmp;
// expose to check whether they are correctly set in tests
extern struct SimulationState **gp_sim_state;
extern struct SimulationEnv **gp_sim_env;
extern struct LoggingState **gp_log_state;
#endif

void set_state(struct SimulationState **p_ss, struct SimulationEnv **p_se,
               struct LoggingState **p_ls);
void initialize_states(struct SimulationState **ss, struct SimulationEnv **se,
                       struct LoggingState **ls);
void clean_and_error(int error);
void call_exit(int errnum);

// utility function for heap cleanup
void free_if_exists(void **pointer);

#ifdef DUMP_CORE
void create_coredump(void);
#endif

#endif // ERROR_H
