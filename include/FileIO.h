#ifndef FILEIO_H
#define FILEIO_H

#include "State.h"
#include <time.h>

extern const int BUFFER_SIZE;
extern const int ARR_BUFFER_SIZE;

bool simulation_parameters_from_file(char *filename, struct SimulationState *ss,
                                     struct SimulationEnv *se, struct LoggingState *ls,
                                     FILE *temp_log, time_t starttime); // need to modify!
bool process_xyz_file(FILE *temp_log, FILE *input_file, struct SimulationState *ss,
                      struct SimulationEnv *se, struct LoggingState *ls); // need to modify
bool process_kmc_file(FILE *temp_log, FILE *input_file, struct SimulationState *ss,
                      struct SimulationEnv *se, struct LoggingState *ls); // need to modify
bool process_in_file(FILE *temp_log, FILE *input_file, struct SimulationState *ss,
                     struct SimulationEnv *se, struct LoggingState *ls);
bool process_kmx_file(FILE *temp_log, FILE *input_file, struct SimulationState *ss,
                      struct SimulationEnv *se,
                      struct LoggingState *ls); // is this actually relevant?

int match_atom_type(char *type, char *types[], int *num_types, FILE *temp_log);

void input_logging(struct SimulationState *sim_state, struct SimulationEnv *sim_env,
                   struct LoggingState *log_state);
bool output_log_file(FILE *sim_log_file, int frame_num, unsigned long int iter,
                     double elapsed_stime, double temperature, double overpotential, int atom_cnt,
                     double total_internal_energy);
bool write_xyz_file(char *xyz_filename, int frame_num, char *suffix, struct SimulationState *ss,
                    struct SimulationEnv *se);
bool output_kmc_file(char *kmc_filename); // need to rework

void write_backlog(FILE *, FILE *);

// exposed for testing
int parse_input(char *line, FILE *temp_log, struct SimulationState *ss, struct SimulationEnv *se,
                struct LoggingState *ls);


#endif // FILEIO_H
