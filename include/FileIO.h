#ifndef FILEIO_H
#define FILEIO_H

#include "State.h"
#include <time.h>

extern const size_t BUFFER_SIZE;
extern const size_t ARR_BUFFER_SIZE;
extern char outFile[260];

void safe_log(FILE *stream, const char *fmt, ...);

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

int match_atom_type(char *type, char *types[], int *num_types);

void input_logging(struct SimulationState *sim_state, struct SimulationEnv *sim_env,
                   struct LoggingState *log_state);
bool output_log_file(FILE *sim_log, int frame_num, unsigned long int iter, double elapsed_stime,
                     double temperature, double overpotential, long int atom_cnt,
                     double total_internal_energy);
void output_kmc_csv_header(FILE *csv_file);
void log_kmc_state_csv(FILE *csv_file, int frame_num, unsigned long int iter, double elapsed_stime,
                       double temperature, double overpotential, long int atom_cnt,
                       double total_internal_energy);
void output_mc_csv_header(FILE *csv_file);
void log_mc_state_csv(FILE *csv_log_file, const int frame_num, const unsigned long int mmc_steps,
                      const unsigned long int iter, const double sys_energy);
bool write_xyz_file(char *xyz_filename, int frame_num, char *suffix, struct SimulationState *ss,
                    struct SimulationEnv *se);
bool output_kmc_file(char *kmc_filename); // need to rework

void write_backlog(FILE *, FILE *);

void output_kmc_iter_header(FILE *csv_file);
void log_kmc_iter(FILE *csv_log_file, const unsigned long int mcss, const double sim_time,
                  const double sys_energy, const int uvw1[3], const int uvw2[3], int is_evap);
void output_mc_iter_header(FILE *csv_file);
void log_mc_iter(FILE *csv_log_file, const unsigned long int mcss, const double sys_energy,
                 const double deltaE, const int performed, const int uvw1[3], const int uvw2[3]);

// exposed for testing
int parse_input(char *line, FILE *temp_log, struct SimulationState *ss, struct SimulationEnv *se,
                struct LoggingState *ls);

#endif // FILEIO_H
