#ifndef FILEIO_H
#define FILEIO_H

#include "State.h"

extern const size_t BUFFER_SIZE;
extern const size_t ARR_BUFFER_SIZE;

typedef struct {
    const char *name;
    CsvFieldFuncPtr get_value;
    unsigned flavor;
} CsvFieldFunc;
extern const CsvFieldFunc csv_field_map[];
extern const size_t CSV_FIELD_FUNCS_COUNT;

void safe_log(FILE *stream, const char *fmt, ...);

bool simulation_parameters_from_file(char *filename, struct SimulationState *ss,
                                     struct SimulationEnv *se, struct LoggingState *ls);
bool process_xyz_file(FILE *input_file, struct SimulationState *ss, struct SimulationEnv *se,
                      struct LoggingState *ls);
void process_in_file(FILE *input_file, struct SimulationState *ss, struct SimulationEnv *se,
                     struct LoggingState *ls);

void input_logging(struct SimulationState *sim_state, struct SimulationEnv *sim_env,
                   struct LoggingState *log_state);
bool output_log_file(FILE *sim_log, int frame_num, unsigned long int iter, double elapsed_stime,
                     double temperature, double overpotential, long int atom_cnt,
                     double total_internal_energy);
// void output_kmc_csv_header(FILE *csv_file);
// void log_kmc_state_csv(FILE *csv_file, int frame_num, unsigned long int iter, double elapsed_stime,
//                        double temperature, double overpotential, long int atom_cnt,
//                        double total_internal_energy);
// void output_mc_csv_header(FILE *csv_file);
// void log_mc_state_csv(FILE *csv_log_file, const int frame_num, const unsigned long int mmc_steps,
//                       const unsigned long int iter, const double sys_energy);
bool write_xyz_file(char *xyz_filename, int frame_num, char *suffix, struct SimulationState *ss,
                    struct SimulationEnv *se);

void output_csv_header(FILE *csv_file, struct LoggingState *ls);
void log_state_csv(FILE *csv_file, struct SimulationState *ss, struct LoggingState *ls);

void write_backlog(FILE *, FILE *);

void output_kmc_steps_header(FILE *csv_file, const bool output_coord);
void log_kmc_steps(FILE *csv_file, const unsigned long int iter, const double sim_time,
                   const double sys_energy, const int uvw1[3], const int uvw2[3], const int is_evap,
                   const int coordination);
void output_mc_steps_header(FILE *csv_file, bool output_coord);
void log_mc_steps(FILE *csv_file, const unsigned long int iter, const double sys_energy,
                  const double deltaE, const int performed, const int uvw1[3], const int uvw2[3],
                  const int coordination);

void write_xyz_suffix(char *suffix, OutputScheduleMode mode, double checkpoint);
void write_logs(int output_csv, int output_xyz, struct SimulationState *ss, struct SimulationEnv *se, struct LoggingState *ls);
void output_if_passed_checkpoint(struct SimulationState *ss, struct SimulationEnv *se, struct LoggingState *ls);


#endif // FILEIO_H
