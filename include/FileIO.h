#ifndef FILEIO_H
#define FILEIO_H

#include "Common.h"
#include <time.h>

bool simulation_parameters_from_file(char* filename, struct SimulationState* ss, struct SimulationEnv* se, struct LoggingState* ls, FILE* temp_log, time_t starttime); //need to modify!
bool process_xyz_file(FILE* temp_log, FILE* input_file, struct SimulationState* ss, struct SimulationEnv* se, struct LoggingState* ls); //need to modify
bool process_kmc_file(FILE* temp_log, FILE* input_file, struct SimulationState* ss, struct SimulationEnv* se, struct LoggingState* ls); //need to modify
bool process_in_file(FILE* temp_log, FILE* input_file, struct SimulationState* ss, struct SimulationEnv* se, struct LoggingState* ls);
bool process_kmx_file(FILE* temp_log, FILE* input_file, struct SimulationState* ss, struct SimulationEnv* se, struct LoggingState* ls); //is this actually relevant?

int parse_input(char* line, FILE* temp_log, struct SimulationState* ss, struct SimulationEnv* se, struct LoggingState* ls);
int parse_boolean(char* str);
int match_atom_type(char* type, char* types[], int* num_types, FILE* temp_log);
int parse_datalog_params(char* params, int cursor, struct LoggingState* ls, FILE* temp_log);

bool output_log_file(FILE* sim_log_file, int frame_num, double elapsed_stime, double temperature, double overpotential, int atom_cnt, double total_internal_energy);
bool write_xyz_file(struct SimulationState* sim_state, char* xyz_filename, int frame_num);
bool output_kmc_file(char *kmc_filename); //need to rework

int get_num_bond_types(int num_elements);
int get_env_index(int nn, int bond_idx, struct SimulationEnv* se);

#endif // FILEIO_H
