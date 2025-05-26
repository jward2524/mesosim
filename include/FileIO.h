#include "Common.h"

// IO
bool simulation_parameters_from_file(char *filename, struct SimulationState *ss); //need to modify!
bool process_xyz_file(FILE* temp_log, FILE* input_file); //need to modify
bool process_kmc_file(FILE* temp_log, FILE* input_file, struct SimulationState *ss); //need to modify
bool process_in_file(FILE* temp_log, FILE* input_file, struct SimulationState *ss);
bool process_kmx_file(FILE* temp_log, FILE* input_file, struct SimulationState *ss); //is this actually relevant?

int parse_input(char* line, struct SimulationState *ss);
int parse_boolean(char* str);
int match_atom_type(char* type, char* types[], int* num_types);
int parse_datalog_params(char* params, int cursor);

bool output_log_file(FILE* sim_log_file, int frame_num, struct SimulationState *ss);
bool write_xyz_file(char* xyz_filename, int frame_num, struct SimulationState *ss);
bool output_kmc_file(char *kmc_filename); //need to rework
