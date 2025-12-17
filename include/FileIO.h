#ifndef FILEIO_H
#define FILEIO_H

#include "Common.h"
#include <time.h>

extern const int BUFFER_SIZE;
bool simulation_parameters_from_file(char* filename, struct SimulationState* ss, struct SimulationEnv* se, struct LoggingState* ls, FILE* temp_log, time_t starttime); //need to modify!
bool process_xyz_file(FILE* temp_log, FILE* input_file, struct SimulationState* ss, struct SimulationEnv* se, struct LoggingState* ls); //need to modify
bool process_kmc_file(FILE* temp_log, FILE* input_file, struct SimulationState* ss, struct SimulationEnv* se, struct LoggingState* ls); //need to modify
bool process_in_file(FILE* temp_log, FILE* input_file, struct SimulationState* ss, struct SimulationEnv* se, struct LoggingState* ls);
bool process_kmx_file(FILE* temp_log, FILE* input_file, struct SimulationState* ss, struct SimulationEnv* se, struct LoggingState* ls); //is this actually relevant?

int parse_input(char* line, FILE* temp_log, struct SimulationState* ss, struct SimulationEnv* se, struct LoggingState* ls);
int parse_boolean(char* str);
int match_atom_type(char* type, char* types[], int* num_types, FILE* temp_log);
int parse_datalog_params(char* params, int cursor, struct LoggingState* ls, FILE* temp_log);

void input_logging(struct SimulationState* sim_state, struct SimulationEnv* sim_env, struct LoggingState* log_state);
bool output_log_file(FILE* sim_log_file, int frame_num, unsigned long int iter, double elapsed_stime, double temperature, double overpotential, int atom_cnt, double total_internal_energy);
bool write_xyz_file(char* xyz_filename, int frame_num, char* suffix, struct SimulationState* ss, struct SimulationEnv* se);
bool output_kmc_file(char *kmc_filename); //need to rework

int get_num_bond_types(int num_elements);
int get_env_index(int nn, int bond_idx, struct SimulationEnv* se);
int get_bond_index(int a, int b, struct SimulationEnv* se);

// key-value pairs
struct KV 
{
	// malloc'd in parse_key_value
	char *key;
	char *value;
};

// Property descriptor
typedef struct {
    char name[32];   // e.g. "species", "pos", "force"
    char type;       // 'S' (string), 'R' (real), 'I' (int)
    int ncols;       // number of columns for this property
} PropertyDesc;

int parse_key_value(const char *kv_str, size_t kv_len, struct KV *kv);
int parse_comment(const char *line, struct KV **outpairs, size_t *outpairs_cnt);
int parse_properties_value(const char *propval, PropertyDesc **out_props, int *out_nprops);
int tokenize_line(char *line, char **tokens, int maxtok);
int fill_atom_from_tokens(Atom *atom, char **tokens, int ntokens, PropertyDesc *props, int nprops);
void fill_atom_from_xyz(Atom *a, char **tokens, int ntok);

#endif // FILEIO_H
