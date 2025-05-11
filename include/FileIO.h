// IO
bool simulation_parameters_from_file(char *filename); //need to modify!
bool process_xyz_file(FILE* temp_log, FILE* input_file); //need to modify
bool process_kmc_file(FILE* temp_log, FILE* input_file); //need to modify
bool process_in_file(FILE* temp_log, FILE* input_file);
bool process_kmx_file(FILE* temp_log, FILE* input_file); //is this actually relevant?

int parse_input(char* line);
int parse_boolean(char* str);
int match_atom_type(char* type, char* types[], int* num_types);
int parse_datalog_params(char* params, int cursor);

bool output_log_file(FILE* sim_log_file, int frame_num);
bool write_xyz_file(char* xyz_filename, int frame_num);
bool output_kmc_file(char *kmc_filename); //need to rework

/* IO related variables */
extern char console_outstring[512]; //keep for now; change when file writing changes
extern FILE *view_save_file;
extern FILE *temp_log;
extern char command_string[1024];
// extern char return_message[512];
extern char outFile[260];
extern char default_extension[];
extern FILE *sim_log_file;