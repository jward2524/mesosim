// IO
bool get_input_file(char *filename); //need to modify!
bool process_xyz_file(FILE* temp_log, char *xyz_filename); //need to modify
bool process_kmc_file(FILE* temp_log, char *kmc_filename); //need to modify
bool process_in_file(FILE* temp_log, char* in_filename);
bool process_kmx_file(FILE* temp_log, char* kmx_filename); //is this actually relevant?

int parse_input(char* line);
int parse_boolean(char* str);
int match_atom_type(char* type, char* types[], int* num_types);

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