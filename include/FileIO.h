// IO
bool get_input_file(char *filename); //need to modify!
bool process_xyz_file(char *xyz_filename); //need to modify
bool process_kmc_file(char *kmc_filename); //need to modify
bool process_in_file(char* in_filename);
bool process_kmx_file(char* kmx_filename); //is this actually relevant?

int parse_input(char* line);
int parse_boolean(char* str);
int match_atom_type(char* type, char* types[], int* num_types);

bool output_log_file(int frame_num);
bool write_xyz_file(char* xyz_filename, int frame_num);
bool output_kmc_file(char *kmc_filename); //need to rework