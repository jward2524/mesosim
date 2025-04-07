// prototypes

// simulation
unsigned long perform_simulation(void); //modified!

void get_system_rw_radius(void);

void getshifts(void);
void initialize_neighbor_offsets(void);
void write_initial_information(FILE *datfil); //defined in Simulation Aux for some reason?
void initialize_jump_offsets(int l_t);
void initialize_zones(void);

void get_system_normal(void);

void check_system(void);

int add_atom(double x, double y, double z, int type, int special);
void set_default_orientation(void); //double check if this is used for anything other than graphicsw

void take_off_transition_list(int at, int vc);
//int get_initial_configuration(int at, int vc, int initial_config[]); //keep both for now, only 2 is used
//int get_final_configuration(int at, int vc, int final_config[]);
int get_initial_configuration2(int at, int vc, int initial_config[]);
int get_final_configuration2(int at, int vc, int final_config[]);
int is_on_transition_list(double r);
int create_new_transition(double r);
void add_to_transition_list(int tl, int at, int vc);
int refresh_transitions(int at);

//all formerly inline methods but gcc is getting angry at me :(
void findzone(int *xz, int *yz, int *zz, double xxx, double yyy, double zzz);
void adjust_pbc(double *x, double *y, double *z);
int atom_at(double cx, double cy, double cz);

void write_final_information(FILE *datfil); //need to change this
//void copy_xyz_to_coord(void); //probably replace with organize
void compute_transition_array(void);
void remove_atom(int at);
void move_atom(int ia, int fa);

void make_buried_atoms_real(void);

int random_reincarnate_atom(double x, double y, double z, int type, int vc); //feels like it'll never be called
int reincarnate_atom(double x, double y, double z, int type, int vc); //will not be called
void bury_atom(int at, int *pos); //no longer relevant

// console
//void sparse(char original_string[], char string1[], char string2[], int comp_length); //keep this seems like a utility?
//int parse (int num, char cps[], int kind[], double defs[], double dest[]); //totally unused

// rotation and molecule orientation
void rotmat(Atom* atm[], int na, double rtmat[3][3]);
void organize(Atom* atm[], int n); //keep for now but need to retool
void orthomol(Atom* atm[], int na, double com[3][3]);
void centerg(Atom* atm[], int na); //not really relevant now

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


// general atom and bond handling routines
void copy_atom(int, int);
void create_default_atom(int na); //can modify this to remove things like color?

void kill_atom(int atom_number);

//void clear_chaff(void); //aspirational doing this

//void symop(double x[3], double smat[3][4], double y[3]);

void initialize_symmetry_elements(void);
void initialize_lattice_geometry(void);

void set_latmat(int lt);

void cell_to_latmat(double c[6], double ltmt[3][3]);
void latmat_to_cell(double ltmt[3][3], double c[6]);


int calculate_surf_diffusion_rate(	int initial_configuration[],
												int final_configuration[],
												int number_of_neighbors,
												int atom_type,
												double nnE[6],
												double temperature,
												double overpotential,
												double *rate);

int	calculate_evaporation_rate( int initial_configuration[],
												int number_of_neighbors,
												int atom_type,
												double nnE[6],
												double temperature,
												double overpotential,
												double *rate);

void calculate_internal_energy(int nat);

//bool do_initialize_simulation_parameters(int simulation_index); //modified AND needs additional arguments
void do_initialize_simulation(int simulation_index);

//void initialize_simulation_parameters(void); //this might only be a file thing?
void general_simulation_initialization(void); //will probably need to change arguments

void initialize_flat_sheet_1(int z);
void initialize_spherical_cluster_1(int radius_of_sphere);
void initialize_from_file(char* filename);

//void initialize_vacancies_1(void); //might be good to have