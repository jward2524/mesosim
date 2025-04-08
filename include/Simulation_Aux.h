void get_system_rw_radius(void);
void getshifts(void);
void initialize_neighbor_offsets(void);
// void write_initial_information(FILE *datfil); //defined in Simulation Aux for some reason?
void initialize_jump_offsets(int l_t);
void initialize_zones(void);
void get_system_normal(void);
void set_default_orientation(void); //double check if this is used for anything other than graphicsw
int get_initial_configuration2(int at, int vc, int initial_config[]);
int get_final_configuration2(int at, int vc, int final_config[]);
void findzone(int *xz, int *yz, int *zz, double xxx, double yyy, double zzz);
void adjust_pbc(double *x, double *y, double *z);
void set_latmat(int lt);

void calculate_internal_energy(int nat);

//bool do_initialize_simulation_parameters(int simulation_index); //modified AND needs additional arguments
void do_initialize_simulation(int simulation_index);

//void initialize_simulation_parameters(void); //this might only be a file thing?
void general_simulation_initialization(void); //will probably need to change arguments

void initialize_flat_sheet_1(int z);
void initialize_spherical_cluster_1(int radius_of_sphere);
void initialize_from_file(char* filename);

//void initialize_vacancies_1(void); //might be good to have