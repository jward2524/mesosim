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
void set_primitive_basis(int lt);

void calculate_internal_energy(int atom_cnt);

//bool do_initialize_simulation_parameters(int simulation_index); //modified AND needs additional arguments
void do_initialize_simulation(int simulation_index);

//void initialize_simulation_parameters(void); //this might only be a file thing?
void general_simulation_initialization(void); //will probably need to change arguments

void initialize_flat_sheet_1(int z);
void initialize_spherical_cluster_1(int radius_of_sphere);
void initialize_from_file(char* filename);

//void initialize_vacancies_1(void); //might be good to have

extern double total_internal_energy;

// [ ]: what are these?
extern int zixshift, ziyshift, zizshift;
extern int ssxshift, ssyshift, sszshift;
extern int zsh, ysh, xsh;							// shifts

// [ ]: what are the units for this? how does it relate to atomic spacing?
extern int ssx, ssy, ssz;					// system size x, y, z
extern double ssr;
extern int zix, ziy, ziz;

extern int lattice_type;
extern int number_of_possible_neighbors;

extern int sheet_thickness;
extern int cluster_radius;
extern char atoms_filename[256];

// [ ]: what is this for?
extern Zone zone_arr[ZONES_IN_X][ZONES_IN_Y][ZONES_IN_Z];

extern double initialoverpotential;
extern double overpotentialramprate;
extern double maxoverpotential;

extern double substrate_percent_a;
extern double substrate_percent_b;

extern double initial_logtime;

extern int analysis_type;

extern double logtime_multiplier;

//double vacancy_density = 0.01; //can always add back in

extern double overpotential_ramp_rate;

//int ncsk = 0;

extern int total_volume_dissolved;

extern double normal_x, normal_y, normal_z; // XXX: likely vistigal