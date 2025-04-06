//extern bool evaporation_flag;

extern FILE *sim_log_file;


extern char coordinate_log_prefix[256];

//extern int number_of_simulation_runs;

extern int zixshift, ziyshift, zizshift;
extern int ssxshift, ssyshift, sszshift;
extern int zsh, ysh, xsh;							

extern int ssx, ssy, ssz;
extern double ssr;
extern int zix, ziy, ziz;

extern double run_time;
extern double data_time_interval;
extern double time_interval_end;

extern int lattice_type;
extern int number_of_possible_neighbors;

extern int sheet_thickness;
extern int cluster_radius;
extern char atoms_filename[256];

extern Zone zone[ZONES_IN_X][ZONES_IN_Y][ZONES_IN_Z];

extern int number_rates;
extern int total_current_transitions;
extern double sum_of_frequencies;

extern double overpotential;
extern double initialoverpotential, overpotentialramprate, maxoverpotential;

//extern double nnEaa, nnEab, nnEac, nnEbb, nnEbc, nnEcc;
extern double nnE[6];
extern double nnnE[6]; //second degree

/*extern double deposition_rate_of_a;
extern double deposition_rate_of_b;
extern double deposition_rate_of_c;*/

extern bool solubility[3]; //determines whether soluble or not

extern double temperature;


extern double substrate_percent_a;
extern double substrate_percent_b;
extern int dissolution;

extern int number_final_configuration_neighbors;
extern int number_intial_configuration_neighbors;

extern Rate rate[];
//extern Transition_List transition_list[MAXIMUM_NUMBER_OF_CONCURRENT_TRANSITIONS];
extern Transition_List *transition_list[MAXIMUM_NUMBER_OF_CONCURRENT_TRANSITIONS];
extern Trans_Prob transition_probability;

extern long int final_iteration;
extern double lastxt, lastyt, lastzt;

//extern bool simulation_is_going;

extern double sum_of_rate_populations;
extern double current_probability;

//extern int no_bond_direction;

//extern double dep_x, dep_y, dep_z;
						

extern Atom_Color atom_color[];

//extern double a_pn[3];

extern double initial_logtime;

extern int analysis_type;

extern double logtime_multiplier;

//extern double vacancy_density; //can always add back in

extern double overpotential_ramp_rate;

//extern int ncsk;

extern int total_volume_dissolved;