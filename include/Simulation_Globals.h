bool evaporation_flag = true;

FILE *sim_log_file = NULL;

char coordinate_log_prefix[256] = "default_simulation_analysis.dat";

//int number_of_simulation_runs = 2000;
		
// [ ]: what are these?
int zixshift, ziyshift, zizshift;
int ssxshift, ssyshift, sszshift;
int zsh, ysh, xsh;							// shifts

// [ ]: what are the units for this? how does it relate to atomic spacing?
int ssx = DSIMSIZE, ssy = DSIMSIZE, ssz = DSIMSIZE;					// system size x, y, z
double ssr;
int zix = TTS, ziy = TTS, ziz = TTS;							// zones in x, y, z

double run_time = 1.e8; //default (in seconds)
double data_time_interval = 0.1;
double time_interval_end;


int lattice_type = FCC;
int number_of_possible_neighbors = 12;

int sheet_thickness = -1;
int cluster_radius = -1;
char atoms_filename[256] = "";

// [ ]: what is this for?
Zone zone[ZONES_IN_X][ZONES_IN_Y][ZONES_IN_Z];

int number_rates;
int total_current_transitions;
double sum_of_frequencies;

double overpotential = 0.0;
double initialoverpotential = DEFAULT_OVERPOTENTIAL;
double overpotentialramprate = 0.0;
double maxoverpotential = DEFAULT_OVERPOTENTIAL;

double nnE[6] = {
				DEFAULT_BOND_ENERGY_AA,
				DEFAULT_BOND_ENERGY_AB,
				DEFAULT_BOND_ENERGY_AC,
				DEFAULT_BOND_ENERGY_BB,
				DEFAULT_BOND_ENERGY_BC,
				DEFAULT_BOND_ENERGY_CC};

double nnnE[6] = {0., 0., 0., 0., 0., 0.};

bool solubility[3] = {false, false, false}; //all elements cannot dissolve by default

/*double deposition_rate_of_a = DEFAULT_DEPOSITION_RATE_OF_A;
double deposition_rate_of_b = DEFAULT_DEPOSITION_RATE_OF_B;
double deposition_rate_of_c = DEFAULT_DEPOSITION_RATE_OF_C;*/

double temperature = DEFAULT_TEMPERATURE;


double substrate_percent_a = DEFAULT_COMPOSITION_A;
double substrate_percent_b = DEFAULT_COMPOSITION_B;
int dissolution = DISSOLUTION;

int number_final_configuration_neighbors;
int number_intial_configuration_neighbors;

// ENHANCE: malloc?
Rate rate[MAXIMUM_NUMBER_OF_ACTIVATION_BARRIERS];

Transition_List *transition_list[MAXIMUM_NUMBER_OF_CONCURRENT_TRANSITIONS];
//Transition_List transition_list[MAXIMUM_NUMBER_OF_CONCURRENT_TRANSITIONS];

Trans_Prob transition_probability;


long int final_iteration = 1e9;
double lastxt, lastyt, lastzt;

//bool simulation_is_going = false; //probably don't need this

double sum_of_rate_populations;
double current_probability;

//int no_bond_direction;

//double dep_x, dep_y, dep_z;
						
Atom_Color atom_color[10];

//double a_pn[3];

double initial_logtime = 1.0e-4;

int analysis_type = REGULAR_TIME_INTERVALS;

double logtime_multiplier;

//double vacancy_density = 0.01; //can always add back in

double overpotential_ramp_rate = 0.0;

//int ncsk = 0;

int total_volume_dissolved;