unsigned long perform_simulation(void); //modified!
void check_system(void);
void take_off_transition_list(int at, int vc);
int is_on_transition_list(double r);
int create_new_transition(double r);
void add_to_transition_list(int tl, int at, int vc);
int refresh_transitions(int at);
void compute_transition_array(void);
int calculate_surf_diffusion_rate(	int initial_configuration[],
    int final_configuration[],
    int number_of_neighbors,
    int atom_type,
    double nnE[6],
    double temperature,
    double overpotential,
    double *rate
);

int	calculate_evaporation_rate( int initial_configuration[],
    int number_of_neighbors,
    int atom_type,
    double nnE[6],
    double temperature,
    double overpotential,
    double *rate
);

extern int simulation_type;
// bool simulation_initialized;
extern Atom *atom_arr[MAXIMUM_NUMBER_OF_ATOMS];
extern Bond bond[MAXIMUM_NUMBER_OF_BONDS];
extern int atom_cnt;
extern double overpotential;
extern int num_sims;
extern double elapsed_time;
extern Atom temp_atom;
extern double default_color[3];
extern bool simulation_should_kill_itself;
extern char atom_names[3][3];

extern double default_color[3];

extern int simulation_type;
extern int atom_cnt;
extern Atom temp_atom;
extern bool simulation_should_kill_itself;
extern Atom* atom_arr[];
extern double elapsed_time;
 
extern bool evaporation_flag;
extern char coordinate_log_prefix[256];
 
extern double run_time; //default (in seconds)
extern double data_time_interval;
extern double time_interval_end;
 
extern int rate_cnt;
extern int transition_cnt;
extern double frequency_sum;

extern double overpotential;
 
extern double nnE[6];
extern double nnnE[6];
 
extern bool solubility[3]; //all elements cannot dissolve by default

extern double temperature;

extern int dissolution;
 
extern int final_config_neighbor_cnt;
extern int intial_config_neighbor_cnt;
 
// ENHANCE: malloc?
extern Rate rate_arr[MAXIMUM_NUMBER_OF_ACTIVATION_BARRIERS];
 
extern Transition *transition_arr[MAXIMUM_NUMBER_OF_CONCURRENT_TRANSITIONS];
//Transition transition_arr[MAXIMUM_NUMBER_OF_CONCURRENT_TRANSITIONS];
 
extern Trans_Prob transition_probability;
 
 
extern long int final_iteration;
extern int lastxt, lastyt, lastzt;
 
//bool simulation_is_going; //probably don't need this
 
extern double sum_of_rate_populations;
extern double current_probability;