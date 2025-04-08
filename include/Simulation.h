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
    double *rate);

int	calculate_evaporation_rate( int initial_configuration[],
    int number_of_neighbors,
    int atom_type,
    double nnE[6],
    double temperature,
    double overpotential,
    double *rate);
