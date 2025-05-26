#ifndef COMMON_H
#define COMMON_H

#include "Geometry.h"
#include <stdio.h>
#include <stdbool.h>


// Simulation.h
// extern Atom *g_atom_arr[MAXIMUM_NUMBER_OF_ATOMS];
// extern int g_atom_cnt;
// extern Rate g_rate_arr[MAXIMUM_NUMBER_OF_ACTIVATION_BARRIERS];
// extern int g_rate_cnt;
// extern Transition *g_transition_arr[MAXIMUM_NUMBER_OF_CONCURRENT_TRANSITIONS];
// extern int g_transition_cnt;
// extern Trans_Prob g_transition_probability;
// extern Bond g_bond[MAXIMUM_NUMBER_OF_BONDS];
// extern int g_simulation_type;
extern int g_num_sims; // XXX
extern Atom g_temp_atom;
// extern char g_atom_names[3][3];
extern double g_default_color[3];
// extern bool g_evaporation_flag;
extern char g_coordinate_log_prefix[256];
// extern unsigned long ss->final_iteration; // max number of iterations
// extern double ss->run_stime; // simulation max runtime default (in seconds)
// extern bool ss->simulation_should_kill_itself;
// extern double ss->elapsed_stime;
// extern int ss->sim_end_type;
extern int g_analysis_type;
extern double g_log_interval;
extern double g_next_log_checkpoint;
extern double* g_log_list;
extern int g_log_list_len;
// extern double ss->frequency_sum;
// extern double ss->overpotential;
// extern double g_nnE[6];
// extern double g_nnnE[6];
// extern bool g_solubility[3]; //all elements cannot dissolve by default
// extern double ss->temperature;
// extern int g_dissolution;
extern int g_final_config_neighbor_cnt;
extern int g_intial_config_neighbor_cnt;
extern int g_lastxt, g_lastyt, g_lastzt;
extern double g_sum_of_rate_populations;
extern double g_current_probability;

// Simulation_Aux.h
// extern double ss->total_internal_energy;
// extern int g_zixshift, g_ziyshift, g_zizshift;
// extern int g_ssxshift, g_ssyshift, g_sszshift;
// extern int g_zsh, g_ysh, g_xsh;
// extern int g_ssx, g_ssy, g_ssz;
// extern double g_ssr;
// extern int g_zix, g_ziy, g_ziz;
// extern int g_lattice_type;
// extern int g_max_neighbors;
// extern int g_sheet_thickness;
// extern int g_cluster_radius;
// extern char g_atoms_filename[256];
// extern Zone g_zone_arr[ZONES_IN_X][ZONES_IN_Y][ZONES_IN_Z];
// extern double g_initial_overpotential;
// extern double g_overpotential_ramp_rate;
// extern double g_max_overpotential;
// extern double g_substrate_percent_a;
// extern double g_substrate_percent_b;
// extern int ss->total_volume_dissolved;
// extern double g_normal_x, g_normal_y, g_normal_z; // XXX: likely vistigal
// extern int g_sblimits_lat[3][2];

extern double lhs[6];
extern double normal_lat[6][3];
extern int translation_vector[6][3];

// FileIO.h
extern char g_console_outstring[512];
extern FILE *g_view_save_file;
extern FILE *g_temp_log;
extern char g_command_string[1024];
extern char g_outFile[260];
extern char g_default_extension[];
extern FILE *g_sim_log_file;


struct SimulationEnv
{
    int simulation_type;
    bool evaporation_flag;

    int zixshift, ziyshift, zizshift;
    int ssxshift, ssyshift, sszshift;
    int zsh, ysh, xsh;

    // [ ]: what are the units for this? how does it relate to atomic spacing?
    int ssx, ssy, ssz; // system size x, y, z
    double ssr;
    int zix, ziy, ziz;

    int lattice_type;
    int max_neighbors;

    int sheet_thickness;
    int cluster_radius;
    char atoms_filename[256];

    double initial_overpotential;
    double overpotential_ramp_rate;
    double max_overpotential;

    double substrate_percent_a;
    double substrate_percent_b;

    double normal_x, normal_y, normal_z; // XXX: likely vistigal

    int simbox_limits_lat[3][2];

    double nnE[6];
    double nnnE[6];
    
    bool solubility[3]; //all elements cannot dissolve by default
    
    // double temperature;
    // double ss->overpotential;

    int dissolution; // flag for whether dissolution events can occur
    char atom_names[3][3];

};

struct SimulationState
{
    Bond bond[MAXIMUM_NUMBER_OF_BONDS];

    Atom *atom_arr[MAXIMUM_NUMBER_OF_ATOMS];
    int atom_cnt;

    // ENHANCE: malloc?
    Rate rate_arr[MAXIMUM_NUMBER_OF_ACTIVATION_BARRIERS];
    int rate_cnt;

    Transition* transition_arr[MAXIMUM_NUMBER_OF_CONCURRENT_TRANSITIONS];
    int transition_cnt;
    
    Zone zone_arr[ZONES_IN_X][ZONES_IN_Y][ZONES_IN_Z]; // TODO: change to malloc
    Trans_Prob transition_probability;

    unsigned long final_iteration; // max number of iterations
    double run_stime; // simulation max runtime default (in seconds)
    bool simulation_should_kill_itself;
    double elapsed_stime;
    int sim_end_type;

    double frequency_sum;
    double total_internal_energy;
    
    double temperature;
    double overpotential;
    int total_volume_dissolved;
    
};

struct LoggingState {
    // char console_outstring[512];
    // // char outFile[260];
    // char default_extension[12];
    // FILE *sim_log_file;
    // char position_log_prefix[256];

    // int analysis_type;
    // double log_interval;
    // double next_log_checkpoint;
    // double* log_list;
    // int log_list_len;

};

#endif // COMMON_H
