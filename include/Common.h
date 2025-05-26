#ifndef COMMON_H
#define COMMON_H

#include "Geometry.h"
#include <stdio.h>
#include <stdbool.h>

struct SimulationEnv
{
    int simulation_type;
    bool evaporation_flag;

    int zixshift, ziyshift, zizshift;
    int ssxshift, ssyshift, sszshift;
    int zsh, ysh, xsh; // shifts

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
    
    double temperature;
    double overpotential;

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
    TransProb transition_probability;

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
    char console_outstring[512];
    // char outFile[260];
    char default_extension[12];
    FILE *sim_log_file;
    char position_log_prefix[256];

    int analysis_type;
    double log_interval;
    double next_log_checkpoint;
    double* log_list;
    int log_list_len;

};

#endif // COMMON_H
