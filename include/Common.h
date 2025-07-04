#ifndef COMMON_H
#define COMMON_H

#include "Geometry.h"
#include <stdio.h>
#include <stdbool.h>

struct SimulationEnv
{
    unsigned long long int max_atoms;
    unsigned long long int max_transitions;

    int simulation_type;
    bool evaporation_flag;

    int zixshift, ziyshift, zizshift;
    int ssxshift, ssyshift, sszshift;
    int zsh, ysh, xsh; // shifts

    // [ ]: cartesian, units of nearest-neighbor spacing?
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
    double *substrate_compotition;
    int num_nn_levels;
    int num_elements;
    int num_element_combos; // C(num_elements,2)
    int num_neighbor_types;

    double normal_x, normal_y, normal_z; // XXX: likely vistigal

    // the lattice coorinate limits in order to accomodate the "simulation box"
    int simbox_limits_lat[3][2];

    double nnE[6];
    double nnnE[6];
    // double **nnE; // nnE[num_nn_levels][num_element_combos]

    
    bool solubility[3]; //all elements cannot dissolve by default
    
    double temperature;
    double overpotential;

    int dissolution; // flag for whether dissolution events can occur
    char atom_names[3][3];

};

struct SimulationState
{
    Bond bond[MAXIMUM_NUMBER_OF_BONDS];

    // array of Atom structs, initialized in ____
    // contains all atoms in the simulation, in the first atom_cnt indices
    Atom **atom_arr;
    int atom_cnt;

    // array of Rate structs, initialized in ____
    // contains all the unique rate constants and number of transitions with that rate const
    Rate *rate_arr;
    int rate_cnt;

    // array of Transition structs, initialized in ____
    // the list of possible transitions with atom and transition direction
    // all transitions of same Rate are next to each other
    // they live between rate.transition_start_idx and
    Transition **transition_arr;
    int transition_cnt;

    // array of zones, to help find atoms based on position
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
    int total_atoms_dissolved;
    
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
