#ifndef COMMON_H
#define COMMON_H

#include "Geometry.h"
#include <stdbool.h> // bool
#include <stdio.h>   // FILE

#define FLAVOR_UNDEFINED 0
#define FLAVOR_KMC 1
#define FLAVOR_MC 2

// general note:
// uvw refer to coordinates in lattice vector space
// xyz refer to coordinates in cartesian (orthogonal basis) space

// variables for simulation conditions
// (generally) initialized and then never edited again, only read
struct SimulationEnv {

    unsigned flavor;
    unsigned rand_seed;

    long int max_atoms;
    long int max_transitions;
    long int max_rates;

    int geometry;
    bool evaporation_flag;

    int zixshift, ziyshift, zizshift;
    int ssxshift, ssyshift, sszshift;
    int zsh, ysh, xsh; // shifts

    // [ ]: cartesian, units of nearest-neighbor spacing?
    int system_size_x, system_size_y, system_size_z; // system size x, y, z
    double ssr;
    int zone_count_u, zone_count_v, zone_count_w;

    int lattice_type;

    int sheet_thickness;
    int cluster_radius;
    char atoms_filename[256];

    double initial_overpotential;
    double overpotential_ramp_rate;
    double max_overpotential;

    double *substrate_composition;
    // TODO: change most of these to size_t, if they are used in mallocs
    int num_nn_levels;
    int num_elements;
    int num_bond_types; // Cr(num_elements,2), combination with replacement [two atoms per bond]
    int num_neighbor_types;
    int num_nn_types; // number of distinct bond types, se->num_nn_levels * se->num_bond_types
    int *atoms_per_nn_level; // number of atoms per nn level, [num_nn_levels]

    double normal_x, normal_y, normal_z; // XXX: likely vistigal

    // the lattice coorinate limits in order to accomodate the "simulation box"
    int simbox_limits_lat[3][2];
    int lat_range[3];

    // cartesian simulation cell vectors
    double simbox_vectors_cart[3][3]; // [[u1 u2 u3], [v1 v2 v3], [w1 w2 w3]]
    double simbox_origin_cart[3];

    // TODO: change this (and get_env_index) to be a 3D array nne[shell][i][j] where i and j are atom types
    double *nn_energy; // nnE[env_idx]

    bool *is_soluble; // [num_elements] (implemented as always [8])

    int dissolution;   // flag for whether dissolution events can occur
    char **atom_names; // [num_elements][BUFFER_SIZE or 3 or something]
    int atom_names_cnt;
    // TODO: remove this, can just use num_elements

    /* symmetry related variables */
    double rmat[3][3];  // visualization?
    double centroid[3]; // coordinates for center of gravity

    // relevant vectors to neighbors
    LatticeVector *transition_vectors;
    int num_transition_vectors;  // TODO: change to unsigned char to match offset
    int num_energy_contributors; // based on num_nn_levels, directions must be hard-coded

    // index in transition_vectors that has the jump in the opposite direction in simulation;
    // opposite_tvectors[0]=11 means the opposite direction of se->transition_vectors[0] is
    // se->transition_vectors[11]
    int *opposite_tvectors;

    // primitive unit cell basis vectors + inverted; primitive_basis[*][0] = basis1,
    // primitive_basis[0][*] = x component of basises is also transformation matrix for [lattice to
    // cartesian coordinates] [cartesian to lattice coordinates] respectively
    // [[u1 v1 w1], [u2 v2 w2], [u3 v3 w3]]
    double primitive_basis[3][3];
    double invert_primitive_basis[3][3];

    // TODO: this is never used except for printing?
    double ucell_params[6]; // unit cell parameters; a b c alpha beta gamma
};

// variables that describe the simulation state
// change on every iteration
struct SimulationState {
    Bond bond[MAXIMUM_NUMBER_OF_BONDS];

    // array of Atom structs, initialized in ____
    // contains all atoms in the simulation, in the first atom_cnt indices
    Atom **atom_arr;
    long int atom_cnt;

    // array of Rate structs, initialized in ____
    // contains all the unique rate constants and number of transitions with that rate const
    Rate *rate_arr;
    long int rate_cnt;

    // array of Transition structs, initialized in ____
    // the list of possible transitions with atom and transition direction
    // all transitions of same Rate are next to each other
    // they live between rate.transition_start_idx and
    Transition **transition_arr;
    long int transition_cnt;

    // array of zones, to help find atoms based on position
    Zone zone_arr[ZONES_IN_X][ZONES_IN_Y][ZONES_IN_Z]; // TODO: change to malloc
    TransProb transition_probability;

    unsigned long iter;
    unsigned long final_iteration; // max number of iterations
    double run_stime;              // simulation max runtime default (in seconds)
    bool simulation_should_kill_itself;
    double elapsed_stime; // current simulation time (in seconds)
    int sim_end_type;

    double frequency_sum;

    // sum of all bond energies in system (no double-counting)
    // currently: incremental modification of this value throughout simulation and storage as a
    // floating-point number in refresh_transitions and remove_atom means it loses some accuracy
    // over iterations return to sum over all atom energies (during log checkpoint) if want to
    // regain some precision
    double total_internal_energy;

    double temperature;
    double overpotential;
    int total_atoms_dissolved;
};

typedef enum {
    OUTPUT_SCHEDULE_NONE = 0,
    OUTPUT_SCHEDULE_INTERVAL_ITERATION,
    OUTPUT_SCHEDULE_INTERVAL_TIME,
    OUTPUT_SCHEDULE_LIST_ITERATION,
    OUTPUT_SCHEDULE_LIST_TIME
} OutputScheduleMode;

typedef struct {
    OutputScheduleMode mode;           // schedule mode
    double interval;                   // interval value (if mode is interval)
    double *list;                      // list of output points (if mode is list)
    int list_len;                      // number of entries in list
} OutputSchedule;

// variables that describe state of logging
// files and when to log
struct LoggingState {
    char console_outstring[512];
    // char outFile[260];
    char default_extension[12];
    FILE *sim_log;
    char position_log_prefix[256];

    bool output_xyz;

    bool output_state_csv;
    FILE *state_csv;

    bool output_iter_csv;
    FILE *iter_csv;

    int analysis_type;
    double log_interval; // interval between log checkpoints, based on analysis_type?
    double next_log_checkpoint;
    double *log_list; // list of log checkpoints
    int log_list_len; // length of log_list
    int framenum;     // counter/id for number of outputs / output files

    int verbose;

    // CSV output configuration
    char csv_filename[256];
    bool csv_filename_provided;
    char **csv_fields;
    int csv_field_count;
    OutputSchedule csv_schedule;

    // XYZ output configuration
    char xyz_prefix[256];
    bool xyz_prefix_provided;
    OutputSchedule xyz_schedule;
    int xyz_framenum;
};

#endif // COMMON_H
