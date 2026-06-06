#ifndef STATE_H
#define STATE_H

#include "Geometry.h"
#include "Random.h"
#include <stdbool.h> // bool
#include <stdio.h>   // FILE

#define FLAVOR_UNDEFINED 0
#define FLAVOR_KMC 1
#define FLAVOR_MC 2

// general note:
// uvw refer to coordinates in lattice vector space
// xyz refer to coordinates in cartesian (orthogonal basis) space

struct SimulationConfig {

    /* Not in simulation structs */
    // initial geometry
    int geometry;
    int geometry_param;
    char atoms_filename[256];
    int lattice_type;

    /* SimulationState */
    int sim_end_type;
    double run_stime;
    unsigned long final_iteration; // max number of iterations
    double temperature;
    unsigned rand_seed;
    RandomState rand_state;

    /* SimulationEnv */
    int system_size_x, system_size_y, system_size_z;

    double initial_overpotential;
    double overpotential_ramp_rate;
    double max_overpotential;

    double *substrate_composition;

    int dissolution;
    bool *is_soluble;
    int atom_names_cnt;
    char **atom_names;

    int num_nn_levels;
    int num_nn_types;
    int num_elements;
    int num_bond_types;
    double *nn_energy;

    unsigned flavor;

    // output formats
};

// TODO: break out user-defined parameters (from input file)
// to separate parameters that set derived parameters from derived parameters
// variables for simulation conditions
// (generally) initialized and then never edited again, only read
struct SimulationEnv {

    /* --------- Simulation bounds --------- */
    unsigned flavor;
    RandomState rand_state;

    long int max_atoms;
    long int max_transitions;
    long int max_rates;

    // TODO: to input?
    size_t zone_count_u, zone_count_v, zone_count_w;

    double overpotential_ramp_rate;
    double max_overpotential;

    /* --------- Simulation space --------- */
    // [ ]: cartesian, units of nearest-neighbor spacing?
    int system_size_x, system_size_y, system_size_z;

    // shifts
    int zixshift, ziyshift, zizshift;
    int ssxshift, ssyshift, sszshift;
    int zsh, ysh, xsh;

    // the lattice coorinate limits in order to accomodate the "simulation box"
    int simbox_limits_lat[3][2];
    int lat_range[3];

    // cartesian simulation cell vectors
    double simbox_vectors_cart[3][3]; // [[u1 u2 u3], [v1 v2 v3], [w1 w2 w3]]
    double simbox_origin_cart[3];

    // relevant vectors to neighbors
    int lattice_type;
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

    /* --------- Atom properties --------- */
    int dissolution;  // flag for whether dissolution events can occur
    bool *is_soluble; // [num_elements] (implemented as always [8])

    char **atom_names; // [num_elements][BUFFER_SIZE or 3 or something]
    int atom_names_cnt;
    // TODO: remove this, can just use num_elements

    /* --------- Atomic environments --------- */
    double *substrate_composition;
    // TODO: change most of these to size_t, if they are used in mallocs
    int num_nn_levels;
    int num_elements;
    int num_bond_types; // Cr(num_elements,2), combination with replacement [two atoms per bond]
    // int num_neighbor_types; // XXX
    int num_nn_types; // number of distinct bond types, se->num_nn_levels * se->num_bond_types
    int *atoms_per_nn_level; // number of atoms per nn level, [num_nn_levels]

    // TODO: change this (and get_env_index) to be a 3D array nne[shell][i][j] where i and j are
    // atom types
    double *nn_energy; // nnE[env_idx]
};

// variables that describe the simulation state
// change on every iteration
struct SimulationState {
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
    Zone ***zone_arr;
    TransProb transition_probability;

    unsigned long iter;
    unsigned long mmc_steps;
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

// TODO: add logarithmic intervals
typedef enum {
    OUTPUT_SCHEDULE_NONE = 0,
    OUTPUT_SCHEDULE_INTERVAL_ITERATION,
    OUTPUT_SCHEDULE_INTERVAL_TIME,
    OUTPUT_SCHEDULE_LIST_ITERATION,
    OUTPUT_SCHEDULE_LIST_TIME
} OutputScheduleMode;

typedef struct {
    OutputScheduleMode mode; // schedule mode
    double interval;         // interval value (if mode is interval)
    double *list;            // list of output points (if mode is list)
    int list_len;            // number of entries in list
    int list_idx;            // current index in list
    double next_checkpoint;
    int frame_num;
} OutputSchedule;

struct CsvLsView {
    int precision;
};
// returns malloced pointer to string
typedef const char *(*CsvFieldFuncPtr)(const struct SimulationState *ss,
                                       const struct CsvLsView *view);

typedef enum { OUTPUT_FORMAT_CSV = 1, OUTPUT_FORMAT_XYZ, OUTPUT_FORMAT_STEPS_CSV } OutputFormatType;

// ENHANCE: can define each struct in the union seperately for better type control in functions?
typedef struct {
    int type;
    bool is_active;
    union {
        struct {
            char filename[256];
            FILE *file;
            OutputSchedule schedule;
            CsvFieldFuncPtr *field_funcs;
            char **field_names;
            int field_count;
            int frame_num;
        } csv;
        struct {
            OutputSchedule schedule;
            char prefix[256];
            char suffix[256];
            int frame_num;
            bool stripped; // only output under-coordinated atoms
        } xyz;
        struct {
            char filename[256];
            FILE *file;
            bool with_coordination;
        } steps;
    };
} OutputFormat;

typedef struct {
    char filename[256];
    unsigned long interval;
    unsigned long next_checkpoint;
    int frame_num;
} CheckpointFormat;

// variables that describe state of logging
// files and when to log
struct LoggingState {
    FILE *sim_log;
    int framenum; // counter/id for number of outputs / output files
    int verbose;
    unsigned long verbose_interval; // interval for printing verbose output to console

    // precision of incremented doubles
    int increment_precision;
    int stime_precision;
    int overpot_precision;

    CheckpointFormat checkpoint;

    OutputFormat *out_formats;
    int out_formats_cnt;
};

#endif // STATE_H
