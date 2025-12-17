#ifndef ATOMS_H
#define ATOMS_H

#include "State.h"

int add_atom(int u, int v, int w, int type, int special, struct SimulationState *ss,
             struct SimulationEnv *se);
int atom_at(int u, int v, int w, Atom **atom_arr, Zone zone_arr[ZONES_IN_X][ZONES_IN_Y][ZONES_IN_Z],
            struct SimulationEnv *se);
void remove_atom(int at, struct SimulationState *ss, struct SimulationEnv *se);
void move_atom(int ia, int fa, Atom **atom_arr, Zone zone_arr[ZONES_IN_X][ZONES_IN_Y][ZONES_IN_Z],
               Transition **transition_arr, struct SimulationEnv *se);

// void make_buried_atoms_real(void);

int random_reincarnate_atom(int u, int v, int w, int type,
                            int vc);                         // feels like it'll never be called
int reincarnate_atom(int u, int v, int w, int type, int vc); // will not be called
void bury_atom(int at, int *pos);                            // no longer relevant

// rotation and molecule orientation
// void rotmat(Atom* atm[], int na, double rtmat[3][3]);
void organize(Atom **atom_arr, int atom_cnt,
              double primitive_basis[3][3]); // keep for now but need to retool
void orthomol(Atom **atom_arr, int atom_cnt, double basis[3][3]);
void centerg(Atom **atom_arr, int atom_cnt,
             double centroid[3]); // not really relevant now (visualization)

// general atom and bond handling routines
void copy_atom(int i, int j, Atom **atom_arr);
void create_default_atom(int atom_idx, Atom **atom_arr,
                         struct SimulationEnv *se); // can modify this to remove things like color?

void kill_atom(int atom_number, struct SimulationState *ss, struct SimulationEnv *se);

int get_initial_configuration(int atom_index, int num_transition_vectors, Atom **atom_arr,
                              int initial_config[]);
int get_final_configuration(int at, int offset_idx, struct SimulationState *ss,
                            struct SimulationEnv *se, int final_config[]);

// TODO: into SimulationEnv
extern const LatticeVector BCC_OFFSET[8];
extern const LatticeVector FCC_OFFSET[12];
extern const LatticeVector SC_OFFSET[6];
extern const LatticeVector SC_OFFSET_2[12];
extern const LatticeVector FCC_OFFSET_2[6];
extern const LatticeVector BCC_OFFSET_2[6];

#endif // ATOMS_H
