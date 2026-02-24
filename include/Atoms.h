#ifndef ATOMS_H
#define ATOMS_H

#include "State.h"

long int add_atom(int u, int v, int w, unsigned char type, int special, struct SimulationState *ss,
                  struct SimulationEnv *se);
long atom_at(int u, int v, int w, Atom **atom_arr,
             Zone zone_arr[ZONES_IN_X][ZONES_IN_Y][ZONES_IN_Z], struct SimulationEnv *se);
long atom_at_offset(int u, int v, int w, int offset, Atom **atom_arr,
                    Zone zone_arr[ZONES_IN_X][ZONES_IN_Y][ZONES_IN_Z], struct SimulationEnv *se);
void remove_atom(long int at, struct SimulationState *ss, struct SimulationEnv *se);
void move_atom(long int ia, long int fa, Atom **atom_arr,
               Zone zone_arr[ZONES_IN_X][ZONES_IN_Y][ZONES_IN_Z], Transition **transition_arr,
               struct SimulationEnv *se);

// void make_buried_atoms_real(void);
void bury_atom(int at, int *pos); // no longer relevant

// rotation and molecule orientation
// void rotmat(Atom* atm[], int na, double rtmat[3][3]);
void organize(Atom **atom_arr, int atom_cnt,
              double primitive_basis[3][3]); // keep for now but need to retool
void orthomol(Atom **atom_arr, int atom_cnt, double basis[3][3]);
void centerg(Atom **atom_arr, int atom_cnt,
             double centroid[3]); // not really relevant now (visualization)

// general atom and bond handling routines
void copy_atom(long int i, long int j, Atom **atom_arr);
void create_default_atom(long int atom_idx, Atom **atom_arr,
                         struct SimulationEnv *se); // can modify this to remove things like color?

void kill_atom(long atom_number, struct SimulationState *ss, struct SimulationEnv *se);

int get_initial_configuration(long atom_index, int num_transition_vectors, Atom **atom_arr,
                              int initial_config[]);
int get_final_configuration(long at, int offset_idx, struct SimulationState *ss,
                            struct SimulationEnv *se, int final_config[]);
int get_coordination(long int atom_idx, struct SimulationState *ss, struct SimulationEnv *se);

#endif // ATOMS_H
