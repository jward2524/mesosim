#ifndef MESOSIM_H
#define MESOSIM_H

#include "Common.h"
#include <time.h>

// void initialize_symmetry_elements(void);
void initialize_lattice_geometry(struct SimulationEnv* sim_env);
void write_backlog(FILE*, FILE*);

#endif // MESOSIM_H
