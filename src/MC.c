#include "Atoms.h"
#include "Simulation.h"
#include "State.h"
#include "Utils.h"
#include <stdlib.h>
#include <math.h>

static double calculate_new_energy(const int atom_idx, const int offset_idx,
                                   struct SimulationState *ss, struct SimulationEnv *se);

/**
 * @brief
 *
 * @param ss
 * @param se
 * @param ls
 * @return unsigned long
 */
unsigned long perform_metropolis_mc(struct SimulationState *ss, struct SimulationEnv *se,
                                    struct LoggingState *ls)
{
    // Metropolis MC steps per particle
    unsigned long int mmc_steps = 0;
    ss->iter = 0;

    while (mmc_steps <= ss->final_iteration) {
        // select transition
        // use multiple calls to rand() to allow selection from all possible transitions, not
        // only limited to RAND_MAX
        long long int mult = LLONG_MAX / (long long int)RAND_MAX;
        long long int llrand = 0;
        for (long long int i = 0; i < mult; i++) {
            llrand += rand();
        }

        long long int selected_idx = ((double)llrand / (RAND_MAX * mult)) * ss->transition_cnt;
        Transition *selected_transition = ss->transition_arr[selected_idx];

        // calculate change in energy of system if transition were performed
        int transitioning_atom_idx = selected_transition->atom_idx;
        int offset_idx = selected_transition->offset_idx;

        double possible_energy = calculate_new_energy(transitioning_atom_idx, offset_idx, ss, se);
        double current_energy = ss->atom_arr[transitioning_atom_idx]->energy;

        // energy per atom contains half the sum of bond energies
        // (atoms on other end of bonds contain the other half)
        // multiply both by two to fully break and (re)form all bonds
        double deltaE = 2 * (possible_energy - current_energy);

        int perform_flag = 0;
        // accept transition if deltaE < 0
        if (deltaE <= 0) {
            perform_flag = 1;
        }
        // otherwise, accept based on Boltzmann probability
        else {
            double boltzmann_prob = exp(-1 * deltaE / (kBoltz * ss->temperature));
            double random = drand();
            if (random <= boltzmann_prob) {
                perform_flag = 1;
            }
        }

        if (perform_flag) {
            // perform transition
            // duplicating what is in perform_simulation
            int old_x = ss->atom_arr[transitioning_atom_idx]->lattice[0];
            int old_y = ss->atom_arr[transitioning_atom_idx]->lattice[1];
            int old_z = ss->atom_arr[transitioning_atom_idx]->lattice[2];

            int new_x = old_x + se->transition_vectors[offset_idx].dx;
            int new_y = old_y + se->transition_vectors[offset_idx].dy;
            int new_z = old_z + se->transition_vectors[offset_idx].dz;

            adjust_pbc(&new_x, &new_y, &new_z, se);

            int atype = ss->atom_arr[transitioning_atom_idx]->type;

            // moves atom?
            remove_atom(transitioning_atom_idx, ss, se);
            int transitioned_atom_idx = add_atom(new_x, new_y, new_z, atype, NORMAL, ss, se);

            // update system
            update_outdated_transitions(old_x, old_y, old_z, transitioned_atom_idx, ss, se);
        }
        ss->iter++;

        if (ss->iter % ss->atom_cnt == 0) {
            mmc_steps++;
            // log
        }
    }
    return 0;
}

/**
 * @brief calculates the energy of a hypotetical move
 *
 * @param atom_idx atom that will be moving
 * @param offset_idx index of the direction of movement (index of se->transition_vector)
 * @param ss SimulationState
 * @param se SimulationEnv
 * @return double energy contribution of the atom after the move
 */
static double calculate_new_energy(const int atom_idx, const int offset_idx,
                                   struct SimulationState *ss, struct SimulationEnv *se)
{
    // a mix of get_final_configuration and refresh_transitions

    // atom position after jump offset_idx
    int new_x = ss->atom_arr[atom_idx]->lattice[0] + se->transition_vectors[offset_idx].dx;
    int new_y = ss->atom_arr[atom_idx]->lattice[1] + se->transition_vectors[offset_idx].dy;
    int new_z = ss->atom_arr[atom_idx]->lattice[2] + se->transition_vectors[offset_idx].dz;

    adjust_pbc(&new_x, &new_y, &new_z, se);

    double energy = 0;
    for (int i = 0; i < se->num_energy_contributors; ++i) {
        // if considering the direction where the atom would have come from,
        // ignore it since it would be an empty site after the move
        if (i == se->opposite_tvectors[offset_idx]) {
            continue;
        }

        // location of neighbor
        int neighbor_idx = atom_at_offset(new_x, new_y, new_z, i, ss->atom_arr, ss->zone_arr, se);

        int atom_type = ss->atom_arr[atom_idx]->type;
        int neighbor_type = ss->atom_arr[neighbor_idx]->type;
        int env_idx = get_env_index(i, atom_type, neighbor_type, se);
        energy += se->nn_energy[env_idx];
    }

    return energy / 2.;
}
