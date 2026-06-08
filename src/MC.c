#include "MC.h"
#include "Atoms.h"
#include "Checkpoint.h"
#include "ErrorM.h"
#include "FileIO.h"
#include "Random.h"
#include "Simulation.h"
#include "Utils.h"
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// static const double FABS_TOL = 1e-6;

static double calculate_new_energy(const long atom_idx, const int offset_idx,
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
    // TODO: forbid analysis_types that aren't iteration

    int old_u, old_v, old_w;
    int new_u, new_v, new_w;
    int coord;
    bool simulation_run = true;
    while (simulation_run) {
        ss->iter++; // simulation is on iteration 1
        if (ss->iter % (unsigned long)ss->atom_cnt == 0) {
            ss->mmc_step++; // "" mmc_step 1
        }
        // flag if this iteration is the last iteration
        // based on iter counter so that final mmc_step can be completed
        // done here to keep close to increment
        simulation_run = ss->iter < (ss->final_iteration * (unsigned long)ss->atom_cnt);

        // select transition
        unsigned long long limit =
            ULLONG_MAX - (ULLONG_MAX % (unsigned long long)ss->transition_cnt);
        unsigned long long rand0;
        do {
            rand0 = ran(&se->rand_state);
        } while (rand0 >= limit);
        long int selected_idx = (long)(rand0 % (unsigned long long)ss->transition_cnt);
        Transition *selected_transition = ss->transition_arr[selected_idx];

        // calculate change in energy of system if transition were performed
        long transitioning_atom_idx = selected_transition->atom_idx;
        int offset_idx = selected_transition->offset_idx;

        double possible_energy = calculate_new_energy(transitioning_atom_idx, offset_idx, ss, se);
        double current_energy = ss->atom_arr[transitioning_atom_idx]->energy;

        // energy per atom contains half the sum of bond energies
        // (atoms on other end of bonds contain the other half)
        // multiply both by two to fully break and (re)form all bonds
        double deltaE = 2 * (possible_energy - current_energy);

        int perform_flag = 0;
        // accept transition if deltaE > 0
        // > 0 because a bond is positive energy, and more bonds is more stable
        // sign is opposite of convention; total bond energy is being maximized
        // ENHANCE: fix sign of energy to follow conventions
        if (deltaE >= 0) {
            perform_flag = 1;
        }
        // otherwise, accept based on Boltzmann probability
        else {
            double boltzmann_prob = exp(deltaE / (kBoltz * ss->temperature));
            double random = dran(&se->rand_state);
            if (random <= boltzmann_prob) {
                perform_flag = 1;
            }
        }

        old_u = ss->atom_arr[transitioning_atom_idx]->lattice[0];
        old_v = ss->atom_arr[transitioning_atom_idx]->lattice[1];
        old_w = ss->atom_arr[transitioning_atom_idx]->lattice[2];

        new_u = old_u + se->transition_vectors[offset_idx].dx;
        new_v = old_v + se->transition_vectors[offset_idx].dy;
        new_w = old_w + se->transition_vectors[offset_idx].dz;

        if (perform_flag) {
            // perform transition
            // duplicating what is in perform_kmc

            coord = get_coordination(transitioning_atom_idx, ss, se);

            adjust_pbc(&new_u, &new_v, &new_w, se);

            unsigned char atype = ss->atom_arr[transitioning_atom_idx]->type;

            // moves atom?
            remove_atom(transitioning_atom_idx, ss, se);
            long transitioned_atom_idx = add_atom(new_u, new_v, new_w, atype, NORMAL, ss, se);

            // update system
            update_outdated_transitions(old_u, old_v, old_w, transitioned_atom_idx, ss, se);
        }

        if (ls->verbose) {
            if (ss->iter % 10000 == 0) {
                // printf("Iteration %ld, energy %le\n", ss->iter, ss->total_internal_energy);
            }
            if (ss->iter % (unsigned long)ss->atom_cnt == 0) {
                printf("MMC %lu, iteration %ld, energy %le\n", ss->mmc_step, ss->iter,
                       ss->total_internal_energy);
            }
        }

        StepData step_data = {
            .flavor = se->flavor,
            .iter = ss->iter,
            .sys_energy = ss->total_internal_energy,
            .uvw1 = {old_u, old_v, old_w},
            .uvw2 = {new_u, new_v, new_w},
            .coord = coord,
            .mc =
                {
                    .deltaE = deltaE,
                    .performed = perform_flag,
                },
        };

        output_on_schedule(&step_data, ss, se, ls);
        checkpoint_on_schedule(ss, se, ls);
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
static double calculate_new_energy(const long atom_idx, const int offset_idx,
                                   struct SimulationState *ss, struct SimulationEnv *se)
{
    // a mix of get_final_configuration and refresh_transitions

    // atom position after jump offset_idx
    int new_u = ss->atom_arr[atom_idx]->lattice[0] + se->transition_vectors[offset_idx].dx;
    int new_v = ss->atom_arr[atom_idx]->lattice[1] + se->transition_vectors[offset_idx].dy;
    int new_w = ss->atom_arr[atom_idx]->lattice[2] + se->transition_vectors[offset_idx].dz;

    adjust_pbc(&new_u, &new_v, &new_w, se);

    double energy = 0;
    for (int i = 0; i < se->num_energy_contributors; ++i) {
        // if considering the direction where the atom would have come from,
        // ignore it since it would be an empty site after the move
        if (i == se->opposite_tvectors[offset_idx]) {
            continue;
        }

        // location of neighbor
        long neighbor_idx = atom_at_offset(new_u, new_v, new_w, i, ss->atom_arr, ss->zone_arr, se);

        if (neighbor_idx >= 0) {
            int atom_type = ss->atom_arr[atom_idx]->type;
            int neighbor_type = ss->atom_arr[neighbor_idx]->type;
            int env_idx = get_env_index(i, atom_type, neighbor_type, se);
            energy += se->nn_energy[env_idx];
        }
    }

    return energy / 2.;
}
