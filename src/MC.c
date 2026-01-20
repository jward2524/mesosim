#include "Atoms.h"
#include "ErrorM.h"
#include "FileIO.h"
#include "Simulation.h"
#include "State.h"
#include "Utils.h"
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static const double FABS_TOL = 1e-6;

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
    // Metropolis MC steps per particle
    unsigned long int mmc_steps = 0;

    // TODO: was this set elsewhere?
    // ss->iter = 0;

    for (int i = 0; i < ss->atom_cnt; ++i) {
        refresh_transitions(i, ss, se);
    }

    // TODO: forbid analysis_types that aren't iteration
    char suffix[BUFFER_SIZE];
    snprintf(suffix, BUFFER_SIZE, "i0");

    // initial state
    output_log_file(ls->sim_log, ls->framenum, mmc_steps, ss->elapsed_stime, ss->temperature,
                    ss->overpotential, ss->atom_cnt, ss->total_internal_energy);
    if (ls->output_state_csv) {
        log_mc_state_csv(ls->state_csv, ls->framenum, mmc_steps, ss->iter,
                         ss->total_internal_energy);
    }
    if (ls->output_xyz) {
        write_xyz_file(ls->position_log_prefix, ls->framenum, suffix, ss, se);
    }
    ls->framenum++;

    int checkpoint_reached = 0;
    int old_u, old_v, old_w;
    int new_u, new_v, new_w;
    int simulation_run = 1;
    while (simulation_run) {
        ss->iter++;
        if (ss->iter % (unsigned long)ss->atom_cnt == 0) {
            mmc_steps++;
        }
        // flag if this iteration is the last iteration
        simulation_run = ss->iter < (ss->final_iteration * (unsigned long)ss->atom_cnt);

        // select transition
        // use multiple calls to rand() to allow selection from all possible transitions, not
        // only limited to RAND_MAX
        long long int extended_rand = 0;
        double drand0 = -1;
        int rand0 = -1;
        if (ss->transition_cnt > RAND_MAX) {
            int div = (int)floor((double)ss->transition_cnt / RAND_MAX);
            int rem = (int)(ss->transition_cnt % RAND_MAX);

            for (int i = 0; i < div; i++) {
                rand0 = rand();
                extended_rand += rand0;
            }
            drand0 = drand();
            extended_rand += (long long)nearbyint(rem * drand0);
        } else {
            drand0 = drand();
            extended_rand = (long long)nearbyint(drand0 * (double)ss->transition_cnt);
            if (extended_rand == ss->transition_cnt) {
                extended_rand -= 1;
            }
        }

        long long int selected_idx = extended_rand;
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
            double random = drand();
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
            // duplicating what is in perform_simulation

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
                printf("MMC %lu, iteration %ld, energy %le\n", mmc_steps, ss->iter,
                       ss->total_internal_energy);
            }
        }

        // csv log
        if (ls->output_iter_csv) {
            int uvw1[] = {old_u, old_v, old_w};
            int uvw2[] = {new_u, new_v, new_w};
            log_mc_iter(ls->iter_csv, ss->iter, ss->total_internal_energy, deltaE, perform_flag,
                        uvw1, uvw2);
        }

        if (ls->analysis_type == ITERATION_INTERVALS) {
            checkpoint_reached = fabs(ls->next_log_checkpoint - (double)mmc_steps) < FABS_TOL;
            if (!checkpoint_reached && (mmc_steps > ls->next_log_checkpoint)) {
                fprintf(stderr, "Iterations (%lu) exceeded log checkpoint (%lf) without noticing",
                        mmc_steps, ls->next_log_checkpoint);
                clean_and_error(1);
            }
        }

        if (checkpoint_reached) {
            if (ls->analysis_type == ITERATION_INTERVALS) {
                snprintf(suffix, BUFFER_SIZE, "i%lu", (unsigned long)ls->next_log_checkpoint);
                ls->next_log_checkpoint += ls->log_interval;
            }

            // record iteration in a file here
            printf("Writing file %d: MC steps = %lu\n", ls->framenum, mmc_steps);
            output_log_file(ls->sim_log, ls->framenum, mmc_steps, ss->elapsed_stime,
                            ss->temperature, ss->overpotential, ss->atom_cnt,
                            ss->total_internal_energy);
            if (ls->output_state_csv) {
                log_mc_state_csv(ls->state_csv, ls->framenum, mmc_steps, ss->iter,
                                 ss->total_internal_energy);
            }
            if (ls->output_xyz) {
                write_xyz_file(ls->position_log_prefix, ls->framenum, suffix, ss, se);
            }

            ++ls->framenum;
        }
    }

    // write elapsed_stime to mark finish
    output_log_file(ls->sim_log, ls->framenum, mmc_steps, ss->elapsed_stime, ss->temperature,
                    ss->overpotential, ss->atom_cnt, ss->total_internal_energy);
    if (ls->output_state_csv) {
        log_mc_state_csv(ls->state_csv, ls->framenum, mmc_steps, ss->iter,
                         ss->total_internal_energy);
    }
    if (ls->output_xyz) {
        snprintf(suffix + strlen(suffix), BUFFER_SIZE - strlen(suffix), "_final");
        write_xyz_file(ls->position_log_prefix, ls->framenum, suffix, ss, se);
    }

    if ((ss->final_iteration > 0) && (mmc_steps >= ss->final_iteration)) {
        fprintf(ls->sim_log, "Reached final iteration and terminated\n");
    }

    printf("Finished simulation\n");

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
