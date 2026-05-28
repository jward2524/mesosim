#include "Simulation.h"
#include "Atoms.h"
#include "ErrorM.h"
#include "FileIO.h"
#include "Mesosim.h"
#include "Utils.h"
#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// static const double FABS_TOL = 1e-6;

// ENHANCE: pass struct with all simulation parameters as argument
// potentially FILE* as arguments
unsigned long perform_simulation(struct SimulationState *ss, struct SimulationEnv *se,
                                 struct LoggingState *ls)
{
    // prev and current iteration simulation times, for overpotential moving
    double cur_stime, prev_stime;

    double transition_type_probability;

    // double vap;
    unsigned char atype;

    prev_stime = ss->elapsed_stime;

    bool simulation_end = false;

    for (int i = 0; i < ss->atom_cnt; ++i) {
        // resets all kinetic paramters
        refresh_transitions(i, ss, se);
    }
    compute_transition_array(ss, se);

    if (ls->output_xyz) {
        write_xyz_suffix(ls->xyz_suffix, ls->xyz_schedule.mode, 0.);
    }
    write_logs(ls->output_state_csv, ls->output_xyz, ss, se, ls);

    // start with everything current to the current state
    // choose a transition
    // perform the transition
    // update time
    // update rates
    // next iteration

    // containers for coordinates
    int lastxt, lastyt, lastzt;
    int old_x, old_y, old_z;
    // int neighbor_x, neighbor_y, neighbor_z, neighbor_idx;

    while (!simulation_end) {
        ss->iter++;

        if (ss->simulation_should_kill_itself) // abort simulation (only happens if atoms overlap)
        {
            ss->simulation_should_kill_itself = false;
            fprintf(stderr, "killed somehow\n");
            return 1; // return 1 b/c error
        }

        // pick the type of transition to occur
        double rand1 = drand();
        transition_type_probability = rand1;
        long rate_skip = ss->rate_cnt / 2;
        long j = rate_skip;

        // used to track diffusion/evaporation vs deposition
        int transition_found = 0;
        int is_evaporation = -1;

        int transition_jump_vector;
        long transitioning_atom_idx;
        long transitioned_atom_idx;
        long selected_transition_idx;
        int coord;

        // binary search to select transition
        while (!transition_found) {
            // use >=/<= because transition_type_probability can be 0 or 1
            if ((transition_type_probability >= ss->transition_probability.lbound[j]) &&
                (transition_type_probability <= ss->transition_probability.ubound[j])) {
                long selected_rate = ss->transition_probability.rate_arr_index[j];

                // pick the lucky atom
                // [ ]: third random number?
                // picks a type of transition and then which atom that has that transition will it
                // act on?
                double rand2 = drand();
                selected_transition_idx =
                    ss->rate_arr[selected_rate].transition_start_idx +
                    (long)(rand2 * (double)(ss->rate_arr[selected_rate].transition_count - 1));

                // selected_transition_idx gives the location of the transition_arr, which gives
                // the info about the specific atom
                transitioning_atom_idx = ss->transition_arr[selected_transition_idx]->atom_idx;
                transition_jump_vector = ss->transition_arr[selected_transition_idx]->offset_idx;

                // if jump_vector == se->num_transition_vectors then the atom is going to evaporate
                transition_found = 1;

                old_x = ss->atom_arr[transitioning_atom_idx]->lattice[0];
                old_y = ss->atom_arr[transitioning_atom_idx]->lattice[1];
                old_z = ss->atom_arr[transitioning_atom_idx]->lattice[2];

                coord = ls->steps_coord ? get_coordination(transitioning_atom_idx, ss, se) : -1;

                // perform transition
                if (transition_jump_vector != se->num_transition_vectors) {
                    // diffusion
                    is_evaporation = 0;

                    // coordinates atom is jumping to
                    lastxt = ss->atom_arr[transitioning_atom_idx]->lattice[0] +
                             se->transition_vectors[transition_jump_vector].dx;
                    lastyt = ss->atom_arr[transitioning_atom_idx]->lattice[1] +
                             se->transition_vectors[transition_jump_vector].dy;
                    lastzt = ss->atom_arr[transitioning_atom_idx]->lattice[2] +
                             se->transition_vectors[transition_jump_vector].dz;

                    adjust_pbc(&lastxt, &lastyt, &lastzt, se);

                    atype = ss->atom_arr[transitioning_atom_idx]->type;

                    // moves atom?
                    remove_atom(transitioning_atom_idx, ss, se);
                    transitioned_atom_idx = add_atom(lastxt, lastyt, lastzt, atype, NORMAL, ss, se);
                } else {
                    // dissolution
                    is_evaporation = 1;

                    // if not soluble, then there was an issue somewhere
                    if (!(se->is_soluble[ss->atom_arr[transitioning_atom_idx]->type])) {
                        fprintf(stderr, "Attempting to dissolve an insoluble atom - Terminating\n");
                        clean_and_error(EXIT_FAILURE);
                    }
                    ++ss->total_atoms_dissolved;
                    remove_atom(transitioning_atom_idx, ss, se); // evaporate the atom
                    transitioned_atom_idx = -1;
                }
                // end of performing transition
            } else if (transition_type_probability < (ss->transition_probability.lbound[j])) {
                // search to the left
                rate_skip = rate_skip / 2;
                if (rate_skip == 0)
                    rate_skip = 1;
                j -= rate_skip;
            } else if (transition_type_probability > (ss->transition_probability.ubound[j])) {
                // search to the right
                rate_skip = rate_skip / 2;
                if (rate_skip == 0)
                    rate_skip = 1;
                j += rate_skip;
                if (j == ss->rate_cnt) // no more options!
                    break;
            } else {
                fprintf(stderr, "Transition decision failed\n");
                clean_and_error(EXIT_FAILURE);
            }
        }

        // ? only happens iff jump_vector == se->num_transition_vectors? (dissolution?)
        if (!transition_found) {
            fprintf(stderr, "for some reason I didn't transition\n");
            clean_and_error(EXIT_FAILURE);
        }

        // increment the elapsed time (how much passed before last transition)
        // time passed, then transition occurred
        // whether this is the 1st, 2nd, or 3rd random number in the loop affects the time (diff
        // random numbers)

        // random number can't be zero, else increment=inf
        double rand3 = ((double)rand() + 1) / ((double)RAND_MAX + 1);
        double stime_increment = -log(rand3) / ss->frequency_sum;
        ss->elapsed_stime += stime_increment;
        ls->stime_precision =
            get_precision(ss->elapsed_stime, stime_increment, ls->increment_precision);
        if (ls->stime_precision < 0) {
            fprintf(
                stderr,
                "Simulation time precision calculation encountered an error: tot %le, inc %le\n",
                ss->elapsed_stime, stime_increment);
            clean_and_error(EXIT_FAILURE);
        }
        if (isinf((float)ss->elapsed_stime)) {
            fprintf(stderr, "Simulation time went infinite\n");
            clean_and_error(EXIT_FAILURE);
        }

        if (se->overpotential_ramp_rate > 0.0) {
            cur_stime = ss->elapsed_stime;
            if (ss->overpotential < se->max_overpotential) {
                double overpot_increment = (cur_stime - prev_stime) * se->overpotential_ramp_rate;
                ss->overpotential += overpot_increment;
                if (ss->overpotential > se->max_overpotential) {
                    // TODO : standardize double - equality checks
                    overpot_increment =
                        round((se->max_overpotential - ss->overpotential) * 1e10) / 1e10;
                    ss->overpotential = se->max_overpotential;
                }
                ls->overpot_precision =
                    get_precision(ss->overpotential, overpot_increment, ls->increment_precision);
                if (ls->overpot_precision < 0) {
                    fprintf(stderr,
                            "Overpotential precision calculation encountered an error: tot "
                            "%le, inc %le\n",
                            ss->overpotential, overpot_increment);
                    clean_and_error(EXIT_FAILURE);
                }
            } else if (ss->overpotential > se->max_overpotential) {
                fprintf(stderr, "Overpotential exceeded maximum\n");
            }
            prev_stime = ss->elapsed_stime;
        }

        // update rates
        update_outdated_transitions(old_x, old_y, old_z, transitioned_atom_idx, ss, se);

        // ENHANCE - not all are necessary if overpotential/temperature don't change?
        compute_transition_array(ss, se);
        // end of updating rates

        // after iteration, log if necessary
        if (ls->verbose && (ss->iter % ls->verbose_interval == 0)) {
            printf("Iteration %ld, time %le\n", ss->iter, ss->elapsed_stime);
        }

        StepData step_data = {
            .flavor = se->flavor,
            .iter = ss->iter,
            .sys_energy = ss->total_internal_energy,
            .uvw1 = {old_x, old_y, old_z},
            .uvw2 = {lastxt, lastyt, lastzt},
            .coord = coord,
            .kmc =
                {
                    .sim_time = ss->elapsed_stime,
                    .is_evap = is_evaporation,
                },
        };

        if (ls->output_steps_csv) {
            int uvw1[] = {old_x, old_y, old_z};
            // lastxyzt will have wrong values if current step was evaporation
            int uvw2[] = {lastxt, lastyt, lastzt};
            log_kmc_steps(ls->steps_csv, ss->iter, ss->elapsed_stime, ls->stime_precision,
                          ss->total_internal_energy, uvw1, uvw2, is_evaporation, coord);
        }

        output_on_schedule(&step_data, ss, se, ls);

        // check if simulation is over
        if (ss->sim_end_type == SIM_END_BY_STIME) {
            simulation_end = (ss->elapsed_stime >= ss->run_stime);
        } else if (ss->sim_end_type == SIM_END_BY_ITERATIONS) {
            simulation_end = (ss->iter >= ss->final_iteration);
        }
    }

    // write elapsed_stime to mark finish
    if ((ss->final_iteration > 0) && (ss->iter >= ss->final_iteration)) {
        if (ls->output_xyz) {
            snprintf(ls->xyz_suffix, BUFFER_SIZE, "i%lu", (unsigned long)ss->iter);
        }
        write_logs(ls->output_state_csv, ls->output_xyz, ss, se, ls);
        safe_log(ls->sim_log, "Reached final iteration and terminated\n");
    }
    if ((ss->run_stime > 0) && (ss->elapsed_stime >= ss->run_stime)) {
        if (ls->output_xyz) {
            snprintf(ls->xyz_suffix, BUFFER_SIZE, "t%lf", ss->elapsed_stime);
        }
        write_logs(ls->output_state_csv, ls->output_xyz, ss, se, ls);
        safe_log(ls->sim_log, "Reached end of simulation time and terminated\n");
    }

    printf("Finished simulation\n");

    return 0;
}

/******************************************************************************/
/******************************************************************************/
// updates ss->transition_probability (weighted rate list, used to choose event)
void compute_transition_array(struct SimulationState *ss, struct SimulationEnv *se)
{

    int nonzero_rate_cnt = 0;
    ss->frequency_sum = 0.0;

    double rate_const;
    // cumulative frequency sum, for computing transition probabilities and float imprecision from
    // adding small probabilities
    double *cum_frequency_sum = (double *)malloc(((size_t)ss->rate_cnt + 1) * sizeof(double));
    cum_frequency_sum[0] = 0.0;
    Rate *r;
    // TODO: only calculate if overpoential changed or if k=-1
    for (int rate_idx = 0; rate_idx < ss->rate_cnt; ++rate_idx) {
        r = &ss->rate_arr[rate_idx];
        if (r->transition_count != 0) {
            if (r->is_evaporation) {
                rate_const =
                    calculate_evaporation_rate(r->atom_env, ss->temperature, ss->overpotential, se);
            } else {
                rate_const = calculate_surf_diffusion_rate(
                    r->atom_env, r->final_config_neighbor_cnt, ss->temperature, se);
            }
            if (isnan(rate_const)) {
                fprintf(stderr, "Rate constant is nan\n");
                clean_and_error(EXIT_FAILURE);
            }
            r->k = rate_const;
            r->frequency = r->k * (double)r->transition_count;
            ss->frequency_sum += r->frequency;
            cum_frequency_sum[nonzero_rate_cnt + 1] = ss->frequency_sum;
            // sum_of_rate_populations += r->transition_count;

            ss->transition_probability.rate_arr_index[nonzero_rate_cnt] = rate_idx;
            nonzero_rate_cnt++;
        }
    }

    // now compute bounds for jump probabilities
    for (int i = 0; i < nonzero_rate_cnt; ++i) {
        double lbound = cum_frequency_sum[i] / ss->frequency_sum;
        ss->transition_probability.lbound[i] = lbound;

        double ubound = cum_frequency_sum[i + 1] / ss->frequency_sum;
        ss->transition_probability.ubound[i] = ubound;

        if (isnan(lbound) || isnan(ubound)) {
            fprintf(stderr, "Probability is nan\n");
            clean_and_error(EXIT_FAILURE);
        }
    }

    free(cum_frequency_sum);

    return;
}

void update_outdated_transitions(int old_x, int old_y, int old_z, long transitioned_atom_idx,
                                 struct SimulationState *ss, struct SimulationEnv *se)
{
    // refresh neighbors of previous site
    for (int i = 0; i < se->num_energy_contributors; i++) {
        // transitioned_atom_index - new atom index (-1 if evaporated)
        long neighbor_idx = atom_at_offset(old_x, old_y, old_z, i, ss->atom_arr, ss->zone_arr, se);

        if (neighbor_idx >= 0)
            refresh_transitions(neighbor_idx, ss, se);
    }

    // if transition was not evaporation
    if (transitioned_atom_idx >= 0) {
        // refresh moved atom - also gets refreshed as neighbor of previous site?
        refresh_transitions(transitioned_atom_idx, ss, se);

        // refresh neighbors of moved-to site
        for (int i = 0; i < se->num_energy_contributors; i++) {
            int x = ss->atom_arr[transitioned_atom_idx]->lattice[0];
            int y = ss->atom_arr[transitioned_atom_idx]->lattice[1];
            int z = ss->atom_arr[transitioned_atom_idx]->lattice[2];
            long neighbor_idx = atom_at_offset(x, y, z, i, ss->atom_arr, ss->zone_arr, se);

            if (neighbor_idx >= 0) {
                refresh_transitions(neighbor_idx, ss, se);
            }
        }
    }
}

// updates [atom_arr[atom_idx], neighbor_atom_idxs, transition_arr[i],
// rate_arr[i].transition_start_idx], initializes rate_arr
int refresh_transitions(long atom_idx, struct SimulationState *ss,
                        struct SimulationEnv *se) // atom_idx = index on atom list
{
    long neighbor_idx;
    long rate_idx; // position of rate rate in rate list rate_arr[]

    int atom_rates_cnt; // this is returned as the number of transitions (se->jump_offsets) this
                        // atom can undergo, excluding evaporation
    int start_config[MAXIMUM_NUMBER_OF_NEIGHBORS]; // -1 if empty, type if filled
    int end_config[MAXIMUM_NUMBER_OF_NEIGHBORS];

    unsigned char *atom_env =
        (unsigned char *)calloc((size_t)se->num_nn_types, sizeof(unsigned char));
    // unsigned char env_hash[se->num_nn_levels * se->num_bond_types];

    // update neighbors
    // cycle through neighbor coordinates, and check if there is an atom there
    for (int i = 0; i < se->num_transition_vectors; ++i) {
        ss->atom_arr[atom_idx]->neighbor_atom_idxs[i] = -1;

        neighbor_idx =
            atom_at_offset(ss->atom_arr[atom_idx]->lattice[0], ss->atom_arr[atom_idx]->lattice[1],
                           ss->atom_arr[atom_idx]->lattice[2], i, ss->atom_arr, ss->zone_arr, se);

        if (neighbor_idx >= 0) {
            ss->atom_arr[atom_idx]->neighbor_atom_idxs[i] = neighbor_idx;
        }

        // // if no atom atom_idx that position but occ_neighbor array says there is, fix it
        // if ((j == -1) && (ss->atom_arr[atom_idx]->neighbor_atom_idxs[i] >= 0))
        // 	ss->atom_arr[atom_idx]->neighbor_atom_idxs[i] = -1;
    }

    // update atom_env and energy
    ss->total_internal_energy -= ss->atom_arr[atom_idx]->energy;
    ss->atom_arr[atom_idx]->energy = 0;

    for (int i = 0; i < se->num_energy_contributors; ++i) {
        neighbor_idx =
            atom_at_offset(ss->atom_arr[atom_idx]->lattice[0], ss->atom_arr[atom_idx]->lattice[1],
                           ss->atom_arr[atom_idx]->lattice[2], i, ss->atom_arr, ss->zone_arr, se);

        // update atom_env
        if (neighbor_idx >= 0) {
            int atom_type = ss->atom_arr[atom_idx]->type;
            int neighbor_type = ss->atom_arr[neighbor_idx]->type;
            int env_idx = get_env_index(i, atom_type, neighbor_type, se);

            atom_env[env_idx]++;
            ss->atom_arr[atom_idx]->energy += se->nn_energy[env_idx];
        }
    }

    ss->atom_arr[atom_idx]->energy /= 2.;
    ss->total_internal_energy += ss->atom_arr[atom_idx]->energy;

    // update transitions

    // remove all mention of this atom from transition list
    // extra 1 for evaporation
    // true = 1; false = 0
    int atom_soluble = se->is_soluble[ss->atom_arr[atom_idx]->type];
    int indices = se->num_transition_vectors + atom_soluble;
    for (int i = 0; i < indices; ++i) {
        // transition can happen in the "i" direction
        if (ss->atom_arr[atom_idx]->transition_indices[i] != -1) {
            remove_from_transition_list(atom_idx, i, ss);
        }
    }

    // cycle through the neighbor sites
    // if there's an empty one, calculate the transition rate to it
    atom_rates_cnt = 0;

    // ENHANCE: redundant - already calculated initial config in j loop
    // k is number of near neighbors
    int intial_config_neighbor_cnt =
        get_initial_configuration(atom_idx, se->num_transition_vectors, ss->atom_arr, start_config);

    // skip calculating a rate of a fully coordinated atom
    if (intial_config_neighbor_cnt != se->num_transition_vectors) {
        // create transitions to each unoccupied neighbor
        int final_config_neighbor_cnt;
        for (unsigned char i = 0; i < se->num_transition_vectors; ++i) {
            // if unoccupied, consider transition
            // [?]: real use of neighbor_atoms_idxs
            if (ss->atom_arr[atom_idx]->neighbor_atom_idxs[i] == -1) {
                // end_config is only used to identify the transitions that are functionally
                // evaporations (no neighbors in end configuration)
                final_config_neighbor_cnt =
                    get_final_configuration(atom_idx, i, ss, se, end_config);

                if (final_config_neighbor_cnt <= 1) {
                    // transition in this direction leaves stray atom
                    // or leads to pair climbs
                    // these types are handled by evaporation
                    continue;
                }

                ++atom_rates_cnt;

                rate_idx = get_rate(atom_env, final_config_neighbor_cnt, 0, ss, se);
                if (rate_idx == -1) {
                    // if rate doesn't already exist, make new one
                    rate_idx = create_new_rate(atom_env, final_config_neighbor_cnt, 0, ss, se);
                }
                add_to_transition_list(rate_idx, atom_idx, i, ss, se);
            }
        }

        // if atom type is soluble, add dissolution transition
        if (atom_soluble) {
            // dissolution / evaporation transition
            final_config_neighbor_cnt = -1;
            rate_idx = get_rate(atom_env, final_config_neighbor_cnt, 1, ss, se);
            if (rate_idx == -1) {
                // if rate doesn't already exist, make new one
                rate_idx = create_new_rate(atom_env, final_config_neighbor_cnt, 1, ss, se);
            }
            // evaporation is considered to be last in se->transition_vectors
            // (not really in array but uses that index number)
            add_to_transition_list(rate_idx, atom_idx, (unsigned char)se->num_transition_vectors,
                                   ss, se);
        }
    }

    free(atom_env);
    // atom_env = NULL;

    return atom_rates_cnt; // gives number of current transitions for that atom
}

/******************************************************************************/
/******************************************************************************/

// is there a better way of doing this?
//  checks if the rate constant is already in the rate list
int get_rate(unsigned char *atom_env, int final_config_neighbor_cnt, unsigned char is_evaporation,
             struct SimulationState *ss, struct SimulationEnv *se)
{
    int env_cmp, evap_cmp, finalcnt_cmp;

    for (int i = 0; i < ss->rate_cnt; ++i) {
        env_cmp = memcmp(ss->rate_arr[i].atom_env, atom_env, (size_t)se->num_nn_types);
        env_cmp = env_cmp == 0;
        evap_cmp = ss->rate_arr[i].is_evaporation == is_evaporation;
        finalcnt_cmp = ss->rate_arr[i].final_config_neighbor_cnt == final_config_neighbor_cnt;
        if (env_cmp && evap_cmp && finalcnt_cmp)
            return i;
    }

    return -1;
}

/******************************************************************************/
/******************************************************************************/
// create new Rate struct in rate array
// updates rate_array[rate_cnt], rate_cnt
long create_new_rate(unsigned char *atom_env, int final_config_neighbor_cnt,
                     unsigned char is_evaporation, struct SimulationState *ss,
                     struct SimulationEnv *se)
{
    // ss->rate_arr[ss->rate_cnt].k = rate;
    Rate *r = &ss->rate_arr[ss->rate_cnt];
    r->transition_start_idx = ss->transition_cnt;
    r->transition_count = 0;
    r->atom_env = (unsigned char *)malloc((size_t)se->num_nn_types * sizeof(unsigned char));
    memcpy(r->atom_env, atom_env, (size_t)se->num_nn_types * sizeof(unsigned char));
    r->is_evaporation = is_evaporation;
    r->final_config_neighbor_cnt = final_config_neighbor_cnt;

    ++ss->rate_cnt;
    if (ss->rate_cnt > se->max_rates) {
        fprintf(stderr, "Number of rates (%ld) is exceeding set maximum (%ld)", ss->rate_cnt,
                se->max_rates);
        clean_and_error(errno);
    }

    return (ss->rate_cnt - 1);
}

/**
 * @brief Set the transition indices index of the atom corresponding to the provided transition
 *
 * @param t
 * @param transition_idx
 * @param atom_arr
 */
static void set_atom_transition_indices(Transition *t, long transition_idx, Atom **atom_arr)
{
    atom_arr[t->atom_idx]->transition_indices[t->offset_idx] = transition_idx;
}

// add to rate_arr[rate_idx] the atom atom_idx going in direction offset_idx
// updates transition_arr, rate_arr[rate_idx].transition_count,
// atom_arr[atom_idx]->transition_indices[offset_idx]
void add_to_transition_list(long rate_idx, long atom_idx, unsigned char offset_idx,
                            struct SimulationState *ss, struct SimulationEnv *se)
{
    // initial and final transition_arr index
    long initial_trans_idx, final_trans_idx;
    Transition **t_arr = ss->transition_arr;

    // make room for the new arrival
    // adds entry to the end of the list
    Transition *new_transition = (Transition *)malloc(sizeof(Transition));
    if (new_transition == NULL) {
        fprintf(stderr, "Couldn't allocate memory for transition %ld: %s\n", ss->transition_cnt,
                strerror(errno));
        clean_and_error(errno);
    }
    new_transition->atom_idx = atom_idx;
    new_transition->offset_idx = offset_idx;

    // make space for new transition in its rate group
    // move transitions further down the array by taking the first transition in a rate group and
    // moving it after its last transition, which will be open as we work backwards through the
    // array
    Rate *r_arr = ss->rate_arr;
    long r_cnt = ss->rate_cnt;
    for (long i = r_cnt - 1; i > rate_idx; --i) {
        initial_trans_idx = r_arr[i].transition_start_idx;
        final_trans_idx = initial_trans_idx + r_arr[i].transition_count;

        t_arr[final_trans_idx] = t_arr[initial_trans_idx];
        t_arr[initial_trans_idx] = NULL;

        // avoid like two lines
        if (initial_trans_idx != final_trans_idx) {
            set_atom_transition_indices(t_arr[final_trans_idx], final_trans_idx, ss->atom_arr);
        }
        ++r_arr[i].transition_start_idx;
    }

    // add new arrival
    long new_trans_idx = r_arr[rate_idx].transition_start_idx + r_arr[rate_idx].transition_count;

    ss->transition_arr[new_trans_idx] = new_transition;
    // manual set_atom_transition_indices to avoid two pointer derefernces
    ss->atom_arr[atom_idx]->transition_indices[offset_idx] = new_trans_idx;
    ++r_arr[rate_idx].transition_count;

    ++ss->transition_cnt;
    if (ss->transition_cnt > se->max_transitions) {
        fprintf(stderr, "Number of transitions (%ld) is exceeding set maximum (%ld)",
                ss->transition_cnt, se->max_transitions);
        clean_and_error(errno);
    }

    return;
}

// updates atom_arr[atom_idx], transition_arr[i], rate_arr[i].transition_start_idx
// removes atom jumping in the se->transition_vectors[offset_idx] direction
void remove_from_transition_list(long atom_idx, int offset_idx, struct SimulationState *ss)
{
    // old position on transition list, to be removed
    long removed_t_idx = ss->atom_arr[atom_idx]->transition_indices[offset_idx];

    Rate *r_arr = ss->rate_arr;
    // find out what Rate in rate_arr this is to update it
    // rate_idx points to the current rate list it's on
    long rate_idx = 0;
    // copy of rate_cnt to avoid multiple dereferences in loop
    long r_cnt = ss->rate_cnt;
    for (long i = 0; i < r_cnt; ++i) {
        if (removed_t_idx < (r_arr[i].transition_start_idx + r_arr[i].transition_count)) {
            rate_idx = i;
            break;
        }
    }

    // remind atom it can no longer jump
    // manual set_atom_transition_indices
    ss->atom_arr[atom_idx]->transition_indices[offset_idx] = -1;

    Transition **t_arr = ss->transition_arr;

    // free the transition
    free(t_arr[removed_t_idx]);
    t_arr[removed_t_idx] = NULL;

    // decrement the number of transitions in its rate list
    // if there are still transitions in that rate, move the last one to the removed spot
    --r_arr[rate_idx].transition_count;
    if (r_arr[rate_idx].transition_count != 0) {
        long transition_end_idx =
            r_arr[rate_idx].transition_start_idx + r_arr[rate_idx].transition_count;

        // if they are the same, no movement is necessary (and doing movement will break things)
        if (removed_t_idx != transition_end_idx) {
            t_arr[removed_t_idx] = t_arr[transition_end_idx];
            t_arr[transition_end_idx] = NULL;
            set_atom_transition_indices(t_arr[removed_t_idx], removed_t_idx, ss->atom_arr);
        }
    }

    // move transitions of rates higher on transition list up one spot
    for (long i = rate_idx + 1; i < r_cnt; ++i) {
        // move rate_arr offsets of larger indicies up one
        --r_arr[i].transition_start_idx;

        long transition_start_idx = r_arr[i].transition_start_idx;
        long transition_end_idx = r_arr[i].transition_start_idx + r_arr[i].transition_count;
        // count is always at least 1

        t_arr[transition_start_idx] = t_arr[transition_end_idx];
        t_arr[transition_end_idx] = NULL;

        // update the transition index in the corresponding atom in atom_arr to have the new
        // (lower) transition index
        set_atom_transition_indices(t_arr[transition_start_idx], transition_start_idx,
                                    ss->atom_arr);
    }

    // if the rate's count is zero, remove it and move rates up one spot
    if (r_arr[rate_idx].transition_count == 0) {
        // free atom_env, everything else will be overwritten
        free(r_arr[rate_idx].atom_env);

        for (long i = rate_idx + 1; i < r_cnt; ++i) {
            // TODO: use pointers for rate array? to make this more logical
            memcpy(&r_arr[i - 1], &r_arr[i], sizeof(Rate));
        }
        // zero last rate (don't free atom_env bc still being used by rate_cnt-2)
        memset(&r_arr[r_cnt - 1], 0, sizeof(Rate));

        --ss->rate_cnt;
        if (ss->rate_cnt < 0) {
            fprintf(stderr, "Number of rates (%ld) has dropped below zero", ss->rate_cnt);
            clean_and_error(errno);
        }
    }

    --ss->transition_cnt;
    if (ss->transition_cnt < 0) {
        fprintf(stderr, "Number of transitions (%ld) has dropped below zero", ss->transition_cnt);
        clean_and_error(errno);
    }

    return;
}

// calculates surface diffusion rate using atom environment
double calculate_surf_diffusion_rate(unsigned char *atom_env, int final_config_neighbor_cnt,
                                     double temperature, struct SimulationEnv *se)
{
    // ENHANCE: pass the number of nearest neighbors? (best practice)
    double energy = 0.0;

    // why is this static, not const? it should also be a simulation input
    // "optional" anistropy factor, indices match se->transition_vectors
    // static double b_anisotropy_factor[12] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};

    int neighbor_cnt_initial = 0;
    for (int i = 0; i < se->num_nn_types; i++) {
        energy += se->nn_energy[i] * atom_env[i];
        // TODO: anisotropy factor would require storing directions in Rate (probably as sum, to
        // be multiplied with energy)
        neighbor_cnt_initial++;
    }

    // these override the previous energy sum
    if (neighbor_cnt_initial == 0) {
        // no neighbors - this condition corresponds to a diffuser walking through a lattice
        // aka a lattice gas
        energy = -1.0;
    }

    // to simulate Ehrlich-Schwoebel barrier
    // final_config_neighbor_cnt <= 1 don't transition, handled by refresh_transitions()
    int barrier_factor = 4;
    if (final_config_neighbor_cnt <= 2) {
        energy = barrier_factor * energy;
    }

    // ENHANCE: replace calculating the exp with memoizing up the value (uhash?) -> speedup?
    // overpotential is only for evaporation
    //*rate = 1e13*exp(-energy/(kBoltz*temperature)) //+
    // 1e-4*exp(-(energy-overpotential)/(kBoltz*temperature));
    return 1e13 * exp(-energy / (kBoltz * temperature));
}

// calculates evaporation rate using atom environment
double calculate_evaporation_rate(unsigned char *atom_env, double temperature, double overpotential,
                                  struct SimulationEnv *se)
{
    double energy = 0.0;

    // static double b_anisotropy_factor[12] = {1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};
    // // optional anistropy factor

    for (int i = 0; i < se->num_nn_types; i++) {
        energy += se->nn_energy[i] * atom_env[i];
        // TODO: anisotropy factor would require storing directions in Rate (probably as sum, to
        // be multiplied with energy)
    }

    // dissolution/evaporation equation
    return 1e4 * exp(-(energy - overpotential) / (kBoltz * temperature));
}
