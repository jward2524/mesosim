#include "KMC.h"
#include "Atoms.h"
#include "Checkpoint.h"
#include "ErrorM.h"
#include "FileIO.h"
#include "Random.h"
#include "Simulation.h"
#include "Utils.h"
#include <math.h>
#include <stdlib.h>

static double calculate_surf_diffusion_rate(unsigned char *atom_env, int final_config_neighbor_cnt,
                                            double temperature, struct SimulationEnv *se);
static double calculate_evaporation_rate(unsigned char *atom_env, double temperature,
                                         double overpotential, struct SimulationEnv *se);

unsigned long perform_kmc(struct SimulationState *ss, struct SimulationEnv *se,
                          struct LoggingState *ls)
{
    // start with everything current to the current state
    // choose a transition
    // perform the transition
    // update time
    // update rates
    // next iteration

    // prev and current iteration simulation times, for overpotential moving
    double cur_stime, prev_stime;
    double transition_type_probability;

    prev_stime = ss->elapsed_stime;

    // containers for coordinates
    int lastxt, lastyt, lastzt;
    int old_x, old_y, old_z;
    // int neighbor_x, neighbor_y, neighbor_z, neighbor_idx;

    bool simulation_run = true;
    while (simulation_run) {
        ss->iter++;

        // pick the type of transition to occur
        double rand1 = dran(&se->rand_state);
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
                double rand2 = dran(&se->rand_state);
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

                coord = get_coordination(transitioning_atom_idx, ss, se);

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

                    unsigned char atype = ss->atom_arr[transitioning_atom_idx]->type;

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
        double rand3 = ((double)ran(&se->rand_state) + 1) / ((double)ULLONG_MAX + 1);
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

        output_on_schedule(&step_data, ss, se, ls);
        checkpoint_on_schedule(ss, se, ls);

        // check if simulation is over
        if (ss->sim_end_type == SIM_END_BY_STIME) {
            simulation_run = (ss->elapsed_stime < ss->run_stime);
        } else if (ss->sim_end_type == SIM_END_BY_ITERATIONS) {
            simulation_run = (ss->iter < ss->final_iteration);
        }
    }

    return 0;
}

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

// calculates surface diffusion rate using atom environment
static double calculate_surf_diffusion_rate(unsigned char *atom_env, int final_config_neighbor_cnt,
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
static double calculate_evaporation_rate(unsigned char *atom_env, double temperature,
                                         double overpotential, struct SimulationEnv *se)
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
