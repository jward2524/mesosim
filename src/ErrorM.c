#include "ErrorM.h"
#include "State.h"
#include <errno.h>
#include <stdlib.h>
#include <string.h>

static struct SimulationState *sim_state = NULL;
static struct SimulationEnv *sim_env = NULL;
static struct LoggingState *log_state = NULL;

static void free_if_exists(void **pointer);

// TODO: make setter functions for state variables
void set_state(struct SimulationState *ss, struct SimulationEnv *se, struct LoggingState *ls)
{
    sim_state = ss;
    sim_env = se;
    log_state = ls;
}

// frees pointer only if it isn't NULL and sets pointer to NULL after free
static void free_if_exists(void **pointer)
{
    if (*pointer == NULL)
    {
        return;
    }
        
    free(*pointer);
    *pointer = NULL;
    return;
}

// emits generic error message to log file, frees allocated memory, and exits
void clean_and_exit(int error)
{
    // errors during: reading input, m/calloc'ing, usage(), making temp file
    
    if (error != 0)
    {
        fprintf(log_state->sim_log_file, "Error encountered - check stderr\n");
        fprintf(log_state->sim_log_file, "%s", strerror(errno));
    }

    // SimulationState
    if (sim_state != NULL)
    {
        for (int i = 0; i < sim_state->atom_cnt; i++)
        {
            free_if_exists((void **)&(sim_state->atom_arr[i]));
        }
        free_if_exists((void **)&sim_state->atom_arr);
        
        for (int i = 0; i < sim_state->transition_cnt; i++)
        {
            free_if_exists((void **)&sim_state->transition_arr[i]);
        }
        free_if_exists((void **)&sim_state->transition_arr);

        free_if_exists((void **)&sim_state->rate_arr);
        free_if_exists((void **)&sim_state); // we know it exists, but still useful
    }

    // SimulationEnv
    if (sim_env != NULL)
    {
        free_if_exists((void **)&sim_env->atom_names);
        free_if_exists((void **)&sim_env->substrate_composition);
        free_if_exists((void **)&sim_env->nn_energy);
        free_if_exists((void **)&sim_env->is_soluble);
        free_if_exists((void **)&sim_env->transition_vectors);
        free_if_exists((void **)&sim_env->opposite_tvectors);
        free_if_exists((void **)&sim_env);
    }

    // LoggingState
    if (log_state != NULL)
    {
        fclose(log_state->sim_log_file);
        free_if_exists((void **)&log_state->log_list);
        free_if_exists((void **)&log_state);
    }

    if (error !=0) 
    {
        exit(error);
    }
}
