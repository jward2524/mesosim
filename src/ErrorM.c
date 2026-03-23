#include "ErrorM.h"
#include "FileIO.h"
#include "State.h"
#include <errno.h>
#include <setjmp.h>
#include <stdlib.h>
#include <string.h>

#ifdef TEST
#define INTERNAL
#include "unity.h"
volatile int exit_flag = 0;
volatile int exit_errno = 0;
volatile int jmp_set = 0;
int expected_exit_errno = 0;
jmp_buf test_exit_jmp;
#else
#define INTERNAL static
#endif

INTERNAL struct SimulationState *sim_state = NULL;
INTERNAL struct SimulationEnv *sim_env = NULL;
INTERNAL struct LoggingState *log_state = NULL;

void set_state(struct SimulationState *ss, struct SimulationEnv *se, struct LoggingState *ls)
{
    sim_state = ss;
    sim_env = se;
    log_state = ls;
}

void initialize_states(struct SimulationState **ss, struct SimulationEnv **se,
                       struct LoggingState **ls)
{
    *ss = calloc(1, sizeof(struct SimulationState));
    *se = calloc(1, sizeof(struct SimulationEnv));
    *ls = calloc(1, sizeof(struct LoggingState));
    set_state(*ss, *se, *ls);

    if (*ss == NULL) {
        perror("Couldn't allocate memory for simulation state");
        clean_and_error(errno);
    }
    if (*se == NULL) {
        perror("Couldn't allocate memory for simulation environment");
        clean_and_error(errno);
    }
    if (*ls == NULL) {
        perror("Couldn't allocate memory for logging state");
        clean_and_error(errno);
    }
}

// frees pointer only if it isn't NULL and sets pointer to NULL after free
INTERNAL void free_if_exists(void **pointer)
{
    if (*pointer == NULL) {
        return;
    }

    free(*pointer);
    *pointer = NULL;
    return;
}

#ifdef DUMP_CORE
#include <sys/types.h> // for pid_t
#include <sys/wait.h> // for waitpid
#include <unistd.h> // for getpid, fork, execlp
void create_coredump(void)
{
    pid_t parent_pid = getpid();

    // fork() returns zero in child process, child PID in parent process
    pid_t child_pid = fork();

    // to avoid issues with yama ptrace_scope, dump core of the child from the parent
    // child continues, parent ends
    if (child_pid == 0) {
        char attach_str[32];
        char gcore_str[32];
        snprintf(attach_str, sizeof(attach_str), "attach %d", parent_pid);
        snprintf(gcore_str, sizeof(gcore_str), "generate-core-file core%d.%d", log_state->xyz_framenum, parent_pid);

        // execute gcore from PATH, with argv[] parameters
        // execlp(const char *file, const char *arg, ..., NULL);
        // file will be searched for in PATH if no slash
        // first `arg` must be filename with file being executed, and last arg must be in NULL
        // (bc variadic parameters)
        int ret = execlp("gdb", "gdb", "-batch-silent", "-ex", attach_str, "-ex", gcore_str, "-ex", "detach", "-ex", "quit", NULL);
        // "No such file or directory" errors are normal for shared libraries that aren't compiled with debug symbols
        // "target file /proc/1197/cmdline contained unexpected null characters" is normal and can be ignored

        // exec replaces current process image with new process image
        // only returns if an error occurs, and will set ret = -1
        perror("execlp failed");
        call_exit(EXIT_FAILURE);
    } else {
        // parent process - wait for child to finish dumping core
        int status;
        if (waitpid(child_pid, &status, 0) == -1) {
            perror("waitpid failed");
            call_exit(EXIT_FAILURE);
        }
    }
}
#endif

void call_exit(int error_num)
{
#ifdef TEST
    exit_flag = 1;
    exit_errno = error_num;
    if (!jmp_set) {
        jmp_set = 0;
        fprintf(stderr, "Error: setjmp was not called - exit call was unexpected\nIf expected, wrap test with EXPECT_EXIT\n");
        TEST_FAIL_MESSAGE("Exit call was unexpected");
    }
    jmp_set = 0;
    longjmp(test_exit_jmp, 1);
#else
    exit(error_num);
#endif
}

// emits generic error message to log file, frees allocated memory, and exits
void clean_and_error(int exit_error)
{
    if (exit_error != 0) {
        FILE *fp = log_state->sim_log ? log_state->sim_log : stderr;
        fprintf(fp, "Error encountered - check stderr\n");
        safe_log(fp, "%d: %s\n", errno, strerror(errno));

#if !defined(NDEBUG) && !defined(TEST)
        // if in debug mode, abort to get a core dump
        fprintf(stderr, "Creating core dump\n");
        abort();
#endif
    }

    // if in debug mode, free all allocated memory and exit gracefully for easier Memcheck usage

    // SimulationState
    if (sim_state != NULL) {
        for (int i = 0; i < sim_state->atom_cnt; i++) {
            free_if_exists((void **)&(sim_state->atom_arr[i]));
        }
        free_if_exists((void **)&sim_state->atom_arr);

        for (int i = 0; i < sim_state->transition_cnt; i++) {
            free_if_exists((void **)&sim_state->transition_arr[i]);
        }
        free_if_exists((void **)&sim_state->transition_arr);

        for (int i = 0; i < sim_state->rate_cnt; i++) {
            // free_if_exists((void **)&(sim_state->rate_arr[i].atom_env));
            free(sim_state->rate_arr[i].atom_env);
        }
        free_if_exists((void **)&sim_state->rate_arr);

        free_if_exists((void **)&(sim_state->transition_probability.rate_arr_index));
        free_if_exists((void **)&(sim_state->transition_probability.lbound));
        free_if_exists((void **)&(sim_state->transition_probability.ubound));

        free_if_exists((void **)&sim_state); // we know it exists, but still useful
    }

    // SimulationEnv
    if (sim_env != NULL) {
        for (int i = 0; i < sim_env->atom_names_cnt; i++) {
            free_if_exists((void **)&(sim_env->atom_names[i]));
        }
        free_if_exists((void **)&sim_env->atom_names);
        free_if_exists((void **)&sim_env->substrate_composition);
        free_if_exists((void **)&sim_env->nn_energy);
        free_if_exists((void **)&sim_env->is_soluble);
        free_if_exists((void **)&sim_env->transition_vectors);
        free_if_exists((void **)&sim_env->opposite_tvectors);
        free_if_exists((void **)&sim_env->atoms_per_nn_level);
        free_if_exists((void **)&sim_env);
    }

    // LoggingState
    if (log_state != NULL) {
        if (log_state->sim_log)
            fclose(log_state->sim_log);
        if (log_state->steps_csv)
            fclose(log_state->steps_csv);
        if (log_state->state_csv)
            fclose(log_state->state_csv);
        if (log_state->csv_fields) {
            for (int i = 0; i < log_state->csv_field_count; i++) {
                free_if_exists((void **)&log_state->csv_fields[i]);
            }
        }
        free_if_exists((void **)&log_state->csv_fields);
        free_if_exists((void **)&log_state->csv_field_funcs);
        free_if_exists((void **)&log_state->csv_schedule.list);
        free_if_exists((void **)&log_state->xyz_schedule.list);
        free_if_exists((void **)&log_state);
    }

    // if no error, don't exit (used in tearDown in tests)
    if (exit_error != 0) {
        call_exit(exit_error);
    }
}
