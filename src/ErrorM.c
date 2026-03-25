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

INTERNAL struct SimulationState **gp_sim_state = NULL;
INTERNAL struct SimulationEnv **gp_sim_env = NULL;
INTERNAL struct LoggingState **gp_log_state = NULL;

void set_state(struct SimulationState **p_ss, struct SimulationEnv **p_se,
               struct LoggingState **p_ls)
{
    gp_sim_state = p_ss;
    gp_sim_env = p_se;
    gp_log_state = p_ls;
}

void initialize_states(struct SimulationState **p_ss, struct SimulationEnv **p_se,
                       struct LoggingState **p_ls)
{
    *p_ss = calloc(1, sizeof(struct SimulationState));
    *p_se = calloc(1, sizeof(struct SimulationEnv));
    *p_ls = calloc(1, sizeof(struct LoggingState));
    set_state(p_ss, p_se, p_ls);

    if (*p_ss == NULL) {
        perror("Couldn't allocate memory for simulation state");
        clean_and_error(errno);
    }
    if (*p_se == NULL) {
        perror("Couldn't allocate memory for simulation environment");
        clean_and_error(errno);
    }
    if (*p_ls == NULL) {
        perror("Couldn't allocate memory for logging state");
        clean_and_error(errno);
    }
}

// frees pointer only if it isn't NULL and sets pointer to NULL after free
void free_if_exists(void **pointer)
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
#include <sys/wait.h>  // for waitpid
#include <unistd.h>    // for getpid, fork, execlp
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
        snprintf(gcore_str, sizeof(gcore_str), "generate-core-file core%d.%d",
                 (*gp_log_state)->xyz_framenum, parent_pid);

        // execute gcore from PATH, with argv[] parameters
        // execlp(const char *file, const char *arg, ..., NULL);
        // file will be searched for in PATH if no slash
        // first `arg` must be filename with file being executed, and last arg must be in NULL
        // (bc variadic parameters)
        int ret = execlp("gdb", "gdb", "-batch-silent", "-ex", attach_str, "-ex", gcore_str, "-ex",
                         "detach", "-ex", "quit", NULL);
        (void)ret;
        // "No such file or directory" errors are normal for shared libraries that aren't compiled
        // with debug symbols "target file /proc/1197/cmdline contained unexpected null characters"
        // is normal and can be ignored

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
        fprintf(stderr, "Error: setjmp was not called - exit call was unexpected\nIf expected, "
                        "wrap test with EXPECT_EXIT\n");
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
    // convenience pointers
    struct SimulationState *ss = gp_sim_state ? *gp_sim_state : NULL;
    struct SimulationEnv *se = gp_sim_env ? *gp_sim_env : NULL;
    struct LoggingState *ls = gp_log_state ? *gp_log_state : NULL;

    if (exit_error != 0) {
        FILE *fp = (ls && ls->sim_log) ? ls->sim_log : stderr;
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
    if (ss != NULL) {
        for (int i = 0; i < ss->atom_cnt; i++) {
            free_if_exists((void **)&(ss->atom_arr[i]));
        }
        free_if_exists((void **)&ss->atom_arr);

        for (int i = 0; i < ss->transition_cnt; i++) {
            free_if_exists((void **)&ss->transition_arr[i]);
        }
        free_if_exists((void **)&ss->transition_arr);

        for (int i = 0; i < ss->rate_cnt; i++) {
            // free_if_exists((void **)&(ss->rate_arr[i].atom_env));
            free(ss->rate_arr[i].atom_env);
        }
        free_if_exists((void **)&ss->rate_arr);

        free_if_exists((void **)&(ss->transition_probability.rate_arr_index));
        free_if_exists((void **)&(ss->transition_probability.lbound));
        free_if_exists((void **)&(ss->transition_probability.ubound));

        free_if_exists((void **)gp_sim_state);
        gp_sim_state = NULL;
    }

    // SimulationEnv
    if (se != NULL) {
        for (int i = 0; i < se->atom_names_cnt; i++) {
            free_if_exists((void **)&(se->atom_names[i]));
        }
        free_if_exists((void **)&se->atom_names);
        free_if_exists((void **)&se->substrate_composition);
        free_if_exists((void **)&se->nn_energy);
        free_if_exists((void **)&se->is_soluble);
        free_if_exists((void **)&se->transition_vectors);
        free_if_exists((void **)&se->opposite_tvectors);
        free_if_exists((void **)&se->atoms_per_nn_level);

        free_if_exists((void **)gp_sim_env);
        gp_sim_env = NULL;
    }

    // LoggingState
    if (ls != NULL) {
        if (ls->sim_log)
            fclose(ls->sim_log);
        if (ls->steps_csv)
            fclose(ls->steps_csv);
        if (ls->state_csv)
            fclose(ls->state_csv);
        if (ls->csv_fields) {
            for (int i = 0; i < ls->csv_field_count; i++) {
                free_if_exists((void **)&ls->csv_fields[i]);
            }
        }
        free_if_exists((void **)&ls->csv_fields);
        free_if_exists((void **)&ls->csv_field_funcs);
        free_if_exists((void **)&ls->csv_schedule.list);
        free_if_exists((void **)&ls->xyz_schedule.list);

        free_if_exists((void **)gp_log_state);
        gp_log_state = NULL;
    }

    // if no error, don't exit (used in tearDown in tests)
    if (exit_error != 0) {
        call_exit(exit_error);
    }
}
