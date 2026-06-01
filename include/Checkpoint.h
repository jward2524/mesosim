#ifndef CHECKPOINT_H
#define CHECKPOINT_H

#include "State.h"
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#define CHECKPOINT_MAGIC_STRING "MESOSIM2"
#define CHECKPOINT_MAGIC_SIZE 8u
#define CHECKPOINT_FORMAT_VERSION 1u
#define CHECKPOINT_ARRAY_MAGIC 0x4D32 // 'M2' in hex
#define CHECKPOINT_ARRAY_MAGIC_SIZE 2u

typedef enum {
    CHECKPOINT_OK = 0,
    CHECKPOINT_UNIMPLEMENTED = -1,
    CHECKPOINT_ERROR = -2
} CheckpointStatus;

#pragma pack(push, 1)
typedef struct {
    unsigned long iter;
    unsigned long mmc_steps;
    unsigned long final_iteration;
    double run_stime;
    bool simulation_should_kill_itself;
    double elapsed_stime;
    int sim_end_type;
    double frequency_sum;
    double total_internal_energy;
    double temperature;
    double overpotential;
    int total_atoms_dissolved;
} CheckpointStatePayload;
#pragma pack(pop)

#pragma pack(push, 1)
typedef struct {
    unsigned flavor;
    unsigned rand_seed;

    double overpotential_ramp_rate;
    double max_overpotential;

    int system_size_x;
    int system_size_y;
    int system_size_z;

    int num_elements;
    int num_nn_levels;
    int num_bond_types;
    int num_nn_types;
    int num_transition_vectors;

    int dissolution;
    int atom_names_cnt;
    int lattice_type;
} CheckpointEnvPayload;
#pragma pack(pop)

typedef struct SEAP {
    double *substrate_composition;
    int n_substrate_composition;
    double *nn_energy;
    int n_nn_energy;
    bool *is_soluble;
    int n_is_soluble;
    char **atom_names;
    int n_atom_names;
    int *n_atom_names_str; // size of each string in atom_names
} CheckpointEnvArrPayload;

#pragma pack(push, 1)
typedef struct {
    unsigned char type;
    double energy;
    int lattice_u;
    int lattice_v;
    int lattice_w;
    double bsradius;
} CheckpointAtomPayload;
#pragma pack(pop)

typedef enum CheckpointArrayFlag {
    CAF_ATOMS = 1,
    CAF_SUBSTRATE_COMPOSITION,
    // CAF_ATOMS_PER_NN_LEVEL,
    CAF_NN_ENERGY,
    CAF_IS_SOLUBLE,
    CAF_ATOM_NAMES,
    CAF_ATOM_NAME_STR,
    CAF_OUTPUT_FORMATS,
    CAF_SENTINEL,
} CheckpointArrayFlag;

// pack/pop makes the struct have 1-byte alignment instead of default (4?)
#pragma pack(push, 1)
typedef struct {
    char magic[CHECKPOINT_MAGIC_SIZE];
    uint32_t format_version;
    uint32_t header_bytes;
    uint32_t payload_bytes;
    uint32_t checksum;
    uint8_t has_ss;
    uint8_t has_se;
    uint8_t has_atoms;
    uint8_t has_ls;
    // end of header or placeholder for future features
    uint32_t reserved0;
} CheckpointHeader;
#pragma pack(pop)

// ENHANCE: macro for functions exposed only for testing
CheckpointStatus checkpoint_save(const char *path, const struct SimulationState *ss,
                                 const struct SimulationEnv *se, const struct LoggingState *ls);
CheckpointStatus checkpoint_load(const char *path, struct SimulationState *ss,
                                 struct SimulationEnv *se, const struct LoggingState *ls);

uint32_t checkpoint_checksum32(const void *data, size_t size);
void checkpoint_header_set_magic(CheckpointHeader *header);
void checkpoint_header_init(CheckpointHeader *header);
void checkpoint_header_finalize(CheckpointHeader *header, uint32_t payload_bytes,
                                uint32_t checksum);
bool checkpoint_header_has_valid_magic(const CheckpointHeader *header);
bool checkpoint_header_has_valid_structure(const CheckpointHeader *header);
bool checkpoint_header_is_valid(const CheckpointHeader *header);
bool checkpoint_header_has_valid_checksum(const CheckpointHeader *header, const void *payload,
                                          size_t payload_size);
CheckpointStatus checkpoint_header_write(FILE *file, const CheckpointHeader *header);
CheckpointStatus checkpoint_header_read(FILE *file, CheckpointHeader *header);

void fill_state_payload(CheckpointStatePayload *payload, const struct SimulationState *ss);
void apply_state_payload_to_simstate(const CheckpointStatePayload *payload,
                                     struct SimulationState *ss);

typedef struct {
    int framenum;
    int verbose;
    unsigned long verbose_interval;
    int increment_precision;
    int stime_precision;
    int overpot_precision;
} CheckpointLoggingPayload;

void fill_logging_payload(CheckpointLoggingPayload *payload, const struct LoggingState *ls);
void apply_logging_payload_to_loggingstate(const CheckpointLoggingPayload *payload,
                                           struct LoggingState *ls);

void fill_env_payload(CheckpointEnvPayload *payload, const struct SimulationEnv *se);
void apply_env_payload_to_config(const CheckpointEnvPayload *payload,
                                 struct SimulationConfig *config);
void fill_env_array_payload(CheckpointEnvArrPayload *arr_payload, const struct SimulationEnv *se);
void apply_env_arrays(CheckpointEnvArrPayload *arr_payload, struct SimulationEnv *se);
void write_env_arrays(const CheckpointEnvArrPayload *arr_payload, uint8_t **p_payload,
                      uint32_t *p_payload_bytes);
void read_env_arrays(uint8_t *payload, size_t *total_bytes_read,
                     CheckpointEnvArrPayload *arr_payload);

CheckpointStatus append_to_payload(const void *data, size_t data_size, uint8_t **p_payload,
                                   uint32_t *p_payload_bytes);
CheckpointStatus write_array_magic(uint8_t **p_payload, uint32_t *p_payload_bytes);
CheckpointStatus rebuild_rates_and_transitions(struct SimulationState *ss,
                                               struct SimulationEnv *se);

void fill_atom_payload(CheckpointAtomPayload *payload, const Atom *atom);
void apply_atom_payload(const CheckpointAtomPayload *payload, Atom *atom);
void write_atom_array(const Atom **atom_arr, uint32_t n_atoms, uint8_t **p_payload,
                      uint32_t *p_payload_bytes);
CheckpointStatus read_atom_array(const uint8_t *payload, size_t *total_bytes_read, uint32_t *out_n,
                                 CheckpointAtomPayload **out_atom_payload_arr);

CheckpointStatus validate_array_header(const uint8_t *header, uint16_t *out_flag, uint32_t *out_n);
CheckpointStatus write_array_header_into_payload(uint16_t flag, uint32_t n, uint8_t **p_payload);
CheckpointStatus write_array(uint16_t flag, uint32_t n, const void *arr, size_t elem_size,
                             uint8_t **p_payload, uint32_t *p_payload_bytes);
CheckpointStatus read_array(const uint8_t *payload, size_t *bytes_read, uint16_t *out_flag,
                            uint32_t *out_n, void **out_arr);
CheckpointStatus verify_payload_size(FILE *file, CheckpointHeader *header);

#pragma pack(push, 1)
typedef struct {
    int mode;
    double interval;
    int list_len;
    int list_idx;
    double next_checkpoint;
    int frame_num;
} CheckpointOutputSchedulePayload;
#pragma pack(pop)

#pragma pack(push, 1)
typedef struct {
    CheckpointOutputSchedulePayload schedule;
    char filename[256];
    int field_count;
    int frame_num;
} CheckpointOutputFormatCsvPayload;
typedef struct {
    CheckpointOutputSchedulePayload schedule;
    char prefix[256];
    char suffix[256];
    int frame_num;
    uint8_t stripped;
} CheckpointOutputFormatXyzPayload;
typedef struct {
    char filename[256];
    uint8_t with_coordination;
} CheckpointOutputFormatStepsPayload;
typedef struct {
    uint8_t type;
    uint8_t is_active;
    union {
        CheckpointOutputFormatCsvPayload csv;
        CheckpointOutputFormatXyzPayload xyz;
        CheckpointOutputFormatStepsPayload steps;
    } data;
} CheckpointOutputFormatPayload;
#pragma pack(pop)

typedef struct {
    CheckpointOutputFormatPayload *formats;
    int n_out_formats;
} CheckpointOutputFormatArrPayload;

/* Output format array helpers */
CheckpointStatus fill_output_format_array_payload(CheckpointOutputFormatArrPayload *arr_payload,
                                                  const struct LoggingState *ls);
void apply_output_format_array_to_loggingstate(const CheckpointOutputFormatArrPayload *arr_payload,
                                               struct LoggingState *ls);
CheckpointStatus write_output_format_array(const CheckpointOutputFormatArrPayload *arr_payload,
                                           uint8_t **p_payload, uint32_t *p_payload_bytes);
CheckpointStatus read_output_format_array(uint8_t *payload, size_t *total_bytes_read,
                                          CheckpointOutputFormatArrPayload *arr_payload);

#endif // CHECKPOINT_H
