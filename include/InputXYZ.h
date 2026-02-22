#ifndef XYZ_PARSER_H
#define XYZ_PARSER_H
#include "State.h"
#include <ctype.h>

// key-value pairs
struct KV {
    // malloc'd in parse_key_value
    char *key;
    char *value;
};

// Property descriptor
typedef struct {
    char name[32]; // e.g. "species", "pos", "force"
    char type;     // 'S' (string), 'R' (real), 'I' (int)
    int ncols;     // number of columns for this property
} PropertyDesc;

int parse_key_value(const char *kv_str, size_t kv_len, struct KV *kv);
int parse_comment(const char *line, struct KV **outpairs, size_t *outpairs_cnt);
int parse_properties_value(const char *propval, PropertyDesc **out_props, int *out_nprops);
int fill_atom_from_tokens(Atom *atom, char **tokens, int ntokens, char **atom_names,
                          int atom_names_cnt, PropertyDesc *props, int nprops);
void fill_atom_from_xyz(Atom *a, char **tokens, int ntok, char **atom_names, int atom_names_cnt);
void clean_xyz_structs(struct KV *kvpairs, size_t kvpairs_cnt, PropertyDesc *properties);
bool parse_xyz_file(FILE *input_file, struct SimulationState *ss, struct SimulationEnv *se,
                    struct LoggingState *ls);

#endif // XYZ_PARSER_H
