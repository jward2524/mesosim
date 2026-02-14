#include "InputXYZ.h"
#include "Utils.h"
#include <stdlib.h>
#include <string.h>

// mallocs strings in KVs
/**
 * @brief parses a key-value string to identify key and value
 *
 * @param kv_str [lowercase] null-terminated string containing key-value pair separated by '='
 * @param kv_len length of kv_str
 * @param kv struct to store the key and value
 * @return int status
 */
int parse_key_value(const char *kv_str, size_t kv_len, struct KV *kv)
{
    const char *key_start = kv_str;
    const char *value_start;
    // size_t key_len, value_len;
    // int has_eq = 0;

    const char *eq = memchr(kv_str, '=', kv_len);
    if (!eq) {
        return 1;
    }
    value_start = eq + 1;

    // start to '='
    size_t key_len = (size_t)(eq - key_start);
    char *key = (char *)malloc(key_len + 1);
    if (!key) {
        return 2;
    }
    memcpy(key, key_start, key_len);
    key[key_len] = '\0';

    // '=' to end
    // remove quotes if escaped
    size_t value_len;
    char c = (char)*value_start;
    if ((c == '\"') || (c == '\'')) {
        value_start++;
        const char *q = strchr(value_start, c);
        if (!q)
            return 4;
        value_len = (size_t)(q - value_start);
    } else {
        value_len = (size_t)(kv_len - key_len - 1);
    }
    char *value = (char *)malloc(value_len + 1);
    if (!value) {
        return 3;
    }
    memcpy(value, value_start, value_len);
    value[value_len] = '\0';

    kv->key = key;
    kv->value = value;
    return 0;
}

// pointer to array of KVs
int parse_comment(const char *line, struct KV **outpairs, size_t *outpairs_cnt)
{
    // if (!line || !out_pairs || !out_count) return -1;

    unsigned char pairs_max_cnt = 64;
    // array of KVs
    struct KV *pairs = (struct KV *)malloc(sizeof(struct KV) * pairs_max_cnt);
    if (!pairs) {
        fprintf(stderr, "Out of memory allocating %zu bytes\n", sizeof(struct KV) * pairs_max_cnt);
        return 1;
    }

    unsigned char pairs_cnt = 0;

    // loop over characters in line
    int has_eq = 0;
    const char *p = line;
    while (*p) {
        if (pairs_cnt >= pairs_max_cnt) {
            fprintf(
                stderr,
                "There too many key-value pairs in comment line - max %d. Last parsed key was %s\n",
                pairs_max_cnt, pairs[pairs_cnt].key);
            // free(pairs);
            // for (int i = 0; i < pairs_cnt; i++)
            // {
            // 	if (pairs[i].key) free(pairs[i].key);
            // 	if (pairs[i].value) free(pairs[i].value);
            // }
            return 1;
        }

        // trim leading whitespace
        while (*p && isspace((char)*p))
            p++;
        if (!*p)
            break;

        const char *start = p;
        int in_quote = 0; // tracks first quote with 0, or '\'' or '"'
        const char *end = NULL;

        // scan until next unquoted whitespace -> token
        while (*p) {
            char c = p[0];
            if ((c == '\'') || (c == '"')) {
                if (!in_quote) {
                    in_quote = c;
                } else if (in_quote == c) {
                    in_quote = 0; /* closing quote */
                }
                p++;
                continue;
            }
            if (!in_quote && (c == '=')) {
                has_eq = 1;
            }
            if (!in_quote && isspace(c)) {
                end = p;
                break;
            }
            p++;
        }
        // if (!end) end = p; // token goes to end of string

        // only keep tokens that contain '='
        if (!has_eq)
            continue;

        // trim whitespace from token boundaries
        const char *t_start = start;
        while (t_start < end && isspace((unsigned char)*t_start))
            t_start++;
        const char *t_end = end;
        while (t_end > t_start && isspace((unsigned char)*(t_end - 1)))
            t_end--;

        if (t_end <= t_start) {
            // empty token after trimming
            continue;
        }

        // copy token into a new string
        size_t len = (size_t)(t_end - t_start);
        char *token = (char *)malloc(len + 1);
        if (!token) {
            return -2;
        }
        memcpy(token, t_start, len);
        token[len] = '\0';
        for (int i = 0; i < (int)len; i++) {
            token[i] = (char)tolower(token[i]);
        }

        // parse lowercase token into key-value pair
        parse_key_value(token, len, &pairs[pairs_cnt]);
        free(token);
        pairs_cnt++;
    }
    *outpairs = pairs;
    *outpairs_cnt = pairs_cnt;
    return 0;
}

// Parse Properties value like: "species:S:1:pos:R:3:force:R:3"
// pointer to array of PropertyDesc
// mallocs array of PropertyDesc
int parse_properties_value(const char *propval, PropertyDesc **out_props, int *out_nprops)
{
    size_t nprops_max = 64;
    size_t nprops = 0;
    PropertyDesc *props = (PropertyDesc *)malloc(sizeof(PropertyDesc) * nprops_max);
    if (!props) {
        fprintf(stderr, "Out of memory allocating %zu bytes\n", sizeof(PropertyDesc) * nprops_max);
        return 1;
    }

    const char *p = propval;
    while (*p) {
        if (nprops >= nprops_max) {
            fprintf(stderr, "There too many key-value pairs in comment line - max %zu\n",
                    nprops_max);
            free(props);
            return 1;
        }

        // trim starting whitespace
        while (*p && isspace((unsigned char)*p))
            ++p;
        if (!*p)
            break;

        // read name up to ':'
        const char *name_start = p;
        while (*p && (*p != ':') && !isspace((unsigned char)*p))
            ++p;
        if (*p != ':') {
            fprintf(stderr, "Error in 'properties' value: %s", name_start);
            break;
        }

        size_t name_len = (size_t)(p - name_start);

        // if name is zero length or is too long for PropertyDesc struct, error
        // .name includes a byte for null-terminating char
        if (name_len == 0 || name_len >= sizeof(props[nprops].name)) {
            fprintf(stderr, "Property name '%.*s' is too long (>%zu characters)\n",
                    (int)name_len - 1, name_start, sizeof(props[nprops].name));
            break;
        }

        // copy name into struct and terminate with null char
        memcpy(props[nprops].name, name_start, name_len);
        props[nprops].name[name_len] = '\0';
        ++p; // skip ':'

        // read type char
        if (!*p)
            break;
        char type = *p;
        props[nprops].type = type;
        ++p;
        if (*p != ':') {
            fprintf(stderr, "Error in 'properties' value: %s", name_start);
            break;
        }
        ++p;

        // read ncols integer
        // checks for all consecutive digits until end of buffer size
        char buf[32];
        int bi = 0;
        while (*p && isdigit((unsigned char)*p)) {
            if (bi < (int)sizeof(buf) - 1) {
                buf[bi] = *p;
                bi++;
            }
            ++p;
        }
        buf[bi] = '\0';
        if (bi == 0) {
            fprintf(stderr, "Number of columns not found for %s", props[nprops].name);
            break;
        }
        props[nprops].ncols = atoi(buf);

        if (*p && isdigit((unsigned char)*p)) {
            fprintf(stderr, "Number of columns for %s truncated to %s", props[nprops].name, buf);
            break;
        }

        // If there's a ':' next, the next token begins; skip one ':' and continue.
        if (*p == ':')
            p++;

        nprops++;
        // Skip any spaces before next descriptor
        // while (*p && isspace((unsigned char)*p)) ++p;
    }

    if (nprops == 0) {
        free(props);
        return -1;
    }
    *out_props = props;
    *out_nprops = (int)nprops;
    return 0;
}

// Map property tokens into Atom fields
// check before use that token count is correct
int fill_atom_from_tokens(Atom *atom, char **tokens, int ntokens, char **atom_names,
                          int atom_names_cnt, PropertyDesc *props, int nprops)
{
    // Defaults
    // memset(atom, 0, sizeof(*atom));
    // atom->bsradius = 1.0; // default unless overridden by symbol mapping below

    // Check we have enough tokens for declared properties
    int expected = 0;
    for (int k = 0; k < nprops; ++k) {
        expected += props[k].ncols;
    }
    if (ntokens < expected) {
        return expected;
    }

    int tok_idx = 0;
    for (int i = 0; i < nprops; ++i) {
        PropertyDesc *prop = &props[i];
        if (prop->type == 's') {
            // string columns
            if (prop->ncols >= 1 && tokens[tok_idx]) {
                if (strncmp(prop->name, "species", 7) == 0) {
                    int res = get_type_from_name(tokens[tok_idx], atom_names, atom_names_cnt,
                                                 &(atom->type));
                    if (res > 0)
                        return res;
                }
                tok_idx += prop->ncols;
            } else {
                tok_idx += prop->ncols; // skip unimplemented properties
            }
        } else if (prop->type == 'r') {
            // real columns
            if ((prop->ncols == 3) && (!strncmp(prop->name, "pos", 3))) {
                for (int k = 0; k < 3; ++k) {
                    atom->cartesian[k] = atof(tokens[tok_idx]);
                    tok_idx++;
                }
            } else if ((prop->ncols == 1) && !strncmp(prop->name, "energy", 6)) {
                atom->energy = atof(tokens[tok_idx]);
                tok_idx += 1;
            } else {
                tok_idx += prop->ncols; // skip unimplemented properties
            }
        } else if (prop->type == 'i') {
            // no integer-type properties implemented yet
            tok_idx += prop->ncols;
        } else {
            tok_idx += prop->ncols; // skip unimplemented properties
        }
    }
    return 0;
}

// Fallback classic XYZ parser: element + x y z
void fill_atom_from_xyz(Atom *a, char **tokens, int ntok, char **atom_names, int atom_names_cnt)
{
    if (ntok >= 4) {
        int res = get_type_from_name(tokens[0], atom_names, atom_names_cnt, &(a->type));
        if (res > 0)
            return;

        a->cartesian[0] = atof(tokens[1]);
        a->cartesian[1] = atof(tokens[2]);
        a->cartesian[2] = atof(tokens[3]);
    } else {
        // Not enough tokens; leave zeros
    }
}

void clean_xyz_structs(struct KV *kvpairs, size_t kvpairs_cnt, PropertyDesc *properties)
{
    if (kvpairs) {
        for (int j = 0; j < (int)kvpairs_cnt; j++) {
            free(kvpairs[j].key);
            free(kvpairs[j].value);
        }
    }
    if (properties) {
        free(properties);
    }
}
