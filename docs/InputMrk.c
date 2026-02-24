#include "Input.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAX_PARAMS 64
#define MAX_TOKEN  128

typedef struct {
    char name[MAX_TOKEN];
    int optional;
} ParsedParam;

const char *cmd_category_names[] = {
    "Uncategorized",
    "Geometry",
    "Thermodynamics",
    "Output and Logging",
    "Run Control",
};

/**
 * @brief Parse parameters from a command usage string, identifying parameter names and whether they are optional (enclosed in brackets)
 * 
 * @param usage 
 * @param params 
 * @return int 
 */
static int parse_usage_params(const char *usage, ParsedParam *params)
{
    if (!usage) return 0;

    int count = 0;
    int optional_mode = 0;

    char buffer[1024];
    strncpy(buffer, usage, sizeof(buffer) - 1);
    buffer[sizeof(buffer) - 1] = '\0';

    char *token = strtok(buffer, " ");

    // TODO: replace with tokenize_line
    /* First token is command name, skip it */
    token = strtok(NULL, " ");

    while (token && count < MAX_PARAMS) {

        /* Detect optional blocks */
        if (strchr(token, '[')) {
            optional_mode = 1;
        }

        /* Clean token from brackets */
        char cleaned[MAX_TOKEN];
        int j = 0;
        for (int i = 0; token[i] && j < MAX_TOKEN - 1; ++i) {
            if (token[i] != '[' && token[i] != ']')
                cleaned[j++] = token[i];
        }
        cleaned[j] = '\0';

        /* Skip empty */
        if (strlen(cleaned) > 0) {
            strncpy(params[count].name, cleaned, MAX_TOKEN);
            params[count].optional = optional_mode;
            count++;
        }

        if (strchr(token, ']')) {
            optional_mode = 0;
        }

        token = strtok(NULL, " ");
    }

    return count;
}

/* --------------------------------------------------------- */
/* Print example usage line (without brackets)               */
/* --------------------------------------------------------- */

static void print_example_line(FILE *fp, const char *usage)
{
    for (const char *p = usage; *p; ++p) {
        if (*p != '[' && *p != ']')
            fputc(*p, fp);
    }
    fputc('\n', fp);
}

/**
 * @brief Generate markdown documentation for input file format based on command definitions
 *
 * @param fp
 * @param cmds
 */
void generate_input_markdown(FILE *fp, const Command *cmds)
{
    fprintf(fp, "# Simulation Input File Format\n\n");

    /* ============ 1. Example Input File ============= */

    fprintf(fp, "## 1. Example Input File\n\n");
    fprintf(fp, "```text\n");

    /* Print one example per required command */
    for (int i = 0; cmds[i].name != NULL; ++i) {
        // arbitrary choice of kmc
        if ((cmds[i].required_kmc == CMDREQ_REQUIRED) && cmds[i].usage)
            print_example_line(fp, cmds[i].usage);
    }

    fprintf(fp, "```\n\n");

    /* =============== 2. Conceptual Grouping =============== */
    fprintf(fp, "## 2. Conceptual Organization\n\n");

    // Collect unique categories
    for (int i = 0; cmds[i].name != NULL; ++i) {

        const enum CmdCategory cat = cmds[i].category;
        int printed = 0;

        /* Check if category already printed */
        for (int j = 0; j < i; ++j) {
            if ((cmds[j].category >= 0) && (cmds[j].category == cat)) {
                printed = 1;
                break;
            }
        }

        if (!printed && cat) {

            fprintf(fp, "### %s\n\n", cmd_category_names[cat]);
            fprintf(fp, "| Command | Required | Description |\n");
            fprintf(fp, "|----------|----------|-------------|\n");

            for (int k = 0; cmds[k].name != NULL; ++k) {
                if ((cmds[k].category >= 0) && (cmds[k].category == cat)) {

                    fprintf(fp, "| `%s` | %s | %s |\n",
                            cmds[k].name,
                            cmds[k].required_kmc ? "Yes" : "No",
                            cmds[k].description ? cmds[k].description : "");
                }
            }

            fprintf(fp, "\n");
        }
    }

    /* ============== 3. Detailed Command Reference ================= */
    fprintf(fp, "## 3. Detailed Command Reference\n\n");

    for (int i = 0; cmds[i].name != NULL; ++i) {

        fprintf(fp, "### `%s`\n\n", cmds[i].name);

        fprintf(fp, "**Category:** %s  \n", cmd_category_names[cmds[i].category + 1]);

        fprintf(fp, "**Required:** %s\n\n",
                (cmds[i].required_kmc == CMDREQ_REQUIRED) ? "Yes" : "No");

        if (cmds[i].usage) {
            fprintf(fp, "**Syntax**\n\n");
            fprintf(fp, "```text\n%s\n```\n\n", cmds[i].usage);
        }

        /* Parse parameters */
        ParsedParam params[MAX_PARAMS];
        int pcount = parse_usage_params(cmds[i].usage, params);

        if (pcount > 0) {
            fprintf(fp, "**Parameters**\n\n");
            for (int p = 0; p < pcount; ++p) {
                fprintf(fp, "- `%s` — %s\n",
                        params[p].name,
                        params[p].optional ? "optional" : "required");
            }
            fprintf(fp, "\n");
        }

        if (cmds[i].description) {
            fprintf(fp, "**Description**\n\n%s\n\n",
                    cmds[i].description);
        }

        fprintf(fp, "---\n\n");
    }

    fclose(fp);
}

int main() {
    FILE *fp = fopen("README.md", "w");
    if (!fp) {
        perror("Failed to open markdown file");
        return 1;
    }
    generate_input_markdown(fp, commands);
    return 0;
}
