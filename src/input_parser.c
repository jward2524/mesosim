/*
command table
keyword value value value
description


need: 
command
[?] accepted values for the command
description
fxn to perform
*/

struct {
    char name[64];
    char description[64];
    
} param;

struct {
    char keyword[64];
    char parameters[1024];
    char description[1024];
    struct param *parameter_arr;
    int (*command_func)(int tokc, char *tokv);
} cmd;

int cmd_systemsize(int tokc, char *tokv) {

}

struct cmd command_table[] = {
    {
        .keyword = "systemsize";
        .parameters = "u v w";
        .description = "";
        .command_fxn = &cmd_systemsize;
    }
}
