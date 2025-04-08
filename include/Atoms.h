int add_atom(double x, double y, double z, int type, int special);
int atom_at(double cx, double cy, double cz);
void remove_atom(int at);
void move_atom(int ia, int fa);

// void make_buried_atoms_real(void);

int random_reincarnate_atom(double x, double y, double z, int type, int vc); //feels like it'll never be called
int reincarnate_atom(double x, double y, double z, int type, int vc); //will not be called
void bury_atom(int at, int *pos); //no longer relevant

// rotation and molecule orientation
// these prototypes are wrong, redefined somewhere else
// void rotmat(Atom* atm[], int na, double rtmat[3][3]);
void organize(Atom* atm[], int n); //keep for now but need to retool
void orthomol(Atom* atm[], int na, double com[3][3]);
void centerg(Atom* atm[], int na); //not really relevant now

// general atom and bond handling routines
void copy_atom(int, int);
void create_default_atom(int na); //can modify this to remove things like color?

void kill_atom(int atom_number);
void cell_to_latmat(double c[6], double ltmt[3][3]);
void latmat_to_cell(double ltmt[3][3], double c[6]);
