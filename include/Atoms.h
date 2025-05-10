int add_atom(int x, int y, int z, int type, int special);
int atom_at(int cx, int cy, int cz);
void remove_atom(int at);
void move_atom(int ia, int fa);

// void make_buried_atoms_real(void);

int random_reincarnate_atom(int x, int y, int z, int type, int vc); //feels like it'll never be called
int reincarnate_atom(int x, int y, int z, int type, int vc); //will not be called
void bury_atom(int at, int *pos); //no longer relevant

// rotation and molecule orientation
// these prototypes are wrong, redefined somewhere else
// void rotmat(Atom* atm[], int na, double rtmat[3][3]);
void organize(Atom* atom_arr[], int atom_cnt); //keep for now but need to retool
void orthomol(Atom* atom_arr[], int atom_cnt, double basis[3][3]);
void centerg(Atom* atom_arr[], int atom_cnt); //not really relevant now (visualization)

// general atom and bond handling routines
void copy_atom(int, int);
void create_default_atom(int na); //can modify this to remove things like color?

void kill_atom(int atom_number);
void cell_to_latmat(double c[6], double ltmt[3][3]);
void primitive_basis2ucell_params(double ltmt[3][3], double c[6]);

void cartesian2lattice(double ccart[3], double invert_primitive_basis[3][3], int clattice[3]);
void lattice2cartesian(int clattice[3], double primitive_basis[3][3], double ccart[3]);

/* symmetry related variables */
extern double rmat[3][3];
extern double centroid[3];
extern crystal_offset jump_offset[MAXIMUM_NUMBER_OF_NEIGHBORS];
extern int opposite_offset[MAXIMUM_NUMBER_OF_NEIGHBORS];
extern double primitive_basis[3][3], invert_primitive_basis[3][3];
extern double ucell_params[6];
extern double dax, day, daz;
extern crystal_offset lattice_first_offset[24];
extern crystal_offset lattice_second_offset[24];
extern crystal_offset bcc_offset[8];
extern crystal_offset fcc_offset[12];
extern crystal_offset sc_offset[6];
extern crystal_offset sc_second_offsets[12];
extern crystal_offset fcc_second_offsets[6];
extern crystal_offset bcc_second_offsets[6];