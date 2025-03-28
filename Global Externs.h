#include "stdafx.h"
#include "Defs.h"
#include "Geometry.h"

/* variables that get used */

//simulation details
extern Atom *atom[];
extern Bond bond[];
extern int nat;
extern int number_rates;
extern double overpotential;
extern int num_sims;
extern double elapsed_time;
extern double default_color[];
extern bool simulation_initialized;
extern bool simulation_should_kill_itself;
extern char atom_names[3][3];
extern int simulation_type;
extern Atom temp_atom;
extern double total_internal_energy;

//random number generation
extern long seed;
extern long idum2;
extern long iy;
extern long iv[NTAB];

//input-output
extern char console_outstring[];
extern char command_string[];
extern char return_message[];
extern char outFile[];
extern char default_extension[];
extern FILE *temp_log;

//crystal symmetry details
extern double rmat[3][3];
extern double cg[3];
extern crystal_offset jump_offset[];
extern int opposite_offset[];

extern crystal_offset bcc_offset[];
extern crystal_offset fcc_offset[];
extern crystal_offset sc_offset[];

extern crystal_offset sc_second_offsets[];
extern crystal_offset fcc_second_offsets[];
extern crystal_offset bcc_second_offsets[];

extern double dax, day, daz;

extern double latmat[3][3], ilatmat[3][3];
extern double cell[];
extern double normal_x, normal_y, normal_z;

/* variables for curvature, display, symmetry, deposition, etc. */
//extern long long current_iteration;
//extern int *zorder[]; //not used?

//extern unsigned int customcolors[16]; //not used
//extern bool simulation_physics_has_changed;
//extern int generic_flag;
//extern unsigned int simabortcode; //just using return values now
//extern char *allstring;
//extern int number_selected;					// Number of atoms selected
//extern int *selected[];		// Indices of selected atoms
//extern bool rotation_notify_flag;
//extern double imat[4][4];
//extern double newlong[3][3];
//extern Atom_type atom_type[];
//extern Atom_type atom_type_proteins[];
//extern char *alphabetall;
//extern char os[];
//extern char command_file[]; //not used anymore?
//extern int command_location;
//extern FILE *pdb_file;
//extern double pdb_cell[3][3];
//extern char view_save_string[];
//extern FILE *view_save_file;
//extern char command_identifier[];
//extern char command_params_string[];
//extern int command;
//extern int command_status;
//extern Command commands[];
//extern double plotlinewidth;
//extern int plot_format;
//extern double taper;
//extern double plot_sfactor;
//extern crystal_offset roughening_offset;
//extern char space_group_name[256];
//extern int sg_cells_in_x;
//extern int sg_cells_in_y;
//extern int sg_cells_in_z;
//extern double symmat[3][4];
//extern Spacegroup spacegroup[300];
//extern int num_symops;
//extern int deposition_type;
//extern BOOL draw_dropping_atom;
//extern double dax0, day0, daz0;
//extern int niod;
//extern Surface_Atom *surface_atoms[MAXIMUM_NUMBER_OF_SURFACE_ATOMS];
//extern int number_triangles;
//extern int need_to_find_avg_curvature;
//extern Triangle *triangle[MAXIMUM_NUMBER_OF_SURFACE_ATOMS];
//extern int atoms_on_surface, original_atoms_on_surface;
//extern int point_surface_atoms;
//extern int number_surface_atoms;
//extern int patience;
//extern double kappamax, kappamin;
//extern double walker[3];
//extern int picture_count;
//extern double total_surface_energy;