#include "stdafx.h"//change back when needed
#include "Defs.h" //changed for case sensitivity
#include "Geometry.h" //all good here

// global variables that get used

// ENHANCE: use malloc for simulation-size arrays?

/* simulation related variables */
int simulation_type = SIMULATION_TYPE_UNDEFINED;
bool simulation_initialized = false;
Atom *atom[MAXIMUM_NUMBER_OF_ATOMS];
Bond bond[MAXIMUM_NUMBER_OF_BONDS];
int nat = 0;
double overpotential;
int num_sims;
double elapsed_time;
Atom temp_atom;
double default_color[3] = {0., 0., 0.};
bool simulation_should_kill_itself;
char atom_names[3][3]={"1", "2", "3"};
double total_internal_energy;

/* RNG related variables */
long seed = (long)DEFAULT_SEED;
long idum2;
long iy;
long iv[NTAB];

/* IO related variables */
char console_outstring[512]; //keep for now; change when file writing changes
FILE *view_save_file;
FILE *temp_log;
char command_string[1024];
char return_message[512];
char outFile[260] = ""; //MAX_PATH variable Windows related, default 260
char default_extension[] = "out";

/* symmetry related variables */
double rmat[3][3];
double cg[3]={(double)0.,(double)0.,(double)0.};
crystal_offset jump_offset[MAXIMUM_NUMBER_OF_NEIGHBORS];
int opposite_offset[MAXIMUM_NUMBER_OF_NEIGHBORS];
double latmat[3][3], ilatmat[3][3];
double cell[6]={1.,1.,1.,90.,90.,90.};
double dax, day, daz;
double normal_x, normal_y, normal_z;
crystal_offset lattice_first_offset[24];
crystal_offset lattice_second_offset[24];

crystal_offset bcc_offset[8] = 
		{
		{-1, -1, -1},
		{0, 0, -1},
		{1, 0, 0},
		{-1, 0, 0},
		{0, 1, 0},
		{0, -1, 0},
		{0, 0, 1},
		{1, 1, 1}
		};

crystal_offset fcc_offset[12] = 
		{
		{0, 1, -1},
		{1, 0, -1},
		{0, 0, -1},
		{1, 0, 0},
		{-1, 0, 0},
		{0, 1, 0},
		{0, -1, 0},
		{0, 0, 1},
		{-1, 1, 0},
		{1, -1, 0},
		{-1, 0, 1},
		{0, -1, 1}
};

crystal_offset sc_offset[6] = 
		{
		{0, 0, -1},
		{0, 0, 1},
		{0, 1, 0},
		{0, -1, 0},
		{1, 0, 0},
		{-1, 0, 0}};

crystal_offset sc_second_offsets[12] = 
		{{1, 1, 0},
		{1, -1, 0},
		{-1, 1, 0},
		{-1, -1, 0},
		{1, 0, 1},
		{-1, 0, 1},
		{1, 0, -1},
		{-1, 0, -1},
		{0, 1, 1},
		{0, -1, 1},
		{0, 1, -1},
		{0, -1, -1}};

crystal_offset fcc_second_offsets[6] = 
		{{1, -1, 1},
		{-1, 1, 1},
		{-1, -1, 1},
		{-1, 1, -1},
		{1, 1, -1},
		{1, -1, -1}};

crystal_offset bcc_second_offsets[6] = 
		{{1, 0, 1},
		{-1, 0, -1},
		{1, 1, 0},
		{-1, -1, 0},
		{0, 1, 1},
		{0, -1, -1}};

//unused variables
//double imat[4][4]={{1.,0.,0.,0.},{0.,1.,0.,0.},{0.,0.,1.,0.},{0.,0.,0.,1.}};
//double newlong[3][3] = {1.,0.,0.,0.,1.,0.,0.,0.,1.};
//crystal_offset roughening_offset;
//bool simulation_physics_has_changed = false; //changed from Globals
//int *zorder[MAXIMUM_NUMBER_OF_ATOMS]; //not used?
//int generic_flag = 0;
//long long current_iteration; //no longer needed
//int number_rates;
//unsigned int customcolors[16]; //not used
//unsigned int simabortcode;
//BOOL pause_simulation = FALSE;
//BOOL simulation_paused = FALSE;
//BOOL simulation_is_inside_loop = FALSE;
//char *allstring = "all"; //keep this it seems important

// selection and rotation variables
//int number_selected = 0;					// Number of atoms selected
//int *selected[MAXIMUM_NUMBER_OF_ATOMS];		// Indices of selected atoms
//bool rotation_notify_flag = false;  //not used but in comments

// periodic table stuff
/*Atom_type atom_type[NUMBER_OF_ATOM_TYPES]=
        {
			{"C",0.25,1.7,0.,0.,0.},
			{"H",0.1,1.2,0.,0.,0.},
			{"S",0.25,1.8,1.,1.,0.},
			{"F",0.25,1.47,0.,1.,0.},
			{"O",0.25,1.52,1.,0.,0.},
			{"N",0.25,1.55,0.,0.,1.},
			{"P",0.25,1.80,0.,1.,0.},
			{"Ag",0.25,1.72,0.,1.,0.},
			{"Ar",0.25,1.88,0.,1.,0.},
			{"As",0.25,1.85,0.,1.,0.},
			{"Au",0.25,1.66,0.,1.,0.},
			{"Br",0.25,1.85,0.,1.,0.},
			{"Cd",0.25,1.58,0.,1.,0.},
			{"Cl",0.25,1.75,0.,1.,0.},
			{"Cu",0.25,1.40,0.,1.,0.},
			{"F",0.25,1.47,0.,1.,0.},
			{"Ga",0.25,1.87,0.,1.,0.},
			{"He",0.25,1.40,0.,1.,0.},
			{"Hg",0.25,1.55,0.,1.,0.},
			{"I",0.25,1.98,0.,1.,0.},
			{"In",0.25,1.93,0.,1.,0.},
			{"K",0.25,2.75,0.,1.,0.},
			{"Kr",0.25,2.02,0.,1.,0.},
			{"Li",0.25,1.82,0.,1.,0.},
			{"Mg",0.25,1.73,0.,1.,0.},
			{"N",0.25,1.55,0.,1.,0.},
			{"Na",0.25,2.27,0.,1.,0.},
			{"Ne",0.25,1.54,0.,1.,0.},
			{"Ni",0.25,1.63,0.,1.,0.},
			{"Pb",0.25,2.02,0.,1.,0.},
			{"Pd",0.25,1.63,0.,1.,0.},
			{"Pt",0.25,1.72,0.,1.,0.},
			{"Se",0.25,1.90,0.,1.,0.},
			{"Si",0.25,2.10,0.,1.,0.},
			{"Sn",0.25,2.17,0.,1.,0.},
			{"Te",0.25,2.06,0.,1.,0.},
			{"Tl",0.25,1.96,0.,1.,0.},
			{"U",0.25,1.86,0.,1.,0.},
			{"Xe",0.25,2.16,0.,1.,0.},
			{"Zn",0.25,1.39,0.,1.,0.}
		};*/
/*Atom_type atom_type_proteins[NUMBER_OF_ATOM_TYPES_PROTEINS]=
        {
			{"C",0.25,1.7,0.,0.,0.},
			{"H",0.1,1.2,0.,0.,0.},
			{"S",0.25,1.8,1.,1.,0.},
			{"F",0.25,1.47,0.,1.,0.},
			{"O",0.25,1.52,1.,0.,0.},
			{"N",0.25,1.55,0.,0.,1.},
			{"P",0.25,1.80,0.,1.,0.}
		};*/

// command line - maybe don't get rid of these for now
//char *alphabetall="1234567890=`~!@#$%^&*()_+[]\\{}|;':\",./<>?abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
//char os[512];
//char command_file[256];
//int command_location = COMMAND_SCREEN;
//FILE *pdb_file;
//double pdb_cell[3][3]={{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};
//char view_save_string[256] = "defview.xyz";
//char command_identifier[10];
//char command_params_string[1024];

//int command;
//int command_status;

//Command is defined in geometry; don't need to change
/*Command commands[NUMBER_OF_COMMANDS]=
	{{"",0},{"",1},{"",2},{"",3},{"",4},{"",5},{"",6},{"",7},
	{"",8},{"",9},{"cell",10},{"atom",11},{"file",12},{"rnge",13},
	{"list",14},{"geom",15},{"axi*",16},{"axis",17},{"bond",18},
    {"plan",19},{"rotr",20},{"save",21},{"plot",22},{"show",23},
    {"scal",24},{"size",25},{"help",26},{"curr",27},{"dist",28},
    {"angl",29},{"trns",30},{"clra",31},{"clrb",32},{"symm",33},
    {"bild",65},{"halt",36},{"quit",36},{"stop",36},{"disp",38},
    {"dele",39},{"titl",40},{"pers",41},{"link",42},{"bndp",43},
    {"slab",44},{"wind",45},{"end",36},{"bij",47}, {"title",40},
    {"scale",24},{"atype",48},{"btype", 49},{"boxon", 50},{"boxoff", 51},
	{"box", 52},{"rmat",53},{"cmat",54}, {"vatom", 55},{"delbond",56},
	{"arad",57},{"smartbond", 58}, {"cryst1", 10}, {"pdb", 59},{"alabel", 60},
	{"ellipsoid", 61},{"blabel", 62}, {"lbond", 63}, {"setsymm", 64},
	{"build", 65}, {"group", 66}, {"chaff", 67},{"range",13},{"geometry",15},
	{"plane",19},{"rotate",20},{"scale",24},{"current",27},{"distance",28},
	{"distances",28},{"angle",29},{"angles",29},{"torsion",30},
	{"clearatoms",31},{"clearbonds",32},{"symmetry",33},{"exit", 36},
	{"display",38},{"delete",39},{"del",39},{"lineplot", 68},{"lplot",68},
	{"taper",69},{"lwidth",70},{"atomtype",48},{"bondtype",49},{"linewidth",70},
	{"scale1",71},{"scale2",72},{"scale3",73}, {"select", 74}, {"simu", 75},
    {"pmode", 76},{"dmode", 77}, {"nat", 78}};*/

// plotting
//double plotlinewidth = 0.005;
//int plot_format = 0;			// 0 is bonds and stick, 1 is space filling, 2 is line and point
//double taper = 0.3;
//double plot_sfactor = 0.;							// scale factor for plotting

//char space_group_name[256] = "P1";
//int sg_cells_in_x = 1;
//int sg_cells_in_y = 1;
//int sg_cells_in_z = 1;

//double symmat[3][4]={{1.,0.,0.,0.},{0.,1.,0.,0.},{0.,0.,1.,0.}};

//int num_symops;
//Spacegroup spacegroup[300];

//int deposition_type = DEPOSITION_TYPE_RAINFALL;
//DEPOSITION_TYPE_RANDOM_WALKER;

//bool draw_dropping_atom = FALSE; honestly doesn't look like it's used

//double dax0, day0, daz0;

//double kappamax, kappamin;

//int niod;

//Surface_Atom *surface_atoms[MAXIMUM_NUMBER_OF_SURFACE_ATOMS];

//int number_triangles = 0;
//int need_to_find_avg_curvature = 0;

//int atoms_on_surface = 0;
//int original_atoms_on_surface = 0;

//int point_surface_atoms = 0;

//int number_surface_atoms = 0;

//Triangle *triangle[2*MAXIMUM_NUMBER_OF_SURFACE_ATOMS];

//int patience = 0;

//double walker[3];

//int picture_count;

//double total_surface_energy;
