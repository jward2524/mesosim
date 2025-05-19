#pragma once

// The definitions in this file primarily relate specifically to the simulation module
// TODO: change to heap memory so not limited by stack
#define MAXIMUM_NUMBER_OF_ATOMS 500000					// maximum number of atoms; each atom is allocated dynamically
//#define MAXIMUM_NUMBER_OF_SURFACE_ATOMS 2000000		// TODO: move this to input file or calculated from system size

#define MAXIMUM_NUMBER_OF_ACTIVATION_BARRIERS 100		//used! // maximum number of distinct activation barriers for different
														// processes.  This number is large so as to accomodate a
														// time dependent activation barrier, such as in dissolution.

// this is (diffusion)+(deposition)

#define MAXIMUM_NUMBER_OF_CONCURRENT_TRANSITIONS ((MAXIMUM_NUMBER_OF_NEIGHBORS+DISSOLUTION)*((unsigned long long int) MAXIMUM_NUMBER_OF_ATOMS))+10

														// The maximum number of transitions going on simulataneously
														// is equal to (the number of atoms) * (the number of things
														// each atom can do)

#define MAXIMUM_NUMBER_OF_BONDS 1000					// total maximum number of bonds that can be drawn.  Reducing
														// this number is a big memory saver.
#define DISSOLUTION 1 //used

// simulation initial geometries

#define SIMULATION_TYPE_UNDEFINED -1 //used by default!
#define SIMULATION_TYPE_FLAT_SHEET 1
#define SIMULATION_TYPE_CLUSTER 2
#define SIMULATION_TYPE_FROM_FILE 3

// the default number of zones is tts

#define DSIMSIZE 256 //used - system size

#define TTS 256 // zone size

#define ZONES_IN_X TTS		
#define ZONES_IN_Y TTS	
#define ZONES_IN_Z TTS

// numerical and fundamental constants

#define kBoltz 8.617e-5  // eV per Kelvin
#define PI 3.141592654
#define RT 1.0//0.70710678	// 1/sqrt(2); gives the nn spacing in fcc to be exactly 1

#define NUMBER_OF_ZONES ZONES_IN_X*ZONES_IN_Y*ZONES_IN_Z

// different types of default atom configurations

#define NORMAL 0
#define SPECIFIED 9 //needed for fileIO
/*#define BONDED_BELOW 1
#define NOT_IN_ONE_DIRECTION 2
#define RANDOM_BELOW 3
#define RANDOM_ALL_AROUND 4
#define RANDOM_INWARD 5
#define RANDOM_SURROUND 6
#define BURIED_INWARD 7
#define NORMAL_NOGO 8
#define BURIED 10*/

// bonding configuration definitions

#define MAXIMUM_NUMBER_OF_NEIGHBORS 12		//(1nns) default is just FCC
#define MAXIMUM_NUMBER_OF_NEIGHBORS_FCC 12			
#define MAXIMUM_NUMBER_OF_NEIGHBORS_BCC 8			
#define MAXIMUM_NUMBER_OF_NEIGHBORS_SC 6			

#define MAXIMUM_NUMBER_OF_NEIGHBORS2 12	//(2nns) default is SC to make it large
#define MAXIMUM_NUMBER_OF_NEIGHBORS2_FCC 6
#define MAXIMUM_NUMBER_OF_NEIGHBORS2_BCC 6		
#define MAXIMUM_NUMBER_OF_NEIGHBORS2_SC 12		

// TODO: make these into an enum
// TODO: then use enum to get features of lattice
#define FCC 1
#define BCC 2
#define SC 3
#define DIAMOND 4

// simple cubic lattice vectors in orthogonal space
// unit is lattice parameters (a)
// TODO: add parameter that defines length of lattice parameter
// **these are not necessarily unit vectors**
#define SCXV1 1.0
#define SCYV1 0.0
#define SCZV1 0.0

#define SCXV2 0.0
#define SCYV2 1.0
#define SCZV2 0.0

#define SCXV3 0.0
#define SCYV3 0.0
#define SCZV3 1.0

// fcc primitive lattice vectors in orthogonal space (xyz axes)
// RT should be 1/2, not 1/sqrt(2) or 1, for primitive cell
#define FCCXV1 1.0*RT
#define FCCYV1 1.0*RT
#define FCCZV1 0.0*RT

#define FCCXV2 1.0*RT
#define FCCYV2 0.0*RT
#define FCCZV2 1.0*RT

#define FCCXV3 0.0*RT
#define FCCYV3 1.0*RT
#define FCCZV3 1.0*RT

// bcc primitve lattice vectors in orthogonal space

#define BCCXV1 0.5
#define BCCYV1 -0.5
#define BCCZV1 -0.5

#define BCCXV2 0.5
#define BCCYV2 0.5
#define BCCZV2 0.5

#define BCCXV3 -0.5
#define BCCYV3 0.5
#define BCCZV3 -0.5

// time data tracking

#define REGULAR_TIME_INTERVALS 1
#define LN_TIME_INTERVALS 2
#define ITERATION_INTERVALS 3
#define TIME_LIST -1
#define ITERATION_LIST -3

#define SIM_END_BY_STIME 1
#define SIM_END_BY_ITERATIONS 2

/* DEFINITIONS FOR THE CONSOLE INTERFACE */

//#define NUMBER_OF_COMMANDS 200

//#define COMMAND_FILE 1
//#define COMMAND_SCREEN 2


#define MAXIMUM_NUMBER_OF_COSMETIC_BONDS 4

//#define IMPURITY_CONCENTRATION 0.0

// simulation defaults

#define DEFAULT_TEMPERATURE 300
#define DEFAULT_BOND_ENERGY 0.15

#define DEFAULT_BOND_ENERGY_AA 0.15
#define DEFAULT_BOND_ENERGY_AB 0.15
#define DEFAULT_BOND_ENERGY_AC 0.15
#define DEFAULT_BOND_ENERGY_BB 0.15
#define DEFAULT_BOND_ENERGY_BC 0.15
#define DEFAULT_BOND_ENERGY_CC 0.15

#define DEFAULT_RAMP_RATE 0.0
#define DEFAULT_OVERPOTENTIAL 1.05
#define DEFAULT_COMPOSITION_A 0.70
#define DEFAULT_COMPOSITION_B 0.30

//#define DEFAULT_DEPOSITION_RATE_OF_A 0.0
//#define DEFAULT_DEPOSITION_RATE_OF_B 0.0
//#define DEFAULT_DEPOSITION_RATE_OF_C 0.0

// atom defaults

#define DEFAULT_ATOM_NAME "New Atom"
#define DEFAULT_BS_RADIUS 0.25
#define DEFAULT_SF_RADIUS 0.75
#define DEFAULT_ATOM_STYLE 1 //FILLED

#define DEFAULT_ATOM_COLOR_R 0.0
#define DEFAULT_ATOM_COLOR_G 0.0
#define DEFAULT_ATOM_COLOR_B 0.0

//got rid of command mnemonics :)

/*error messages (not all used)*/
#define NO_INPUT_ERROR 0
#define ATOM_ARRAY_CLEARED 1
#define BOND_ARRAY_CLEARED 2
#define ATOM_ARRAY_FILLED 3
#define BOND_ARRAY_FILLED 4
#define FILE_COMMAND_IGNORED 5
#define SHOW_COMMAND_IGNORED 6
#define COULDNT_DECODE_ATOM 7
#define NOT_ENOUGH_PARAMS 8
#define TOO_MANY_PARAMS 9
#define INVALID_AXIS 10
#define INVALID_FLOAT 11
#define INVALID_INTEGER 12
#define DUPLICATE_ATOM_NAME 13
#define CANT_CREATE_BOND 14
#define INVALID_SYMMETRY 15
#define ATOM_ARRAY_NOT_CLEARED 16
#define BOND_ARRAY_NOT_CLEARED 17
#define BAD_BOND_RADIUS 18
#define BAD_BOND_TYPE 19
#define BAD_BOND_COLOR 20
#define NO_SUCH_ATOM 21
#define INVALID_RANGE 22
#define ATOM_COMMAND_INVALID 23
#define CANT_FIND_SYMMETRY 24
#define NO_ATOMS_DELETED 25
#define CANNOT_OPEN_FILE 26
#define INVALID_PARAMS 27
#define GENERAL_ERROR 28

//#define NUMBER_OF_ATOM_TYPES 40
//#define NUMBER_OF_ATOM_TYPES_PROTEINS 7

//#define DEPOSITION_TYPE_RAINFALL 1
//#define DEPOSITION_TYPE_RANDOM_WALKER 2
