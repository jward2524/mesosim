#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "Defs.h"

// atom array contains all atoms
// atom linked list is multiple shorter lists (zone), split based on spatial location for faster addressing
typedef struct
	{
		char name[24]; 			// name, 24 characters long
		char type;				// up to 255 atom types
		// TODO: remove - only used for printing, so only calculate when printing
		double cart_coord[3]; 		// working [cart]esian (orthogonal) [coord]inates

		int lattice[3];		// multiples of unit cell / lattice vectors

		double bsradius;		// atom ball and stick radius // XXX: visualization, and should depend on atom type, not per atom
		//double sfradius;		// atom spacefilling radius
		// XXX: commented code - related to visualization
		//char visible;			// Visible to draw/plot routine
		//char selected;			// Selected or not

		//char style;				// 0  Outline
      							// 1  Filled
								// 2  Anisotropic ellipsoid
								// 3  Isotropic ellipsoid

		//double color[3];   	// Atom fill color.  RGB: 0-1 for each color index.

		// simulation data

		// atom lives multiple times on the transition list transition_arr
		// index of this list refers to the transition vector in jump_offset[index]
		// in the implementation here, the maximum number of times the atom may live on the transition list
		// equals the number of directions it may diffuse, plus one more special transition such as dissolution
		// or transformation into another kind of atom

		int transition_indices[MAXIMUM_NUMBER_OF_NEIGHBORS + DISSOLUTION];		// add 1 for evaporation

		// If neighbor_atom_idxs[i] >= 0, the site at vector offset i is occupied by an atom indexed by the
		// value of neighbor_atom_idxs[i].  There are three special cases
		//		-1: empty
      	//		-2: buried
		//		-3: random
		// position on this list relates to jump vector in jump_offset list
		int neighbor_atom_idxs[MAXIMUM_NUMBER_OF_NEIGHBORS];
		//int int_or_ext[MAXIMUM_NUMBER_OF_NEIGHBORS];

		// linked list to nearby atoms within the same "zone" (macroscopic region of space)
		// by including this information, finding nearby atoms is greatly speeded up.

		int next_atom;
		int previous_atom;
		// XXX: commented code - related to visualization
		// information about drawing bonds to nearby atoms.  this is purely visual, or "cosmetic" and does not affect
		// the simulation.  In fact, the first thing the simulation will do is kill any existing bonds.

		//int bond[MAXIMUM_NUMBER_OF_COSMETIC_BONDS];			// links to the bond drawing list
		//int nob;											// number_of_bonds

		// things specific to drawing ORTEP-stype x-ray structures

		//double biso;			// isotropic thermal ellipsoid radius
		//double ecos[3][3];		// aniostropic thermal ellipsoid parameters
		//double erms[3];

	} Atom;


// The Bond structure contains information about drawing lines between atoms

typedef struct
	{
		int from;				// Pointer to atom at start of bond
		int to;					// Pointer to atom at end of bond
		double radius;			// Bond radius

		//int style;				// 0 Outline
      							// 1 Filled
								// 2 Dashed
								// 3 Dotted
		//double color[3];
   } Bond;

// Atom_type is used in saving the periodic table.  Also gives default colors for atoms based on first letter (e.g., "O"
// atoms are red.

typedef struct
	{
		char name[3];     	
		double bsradius;
		//double sfradius;
		//double color[3];		
	} Atom_type;

// command line commands // XXX
typedef struct
	{
		char name[25];
		int count;
	} Command;

// the direction of an atomic jump / diffusion move
typedef struct
	{
		int dx;
		int dy;
		int dz;
   } crystal_offset;

typedef struct
	{
		int listnum[MAXIMUM_NUMBER_OF_ACTIVATION_BARRIERS];
		double lbound[MAXIMUM_NUMBER_OF_ACTIVATION_BARRIERS]; // ENHANCE: lbound and ubound are duplicating info
		double ubound[MAXIMUM_NUMBER_OF_ACTIVATION_BARRIERS]; // lbound[i] = ubound[i+1] ?
	} TransProb;

typedef struct
	{
		int offset;				// index of the linked list for this zone
	} Zone;

typedef struct
	{
		double k;						// rate constant
		double frequency;				// rate constant * transition_count, for calculating probability of transition
		int transition_count;           // number (count?) of this type of transition in transition_arr
		int transition_start_idx; 		// index to first item in transition_arr with this rate constant 
		unsigned char *atom_env; 	// can't do variable length arrays, so pointer instead
		// pointer to start of 2D array num_bond_types*idx1 + idx2

		// same format as nnE[][]: atom_env[num_nn_levels][num_bond_types]
		// int env_hash; // hash of atom_env array, for comparison of environments - based on number of distinct environments?
		// can't use memcmp bc atom_env[0] is char*, not char
		int is_evaporation;
	} Rate;

typedef struct
{

} AtomicEnv;

typedef struct
	{
		int atom_idx;		// index in atom_arr of atom that the transition belongs to / is acting on
		unsigned char offset_idx;   	// this represents the translation vector, as the index in jump_offset

	} Transition;

typedef struct
	{
		char initial[MAXIMUM_NUMBER_OF_NEIGHBORS];		// indices of nearby initial neighbors
		char final[MAXIMUM_NUMBER_OF_NEIGHBORS];		// indices of nearby final neighbors
		char in;		// number of neighbors in initial configuration
		char fn;		// number of neighbors in final configuration

		double rate;	

		int previous;
		int next;

	} Configuration;

typedef struct
	{
		double symmat[3][4];
	} Symmetry_Element;

typedef struct
	{
		char name[64];
		int elements;								// number of equivalent positions
		Symmetry_Element symmetry_element[256];		// default to too many symmetry_elements
	} Spacegroup;

typedef struct
{
	double r;
	double g;
	double b;
} Atom_Color;
// XXX: visualization I think, if not its commented code

#endif // GEOMETRY_H
