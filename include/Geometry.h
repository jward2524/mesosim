#include "Defs.h"

#pragma once

// [ ]: what is the difference between the atom linked list and the atom array?
typedef struct
	{
		char name[24]; 			// name, 24 characters long
		char type;				// up to 255 atom types

		double coord[3]; 		// working orthogonal coordinates

		// [ ]: why is this a double? when is it ever fractional?
		double lattice[3];		// multiples of unit cell / lattice vectors

		double bsradius;		// atom ball and stick radius // XXX: visualization, and should depend on atom type, not per atom
		//double sfradius;		// atom spacefilling radius

		//char visible;			// Visible to draw/plot routine
		//char selected;			// Selected or not

		//char style;				// 0  Outline
      							// 1  Filled
								// 2  Anisotropic ellipsoid
								// 3  Isotropic ellipsoid

		//double color[3];   	// Atom fill color.  RGB: 0-1 for each color index.

		// simulation data

		// atom lives multiple times on the transition list
		// index refers to the transition vector
		// in the implementation here, the maximum number of times the atom may live on the transition list
		// equals the number of directions it may diffuse, plus one more special transition such as dissolution
		// or transformation into another kind of atom

		int position_on_transition_list[MAXIMUM_NUMBER_OF_NEIGHBORS + DISSOLUTION];		// add 1 for evaporation

		// If occupied_neighbor_sites[i] >= 0, the site at vector offset i is occupied by an atom indexed by the
		// value of occupied_neighbor_sites[i].  There are three special cases
		//		-1: empty
      	//		-2: buried
		//		-3: random
		// position on this list relates to jump vector in jump_offset list
		int occupied_neighbor_sites[MAXIMUM_NUMBER_OF_NEIGHBORS];
		//int int_or_ext[MAXIMUM_NUMBER_OF_NEIGHBORS];

		// linked list to nearby atoms within the same "zone" (macroscopic region of space)
		// by including this information, finding nearby atoms is greatly speeded up.

		int next_atom;
		int previous_atom;
		
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

// command line commands

typedef struct
	{
		char name[25];
		int number;
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
		double lbound[MAXIMUM_NUMBER_OF_ACTIVATION_BARRIERS];
		double ubound[MAXIMUM_NUMBER_OF_ACTIVATION_BARRIERS];
	} Trans_Prob;

typedef struct
	{
		int offset;				// index of the linked list for this zone
	} Zone;

typedef struct
	{
		double k;						// rate constant?
		double frequency;
		int number;           			// number (count?) of this type of transition in transition list
		int offset; 					// index to first item in transition list with this rate constant 
	} Rate;

typedef struct
	{
		int number_in_list;		// location on atom list of atom that the transition belongs to / is acting on
		char offset_vector;   	// this is the translation vector, index in jump offset?

	} Transition_List;

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

/*typedef struct
	{
		int index;				// number of atom
		double normal[3];		// normal to this atom
		double inside_vector[3];
		int neighbors[20];
		int nn;
		int rnn;
		double coord[3];
		double color[3];
		int triangles[20];
		int num_triangles;
		int mark;

		double k1;
		double k2;

		double area;
		double dx[3];
		double newcoord[3];

		double newnormal[3];

		double residual;
		double residual_direction;

	} Surface_Atom;*/ //used in display/plot/curvature
	
/*typedef struct
{
	int index[3];				// surface atoms comprising triangle
	double normal[3][3];		// normals of atoms comprising triangle
	double coord[3][3];
	double overallnormal[3];
	int mark;

	int nn;						// number of bordering triangles
	int neighbors[15];			// indices to bordering triangles

	double area;

	double centroid[3];

	double color[3];

	double kappa1;
	double kappa2;

//	int order[10][3];

	double angles[3];

} Triangle;*/ //used in curvature/display/plot

/*typedef struct
{
	int spoint[750];
} Spoint;*/ //only used in curvature

/*typedef struct
{
	int snum;
	char incarnated;
	Spoint *sindex;
	//sindex[100];
} TZone;*/ //only used in curvature