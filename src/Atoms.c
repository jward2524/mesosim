#include "Atoms.h"
#include "Vector.h"
#include "Random.h"
#include "Simulation.h"
#include "Simulation_Aux.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

// [ ]: atoms, zones, orientation

/* symmetry related variables */
double rmat[3][3]; // visualization?
double centroid[3]={(double)0.,(double)0.,(double)0.}; // coordinates for center of gravity
crystal_offset jump_offset[MAXIMUM_NUMBER_OF_NEIGHBORS]; // possible atom jumps in simulation
int opposite_offset[MAXIMUM_NUMBER_OF_NEIGHBORS]; // index in jump_offset that has the jump in the opposite direction in simulation; opposite_offset[0]=11 means the opposite direction of jump_offset[0] is jump_offset[11]

// primitive unit cell basis vectors + inverted; primitive_basis[*][0] = basis1, primitive_basis[0][*] = x component of basises
// is also transformation matrix for [lattice to cartesian coordinates] [cartesian to lattice coordinates] respectively
double primitive_basis[3][3], invert_primitive_basis[3][3]; 
double ucell_params[6]={1.,1.,1.,90.,90.,90.}; // unit cell parameters; a b c alpha beta gamma // TODO: this is never used except for printing?
double dax, day, daz;
crystal_offset lattice_first_offset[24];
crystal_offset lattice_second_offset[24];

// Atom_Color atom_color[10];

// direction of possible atom jumps for each crystal lattice type
const crystal_offset BCC_OFFSET[8] = 
		{ // not normalized, in lattice coordinates
		{-1, -1, -1},
		{0, 0, -1},
		{1, 0, 0},
		{-1, 0, 0},
		{0, 1, 0},
		{0, -1, 0},
		{0, 0, 1},
		{1, 1, 1}
		};

const crystal_offset FCC_OFFSET[12] = 
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

const crystal_offset SC_OFFSET[6] = 
		{
		{0, 0, -1},
		{0, 0, 1},
		{0, 1, 0},
		{0, -1, 0},
		{1, 0, 0},
		{-1, 0, 0}};

const crystal_offset SC_OFFSET_2[12] = 
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

const crystal_offset FCC_OFFSET_2[6] = 
		{{1, -1, 1},
		{-1, 1, 1},
		{-1, -1, 1},
		{-1, 1, -1},
		{1, 1, -1},
		{1, -1, -1}};

const crystal_offset BCC_OFFSET_2[6] = 
		{{1, 0, 1},
		{-1, 0, -1},
		{1, 1, 0},
		{-1, -1, 0},
		{0, 1, 1},
		{0, -1, -1}};

// updates atom_arr // [ ]: but doesn't update atom_cnt?
void create_default_atom(int n, Atom** atom_arr) // n = position on atom list
{
	int i,j;
	char errorstring[256]; // XXX: unused

	atom_arr[n] = (Atom *)malloc(sizeof(Atom));

	if (atom_arr[n] == NULL)
	{
		// TODO: free mallocs before exiting
		printf("ERROR! Not enough memory to allocate atom %d\n", n);
		fprintf(stderr, "Couldn't allocate memory for atom %d: %s", n, strerror(errno));
        exit(errno);
	}

	strcpy(atom_arr[n]->name, DEFAULT_ATOM_NAME);
	atom_arr[n]->type = 1;
		
	for (i=0;i<3;++i)
	{
		atom_arr[n]->cart_coord[i] = 0.0;
		atom_arr[n]->lattice[i] = 0.0;
	}

	atom_arr[n]->bsradius = DEFAULT_BS_RADIUS; //set or optional // XXX: vis + commented code
	//atom_arr[n]->sfradius = DEFAULT_SF_RADIUS; //set or optional
		
	//atom_arr[n]->visible = true; //can remove
	//atom_arr[n]->selected = false; //can remove

	//atom_arr[n]->style = DEFAULT_ATOM_STYLE; //can remove

	/*atom_arr[n]->color[0] = DEFAULT_ATOM_COLOR_R; //can remove
	atom_arr[n]->color[1] = DEFAULT_ATOM_COLOR_G; //can remove
	atom_arr[n]->color[2] = DEFAULT_ATOM_COLOR_B; //can remove*/

	for (i=0;i<MAXIMUM_NUMBER_OF_NEIGHBORS + DISSOLUTION;++i)
		atom_arr[n]->transition_indices[i] = -1;

	for (i=0;i<MAXIMUM_NUMBER_OF_NEIGHBORS;++i)
		atom_arr[n]->occupied_neighbor_sites[i] = -1;

	// linked list structure // [ ]: but why

	atom_arr[n]->next_atom = -1;
	atom_arr[n]->previous_atom = -1;
	
	// bonding // XXX: commented code
	/*for (i=0;i<MAXIMUM_NUMBER_OF_COSMETIC_BONDS;++i) //can remove this part?
		atom_arr[n]->bond[i] = -1;								// links to the bond drawing list
	atom_arr[n]->nob = 0;*/										// number_of_bonds

	// things specific to x-ray structures

	/*atom_arr[n]->biso = 0.0; //can remove all of this
	for (i=0;i<3;++i)
		for (j=0;j<3;++j)
			atom_arr[n]->ecos[i][j] = 0.0;
	for (i=0;i<3;++i)
		atom_arr[n]->erms[i] = 0.0;*/

	return;
}

/*******************************************************************************
*******************************************************************************/
// updates [atom_arr], atom_arr[i]->lattice, [atom_arr[n]->cart_coord], atom_cnt, zone_arr[xzone][yzone][zzone].offset, atom_arr[pos]->next_atom/prev_atom, atom_arr[pos]->transition_indices; returns index in atom_arr
int add_atom(int x, int y, int z, int type, int special, struct SimulationState* ss, struct SimulationEnv* se) // lattice coordinates xyz, atom type, special atom conditions (unused)
{ // XXX: special isn't really used
	/*if (x > 60)
		printf("made it in!\n");*/
	int i, j, k, m, n1;

	int xzone, yzone, zzone;
	int checkx, checky, checkz; // position of potential move

	int pos, ct, n2; // position in atom array, presumably; // XXX: ct n2 are unused

	double sp[3], spo[3];
	// [ ]: this is a sanity check? iterating over atom list instead of zone (like atom_at)
	if (atom_at(x,y,z, ss->atom_arr, ss->zone_arr, se) >= 0)
	{
		int num_overlapping = 0;
		for (i=0; i < ss->atom_cnt;++i)
		{
			if ((ss->atom_arr[i]->lattice[0] == x)&&(ss->atom_arr[i]->lattice[1] == y)&&(ss->atom_arr[i]->lattice[2] == z))
				++num_overlapping;
		}

		printf("ERROR! Unable to add atom %d; %d other atoms found at (%lf, %lf, %lf)\n", ss->atom_cnt, num_overlapping, x, y, z);
		ss->simulation_should_kill_itself = true;
		return ss->atom_cnt;
	}

	// allocate memory pointed to by the last element of the atom list
	/*if (x > 60)
		printf("atom of type %d being added at %lf %lf %lf\n", type, x, y, z);*/
	pos = ss->atom_cnt; // position in atom array, presumably
	if (pos > se->max_atoms)
	{
		fprintf(stderr, "More atoms (%d) than allocated in atom array (%d)\n", pos, se->max_atoms);
		exit(1);
	}
	create_default_atom(ss->atom_cnt, ss->atom_arr);

	/*if (x > 60)
		printf("made it past making a default: atom_cnt = %d\n", atom_cnt+1);*/
	// set atom color by type
		
	/*atom[atom_cnt]->color[0] = atom_color[type].r;
	atom[atom_cnt]->color[1] = atom_color[type].g;
	atom[atom_cnt]->color[2] = atom_color[type].b;*/

	++ss->atom_cnt;

	findzone(&xzone, &yzone, &zzone, x, y, z, se); // TODO: this is already done in atom_at - why repeat it
	// XXX: commended code
	/*if (x > 60)
		printf("found zone\n");*/
	// xzone, yzone, zzone now have a position open at the end of the zone
	// pos points to this location.  mark the spot and increment the number of atoms
	// in the zone.
	
	// update the zone.  Increment the number of elements.  If the zone was
	// empty, create a link to the first element in that zone

	if (ss->zone_arr[xzone][yzone][zzone].offset == -1)	// first atom in [zone?] linked list
	{
		ss->zone_arr[xzone][yzone][zzone].offset = pos;

		ss->atom_arr[pos]->next_atom = -1;	// no valid link
		ss->atom_arr[pos]->previous_atom = -1;
	}
	else
	{
		// link this atom to the others in the zone linked list

		j = ss->zone_arr[xzone][yzone][zzone].offset;				// first element of list

		while (ss->atom_arr[j]->next_atom != -1)
			j = ss->atom_arr[j]->next_atom;

		// j points to the previous last atom in the zone linked list and points to nothing

		ss->atom_arr[j]->next_atom = pos;
		ss->atom_arr[pos]->previous_atom = j;
		ss->atom_arr[pos]->next_atom = -1;
	}

	/*if (x > 60) // XXX: commented prints
		printf("other zone logic done\n");*/

	ss->atom_arr[pos]->lattice[0] = x;
	ss->atom_arr[pos]->lattice[1] = y;
	ss->atom_arr[pos]->lattice[2] = z;

	ss->atom_arr[pos]->type = type;
	strcpy(ss->atom_arr[ss->atom_cnt-1]->name, se->atom_names[type-1]); // TODO: use pos instead of atom_cnt-1

	/*if (x > 60) // XXX: commented prints
		printf("copied the name: atom is type %s\n", atom[atom_cnt-1]->name);*/

	// find (or set) the occupied neighbor sites

	// [ ]: saturate all the bonds, except it doesn't?

	for (i=0; i < se->max_neighbors; ++i)
	{
		// mark that this atom cannot yet jump in direction i
		/*if (x > 60) // XXX: commented prints
			printf("i am testing neighbor %d\n", i);*/

		ss->atom_arr[pos]->transition_indices[i] = -1;
		// [ ]: if system size is still 1 (in z direction?) and jump isn't zero, skip it?
		// if ((ssz == 1)&&(jump_offset[i].dz != 0))
		// {
		// 	/*if (x > 60) // XXX: commented prints
		// 		printf("met the corner case\n");*/
		// 	continue;
		// }

		// checkx, checky, checkz point to the neighboring site
		// so update occupied neighbor site according to atom_at(checkx, checky, checkz);

		checkx = x + jump_offset[i].dx;
		checky = y + jump_offset[i].dy;
		checkz = z + jump_offset[i].dz;
		/*if (x > 60) // XXX: commented prints
			printf("before pbc xyz %lf %lf %lf\n", checkx, checky, checkz);*/
		adjust_pbc(&checkx, &checky, &checkz, se);
		/*if (x > 60)
			printf("after pbc xyz %lf %lf %lf\n", checkx, checky, checkz);*/

		switch(special)
		{
			case NORMAL:								// normal bonding considerations
			//case NORMAL_NOGO:							//remainder of cases become irrelevant when removing burial
				// set occupied_neighbor_site[i] to the atom at that site.
				// If there really is an atom there, cross-link it to our new atom.
				/*if (x > 60) {
					printf("I'm in the normal case\n");
					printf("atom at site %lf %lf %lf is number %d\n", checkx, checky, checkz, atom_at(checkx, checky, checkz));
				}*/
				ss->atom_arr[pos]->occupied_neighbor_sites[i] = atom_at(checkx, checky, checkz, ss->atom_arr, ss->zone_arr, se);

				/*if (x>60) // XXX: commented prints
					printf("pos = %d, i = %d, atom[pos]->occupied[i] = %d\n", pos, i, atom[pos]->occupied_neighbor_sites[i]);*/
				// if atom is present at potential jump site, fill position in occupied_neighbor_sites of this atom and the found neighbor atom
				if (ss->atom_arr[pos]->occupied_neighbor_sites[i] >= 0 ) {
					/*if (x>60)
						printf("reverse: index %d, opposite offset %d, pos %d\n", atom[pos]->occupied_neighbor_sites[i], opposite_offset[i], pos);*/
					ss->atom_arr[ss->atom_arr[pos]->occupied_neighbor_sites[i]]->occupied_neighbor_sites[opposite_offset[i]] = pos;
					/*if (x>60)
						printf("uno reverse didn't break me\n");*/
				}

				/*if (x > 60)
					printf("About to leave the normal case\n");*/
				break;
			// XXX: commented code
			/*case RANDOM_SURROUND:								// normal bonding considerations
				// set occupied_neighbor_site[i] to the atom at that site.
				// If there really is an atom there, cross-link it to our new atom.
				atom[pos]->occupied_neighbor_sites[i] = atom_at(checkx, checky, checkz);

				if (atom[pos]->occupied_neighbor_sites[i] >= 0 )
					atom[atom[pos]->occupied_neighbor_sites[i]]->occupied_neighbor_sites[opposite_offset[i]] = pos;
				else
					atom[pos]->occupied_neighbor_sites[i] = -3;
				break;

			case BONDED_BELOW:						// bond to bulk - all bonds downward are to buried atoms.
				atom[pos]->occupied_neighbor_sites[i] = atom_at(checkx, checky, checkz);

				if (atom[pos]->occupied_neighbor_sites[i] >= 0 )
					atom[atom[pos]->occupied_neighbor_sites[i]]->occupied_neighbor_sites[opposite_offset[i]] = pos;
				else
					if (jump_offset[i].dz < 0)
						atom[pos]->occupied_neighbor_sites[i] = -2;
				break;

			case BURIED:
				atom[pos]->occupied_neighbor_sites[i] = -2;
				break;

			case RANDOM_BELOW:						// bond to bulk - all bonds downward are to buried atoms.
				atom[pos]->occupied_neighbor_sites[i] = atom_at(checkx, checky, checkz);

				if (atom[pos]->occupied_neighbor_sites[i] >= 0 )
					atom[atom[pos]->occupied_neighbor_sites[i]]->occupied_neighbor_sites[opposite_offset[i]] = pos;
				else if (jump_offset[i].dz < 0)
					atom[pos]->occupied_neighbor_sites[i] = -3;
				break;

			case RANDOM_ALL_AROUND:						// bond to bulk - all bonds downward are to buried atoms.
				if (i == opposite_offset[niod])
				{
					atom[pos]->occupied_neighbor_sites[i] = atom_at(checkx, checky, checkz);

					if (atom[pos]->occupied_neighbor_sites[i] >= 0 )
						atom[atom[pos]->occupied_neighbor_sites[i]]->occupied_neighbor_sites[opposite_offset[i]] = pos;
					else
						atom[pos]->occupied_neighbor_sites[i] = -1;
				}
				else
				{
					atom[pos]->occupied_neighbor_sites[i] = atom_at(checkx, checky, checkz);

					if (atom[pos]->occupied_neighbor_sites[i] >= 0 )
						atom[atom[pos]->occupied_neighbor_sites[i]]->occupied_neighbor_sites[opposite_offset[i]] = pos;
					else
						atom[pos]->occupied_neighbor_sites[i] = -3;
				}

				break;

			case NOT_IN_ONE_DIRECTION:	  	// bond normally to what's actually nearby, -2 otherwise, except -1 in the no_bond_direction
				atom[pos]->occupied_neighbor_sites[i] = atom_at(checkx, checky, checkz);

				if (atom[pos]->occupied_neighbor_sites[i] >= 0)
					atom[atom[pos]->occupied_neighbor_sites[i]]->occupied_neighbor_sites[opposite_offset[i]] = pos;
				else if (i != no_bond_direction)
					atom[pos]->occupied_neighbor_sites[i] = -3;
				break;

			case RANDOM_INWARD:
				sp[0] = jump_offset[i].dx;
				sp[1] = jump_offset[i].dy;
				sp[2] = jump_offset[i].dz;
					
				vecmul(sp, primitive_basis, spo);
				unit(spo, spo);

				if (fdot(spo, a_pn) > 0.5)
				{
					// normal points inward, make buried

					if ((atom[pos]->occupied_neighbor_sites[i] = atom_at(checkx, checky, checkz)) >= 0)
						atom[atom[pos]->occupied_neighbor_sites[i]]->occupied_neighbor_sites[opposite_offset[i]] = pos;
					else
						atom[pos]->occupied_neighbor_sites[i] = -3;
					break;
				}
				else
				{
					// normal points outward; make NORMAL

					atom[pos]->occupied_neighbor_sites[i] = atom_at(checkx, checky, checkz);
	
					if (atom[pos]->occupied_neighbor_sites[i] >= 0 )
						atom[atom[pos]->occupied_neighbor_sites[i]]->occupied_neighbor_sites[opposite_offset[i]] = pos;
					break;
				}

			case BURIED_INWARD:
				sp[0] = jump_offset[i].dx;
				sp[1] = jump_offset[i].dy;
				sp[2] = jump_offset[i].dz;
					
				vecmul(sp, primitive_basis, spo);
				unit(spo, spo);

				if (fdot(spo, a_pn) > 0.5)
				{
					// normal points inward, make buried

					if ((atom[pos]->occupied_neighbor_sites[i] = atom_at(checkx, checky, checkz)) >= 0)
						atom[atom[pos]->occupied_neighbor_sites[i]]->occupied_neighbor_sites[opposite_offset[i]] = pos;
					else
						atom[pos]->occupied_neighbor_sites[i] = -2;
					break;
				}
				else
				{
					// normal points outward; make NORMAL

					atom[pos]->occupied_neighbor_sites[i] = atom_at(checkx, checky, checkz);
	
					if (atom[pos]->occupied_neighbor_sites[i] >= 0 )
						atom[atom[pos]->occupied_neighbor_sites[i]]->occupied_neighbor_sites[opposite_offset[i]] = pos;
					break;
				}
			*/
			case SPECIFIED: //only used when reading kmc files?
				printf("I'm in specified\n");
				// ss->atom_arr[pos]->occupied_neighbor_sites[i] = temp_atom.occupied_neighbor_sites[i];
				break;
			default:
				//printf("I made it to the default????\n");
				break;
		}
	}

	//printf("bond saturated\n");
	// can't evaporate either
	// [ ]: bc bonds saturated, except they aren't? so don't let it evaporate
	ss->atom_arr[pos]->transition_indices[se->max_neighbors] = -1;

	//printf("did you cause a problem\n");
	// set the hopping rates for this atom

	if (special != NORMAL)		// rates will be refreshed soon
		return pos;

	// now set the transition rates
	// ENHANCE: remove this call from function - instead, do it manually when needed (don't need to call after initializing every atom for the first time)
	refresh_transitions(pos, ss, se); // [ ]: why doing this now, after every atom, when neighbors haven't been created yet?

	//printf("my transition refreshed\n");
	// cycle through the nearest neighbors, refresh their transitions [or bury as necessary]

	for (i=0; i < se->max_neighbors; ++i)
	{
		//printf("trying to refresh neighbor %d\n", i);
		j = ss->atom_arr[pos]->occupied_neighbor_sites[i];

		if (j >= 0)
		{
			k = refresh_transitions(j, ss, se);	// refresh transitions of neighbor
			//we don't want to bury our atoms: this is removed from this point
		}
	}

	//printf("neighbor transition refreshed\n");
	return pos;
}


/********************************************************************************/
/********************************************************************************/

// copies one element ia of atom list to another point fa
// used only within the simulation routines, things like bonds, etc. are not copied.

void move_atom(int ia, int fa, Atom** atom_arr, Zone zone_arr[ZONES_IN_X][ZONES_IN_Y][ZONES_IN_Z], Transition** transition_arr, struct SimulationEnv *se)
{
	int n;
	int n2;
	int xzone, yzone, zzone;

	atom_arr[fa]->lattice[0] = atom_arr[ia]->lattice[0];
	atom_arr[fa]->lattice[1] = atom_arr[ia]->lattice[1];
	atom_arr[fa]->lattice[2] = atom_arr[ia]->lattice[2];

	/*atom[fa]->color[0] = atom[ia]->color[0];
	atom[fa]->color[1] = atom[ia]->color[1];
	atom[fa]->color[2] = atom[ia]->color[2];*/

	atom_arr[fa]->bsradius = atom_arr[ia]->bsradius;

	atom_arr[fa]->type = atom_arr[ia]->type;
	strcpy(atom_arr[fa]->name, atom_arr[ia]->name);

	//atom[fa]->biso = atom[ia]->biso;

	atom_arr[fa]->next_atom = atom_arr[ia]->next_atom;
	atom_arr[fa]->previous_atom = atom_arr[ia]->previous_atom;

	if (atom_arr[fa]->next_atom >= 0)							// update valid link
		atom_arr[atom_arr[fa]->next_atom]->previous_atom = fa;

	if (atom_arr[fa]->previous_atom >= 0)
		atom_arr[atom_arr[fa]->previous_atom]->next_atom = fa;
	else													// fa is first element of a zone
		{
			findzone(&xzone, &yzone, &zzone, atom_arr[fa]->lattice[0], atom_arr[fa]->lattice[1], atom_arr[fa]->lattice[2], se);
			zone_arr[xzone][yzone][zzone].offset = fa;
		}

	for (n=0;n < se->max_neighbors;++n)
		{
			atom_arr[fa]->occupied_neighbor_sites[n] =	atom_arr[ia]->occupied_neighbor_sites[n];

			n2 = atom_arr[ia]->occupied_neighbor_sites[n];
			if (n2 >= 0) atom_arr[n2]->occupied_neighbor_sites[opposite_offset[n]] = fa;

			atom_arr[fa]->transition_indices[n] = atom_arr[ia]->transition_indices[n];

			if (atom_arr[fa]->transition_indices[n] >= 0)
			transition_arr[atom_arr[fa]->transition_indices[n]]->atom_idx = fa;
		}

	atom_arr[fa]->transition_indices[se->max_neighbors]
		= atom_arr[ia]->transition_indices[se->max_neighbors];

	if (atom_arr[fa]->transition_indices[se->max_neighbors] >= 0){
		transition_arr[atom_arr[fa]->transition_indices[se->max_neighbors]]->atom_idx = fa;
	}
		

	return;
}

/******************************************************************************/
/******************************************************************************/

void remove_atom(int at, struct SimulationState* ss, struct SimulationEnv* se)
{

	int i,j;
	int type;
	int xzone, yzone, zzone;

	int number_of_new_atoms, number_of_new_random_atoms;
	int x, y, z;
	int vc;

	double subv;

	struct {
		int x;
		int y;
		int z;
		int vc;
	} new_atom[MAXIMUM_NUMBER_OF_NEIGHBORS], new_random_atom[MAXIMUM_NUMBER_OF_NEIGHBORS];

	int nt[MAXIMUM_NUMBER_OF_NEIGHBORS];
	int nnt = 0;

	int nb[MAXIMUM_NUMBER_OF_NEIGHBORS];
	int nnb = 0;

	int nr[MAXIMUM_NUMBER_OF_NEIGHBORS];
	int nnr = 0;

	number_of_new_atoms = 0;
	number_of_new_random_atoms = 0;

	// destroy cross-references to neighboring atoms.  Mark this spot as empty

  	type = ss->atom_arr[at]->type;

	for (i=0; i < se->max_neighbors; ++i)
	{
		j = ss->atom_arr[at]->occupied_neighbor_sites[i];

		switch(j) //might be irrelevant if burial removed // [ ]: burried
		{
			case -2:
				// re-incarnate the buried atom.  this atom will be of type "type"
				// we'll have to add an atom at this point, but we'll do this
				// only after we remove the existence of the current atom
							

				new_atom[number_of_new_atoms].x = ss->atom_arr[at]->lattice[0] + jump_offset[i].dx;
				new_atom[number_of_new_atoms].y = ss->atom_arr[at]->lattice[1] + jump_offset[i].dy;
				new_atom[number_of_new_atoms].z = ss->atom_arr[at]->lattice[2] + jump_offset[i].dz;

				new_atom[number_of_new_atoms].vc = opposite_offset[i];	// only allowed direction

				adjust_pbc(&new_atom[number_of_new_atoms].x,
							&new_atom[number_of_new_atoms].y,
							&new_atom[number_of_new_atoms].z, se);

				++number_of_new_atoms;
				break;

			case -3:
				// incarnate a random atom

				new_random_atom[number_of_new_random_atoms].x = ss->atom_arr[at]->lattice[0] + jump_offset[i].dx;
				new_random_atom[number_of_new_random_atoms].y = ss->atom_arr[at]->lattice[1] + jump_offset[i].dy;
				new_random_atom[number_of_new_random_atoms].z = ss->atom_arr[at]->lattice[2] + jump_offset[i].dz;

				new_random_atom[number_of_new_random_atoms].vc = opposite_offset[i];	// only allowed direction

				adjust_pbc(&new_random_atom[number_of_new_random_atoms].x,
                     		&new_random_atom[number_of_new_random_atoms].y,
							&new_random_atom[number_of_new_random_atoms].z,
							se);

				++number_of_new_random_atoms;
				break;

			case -1:
				if (ss->atom_arr[at]->transition_indices[i] != -1)		// i.e., this is a spot to jump to 
					take_off_transition_list(at, i, ss);
				break;

			default:	// make other atom see this spot as empty
				ss->atom_arr[j]->occupied_neighbor_sites[opposite_offset[i]] = -1;
				nt[nnt] = j;
				++nnt;
				break;
		}
	}

	take_off_transition_list(at, se->max_neighbors, ss);			// dissolution

	// now get rid of the atom.  This is almost equivalent to burying it.
	// find out what zone we're in

	// remove the atom from the atom list

	i = ss->atom_arr[at]->next_atom;
	j = ss->atom_arr[at]->previous_atom;

	if (j == -1)
    {
        // this is the first atom on this list, so make the zone point to
		// the next element in the list.  Note that if the zone had only
		// one element, i should be -1, which will alert the offset that
		// the zone is empty

		findzone(&xzone, &yzone, &zzone, ss->atom_arr[at]->lattice[0], ss->atom_arr[at]->lattice[1], ss->atom_arr[at]->lattice[2], se);
		ss->zone_arr[xzone][yzone][zzone].offset = i;

		if (i != -1)
            ss->atom_arr[i]->previous_atom = -1;
	}
	else
    {
        if (i == -1)
        {
		    // this is the last element on this list,
			ss->atom_arr[j]->next_atom = -1;
		}
		else
        {
			// atom is embedded in the list, nothing special needs be done
			ss->atom_arr[i]->previous_atom = j;
			ss->atom_arr[j]->next_atom = i;
		}
	}

	// now move atom from the end of the atom list to this spot

	if (at != (ss->atom_cnt-1))
	{
		//copy_atom(nat-1, at);
        move_atom((ss->atom_cnt-1), at, ss->atom_arr, ss->zone_arr, ss->transition_arr, se);

	    for (i=0;i<nnt;++i)
			if (nt[i] == (ss->atom_cnt-1))
				nt[i] = at;
	}

	free(ss->atom_arr[ss->atom_cnt-1]);
	--ss->atom_cnt;

	// re-incarnate the buried atoms of the same type

	nnb = 0;

	for (i=0;i<number_of_new_atoms;++i)
    {
        x = new_atom[i].x;
        y = new_atom[i].y;
        z = new_atom[i].z;

        vc = new_atom[i].vc;

        nb[nnb] = reincarnate_atom(x,y,z,type,vc);
        ++nnb;
    }

	nnr = 0;

	for (i=0;i<number_of_new_random_atoms;++i)
    {

        x = new_random_atom[i].x;
        y = new_random_atom[i].y;
        z = new_random_atom[i].z;

        vc = new_random_atom[i].vc;

		subv = drandj(&rand_seed);

        /*if (subv < substrate_percent_a)
			nr[nnr] = random_reincarnate_atom(x,y,z,1,vc);
		else if (subv < (substrate_percent_a + IMPURITY_CONCENTRATION))
            nr[nnr] = random_reincarnate_atom(x,y,z,3,vc);
		else
			nr[nnr]= random_reincarnate_atom(x,y,z,2,vc);
		*/

		if (subv < se->substrate_percent_a)
	        nr[nnr] = random_reincarnate_atom(x,y,z,1,vc);
		else if (subv < (se->substrate_percent_a + se->substrate_percent_b))
            nr[nnr] = random_reincarnate_atom(x,y,z,2,vc);
		else
			nr[nnr]= random_reincarnate_atom(x,y,z,3,vc);

        ++nnr;
        }

	// all atoms have now been incarnated, so refresh transitions of new atoms

	for (i=0;i<nnt;++i)
		refresh_transitions(nt[i], ss, se);

	for (i=0;i<nnb;++i)
		refresh_transitions(nb[i], ss, se);

	for (i=0;i<nnr;++i)
		refresh_transitions(nr[i], ss, se);

	return;
}

/******************************************************************************/
/******************************************************************************/

// checks if there is an atom at point (cx, cy, cz).
// If so, it returns the index to that atom.  If not, return -1.
int atom_at(int cx, int cy, int cz, Atom** atom_arr, Zone zone_arr[ZONES_IN_X][ZONES_IN_Y][ZONES_IN_Z], struct SimulationEnv *se) // lattice coordinate xyz
{
	int i;
	int zx, zy, zz;

	// find proper zone
	findzone(&zx, &zy, &zz, cx, cy, cz, se);

	// cycle through the zone linked list
	i = zone_arr[zx][zy][zz].offset;
	// TODO: if doubles for lattice coordinates are necessary, then should probably change from == to fabs(a-b) < epsilon
	while (i != -1)
	{
		if ((atom_arr[i]->lattice[0] == cx)&&
			(atom_arr[i]->lattice[1] == cy)&&
			(atom_arr[i]->lattice[2] == cz))
			return i;
		else
			i = atom_arr[i]->next_atom;
	}

	return -1;	// no atom
}

int reincarnate(int x, int y, int z, int type, int vc, int buried) {
	// int i,j;

	// int xzone, yzone, zzone;
	// int checkx, checky, checkz;

	// create_default_atom(atom_cnt);

	// // first, find out what zone we're in

	// findzone(&xzone, &yzone, &zzone, x,y,z);

	// if (zone_arr[xzone][yzone][zzone].offset == -1)		// first atom in linked list
	// {
    //     zone_arr[xzone][yzone][zzone].offset = atom_cnt;
	//     atom_arr[atom_cnt]->next_atom = -1;							// no valid link
	// 	atom_arr[atom_cnt]->previous_atom = -1;
	// }
	// else		// link this atom to the others in the zone linked list
	// {

	// 	j = zone_arr[xzone][yzone][zzone].offset;				// first element of list

	// 	while (atom_arr[j]->next_atom != -1)
    //   		j = atom_arr[j]->next_atom;

	// 	// j points to the previous last atom in the zone linked list
	// 	// and points to nothing

	// 	atom_arr[j]->next_atom = atom_cnt;
	// 	atom_arr[atom_cnt]->previous_atom = j;
	// 	atom_arr[atom_cnt]->next_atom = -1;
	// }

	// atom_arr[atom_cnt]->lattice[0] = x;
	// atom_arr[atom_cnt]->lattice[1] = y;
	// atom_arr[atom_cnt]->lattice[2] = z;

	// atom_arr[atom_cnt]->type = type;
	// strcpy(atom_arr[atom_cnt]->name, atom_names[type-1]);

	// /*atom[atom_cnt]->color[0] = atom_color[type].r; //TODO: might not need this part
	// atom[atom_cnt]->color[1] = atom_color[type].g;
	// atom[atom_cnt]->color[2] = atom_color[type].b;*/


	// // find (or set) the occupied neighbor sites
	// // except for the vc direction, all other spots should be occupied

	// for (i=0;i < max_neighbors; ++i)
	// {
	// 	checkx = x + jump_offset[i].dx;
	// 	checky = y + jump_offset[i].dy;
	// 	checkz = z + jump_offset[i].dz;

	// 	adjust_pbc(&checkx, &checky, &checkz, se);

	// 	atom_arr[atom_cnt]->occupied_neighbor_sites[i] = atom_at(checkx, checky, checkz);

	// 	if ((atom_arr[atom_cnt]->occupied_neighbor_sites[i] == -1)&&(i != vc))
	// 		atom_arr[atom_cnt]->occupied_neighbor_sites[i] = buried; //-2 for non-random, -3 for random
	// }

	// // we'll have to add this atom to the transition lists, and also
	// // update the configuration of the atoms contained in n[i]

	// for (i=0; i < max_neighbors; ++i)
	// {
	// 	if (atom_arr[atom_cnt]->occupied_neighbor_sites[i] >= 0)
	// 		atom_arr[atom_arr[atom_cnt]->occupied_neighbor_sites[i]]->occupied_neighbor_sites[opposite_offset[i]] = atom_cnt;

	// 	atom_arr[atom_cnt]->transition_indices[i] = -1;		// initialization
	// }

	// atom_arr[atom_cnt]->transition_indices[max_neighbors] = -1;

	// refresh_transitions(atom_cnt, ss, se);

	// return atom_cnt++; //or do in 2 lines if this doesn't work
	printf("reincarnated something");
	exit(0);
}


/* Maybe want to combine the following 2 methods to just the 1 reincarnate */
int reincarnate_atom(int x, int y, int z, int type, int vc)
{
	return reincarnate(x, y, z, type, vc, -2);
}

int random_reincarnate_atom(int x, int y, int z, int type, int vc)
{
	return reincarnate(x, y, z, type, vc, -3);
}

/******************************************************************************/
/******************************************************************************/

// Remove the given atom from the list of atoms

void kill_atom(int atom_number, struct SimulationState *ss, struct SimulationEnv *se)
{
	int i=0, j;
	int xzone, yzone, zzone;

	//section about cosmetic bonds removed

	// If atom is selected, gracefully remove it from the list of selected atoms: TODO: remove me?

	/*if (atom[atom_number]->selected)
	{
		for (i=0;i<number_selected;++i)
			if (*selected[i] == atom_number)
			{
				*selected[i] = *selected[number_selected-1];
				break;
			}

		free(selected[number_selected]);
		--number_selected;
	}*/

	// copy the last atom in the atom list to this spot, and then free up the memory of the last atom

	i = ss->atom_arr[atom_number]->next_atom;
	j = ss->atom_arr[atom_number]->previous_atom;

	if (j == -1)
    {
        // this is the first atom on this list, so make the zone point to
		// the next element in the list.  Note that if the zone had only
		// one element, i should be -1, which will alert the offset that
		// the zone is empty

		findzone(&xzone, &yzone, &zzone, ss->atom_arr[atom_number]->lattice[0], ss->atom_arr[atom_number]->lattice[1], ss->atom_arr[atom_number]->lattice[2], se);
		ss->zone_arr[xzone][yzone][zzone].offset = i;

		if (i != -1)
            ss->atom_arr[i]->previous_atom = -1;
	}
	else
    {
        if (i == -1)
        {
		    // this is the last element on this list,
			ss->atom_arr[j]->next_atom = -1;
		}
		else
        {
			// atom is embedded in the list, nothing special needs be done

			ss->atom_arr[i]->previous_atom = j;
			ss->atom_arr[j]->next_atom = i;
		}
	}

	// make all atoms whose neighbor this was see an empty spot

	for (i=0; i < se->max_neighbors; ++i)
	{
		if (ss->atom_arr[atom_number]->occupied_neighbor_sites[i] >= 0)
		{
			j = ss->atom_arr[atom_number]->occupied_neighbor_sites[i];
			ss->atom_arr[j]->occupied_neighbor_sites[opposite_offset[i]] = -1;
		}
	}

	// now move atom from the end of the atom list to this spot

	if (atom_number != (ss->atom_cnt-1))
		move_atom((ss->atom_cnt-1), atom_number, ss->atom_arr, ss->zone_arr, ss->transition_arr, se);

	copy_atom(atom_number, ss->atom_cnt-1, ss->atom_arr); //TODO: does this realistically need to happen?

	free(ss->atom_arr[ss->atom_cnt-1]);
	--ss->atom_cnt;

	return;
}

/******************************************************************************/
/******************************************************************************/

// Copy all aspects of atom j into atom i

void copy_atom(int i, int j, Atom** atom_arr)
{
	int m, bn;

	strncpy(atom_arr[i]->name, atom_arr[j]->name, 24); //limited to 24 bc buffer size
	atom_arr[i]->type = atom_arr[j]->type;

	for (m=0;m<3;++m)
      	{
	        atom_arr[i]->cart_coord[m] = atom_arr[j]->cart_coord[m];
			atom_arr[i]->lattice[m] = atom_arr[j]->lattice[m];
		}

	atom_arr[i]->bsradius = atom_arr[j]->bsradius;
	//atom[i]->sfradius = atom[j]->sfradius;

	//atom[i]->visible = atom[j]->visible;
	//atom[i]->selected = atom[j]->selected;

	//atom[i]->style = atom[j]->style;

	/*for (m=0;m<3;++m)
		atom[i]->color[m] = atom[j]->color[m];*/

	for (m=0;m<MAXIMUM_NUMBER_OF_NEIGHBORS+DISSOLUTION;++m)
		atom_arr[i]->transition_indices[m] = atom_arr[j]->transition_indices[m];

	// should there be something like transition_arr[] = new atom number?

	//copy_atom, unlike move_atom, does not preserve the simulation linked list

	atom_arr[i]->next_atom = atom_arr[j]->next_atom;
	atom_arr[i]->previous_atom = atom_arr[j]->previous_atom;

	/*if (generic_flag == 1) return;

	for (m=0;m<atom[j]->nob;++m)
		{
			atom[i]->bond[m] = atom[j]->bond[m];

			bn = atom[i]->bond[m];

			if (bond[bn].from == j) bond[bn].from = i;
			if (bond[bn].to == j) bond[bn].to = i;
		}

	atom[i]->nob = atom[j]->nob;
	atom[i]->biso = atom[j]->biso;*/ //all this seems irrelevant

	return;
}

/******************************************************************************/
/******************************************************************************/
// conversion from aotm->lattice to cartesian and store in atom->cart_coord
void organize(Atom** atom_arr, int atom_cnt) // atom_cnt=number of atoms
{
	// like copy_xyz_to_coord, [but adjusts center of gravity, too, if not commented out]

	orthomol(atom_arr, atom_cnt, primitive_basis);	// use the cell dimensions to orthogonalize
	//centerg(a, atom_cnt); //don't do this!
}

// Orthogonalize all the lattice coordinates according to the cell orthogonalization matrix (com)
// Convert lattice coordinates into cartesian coordinates
void orthomol(Atom** atom_arr, int atom_cnt, double basis[3][3])
{
	int k;

	for (k=0;k<atom_cnt;++k)
    	lattice2cartesian(atom_arr[k]->lattice, basis, atom_arr[k]->cart_coord);

	return;
}

// Translate the coordinates of the atoms so that their center of gravity is on the origin

void centerg(Atom** atom_arr, int atom_cnt)
{
	int i,j;

	for (i=0;i<3;++i)
		centroid[i] = (double)0.0;

	if (atom_cnt == 0)
      	return;

	for(j=0;j<atom_cnt;++j)
		for(i=0;i<3;++i)
			centroid[i] = centroid[i] + atom_arr[j]->cart_coord[i];

	for(i=0;i<3;++i)
		centroid[i] = centroid[i] / (double)atom_cnt;

	for(j=0;j<atom_cnt;++j)
		for(i=0;i<3;++i)
			atom_arr[j]->cart_coord[i] = atom_arr[j]->cart_coord[i]-centroid[i];

	dax -= centroid[0];
	day -= centroid[1];
	daz -= centroid[2];

	return;
}

// converts lattice basis vectors to unit cell parameters
void primitive_basis2ucell_params(double primitive_basis[3][3], double ucell_params[6]) // primitive_basis = basis vectors (rows/first index), ucell_params = unit cell parameters
{
	double rad2deg = 180.0/PI; // radians to degrees conversion factor
	// a b c - magnitude of basis0 basis1 basis2 vectors
	// ENHANCE: if a vector was primitive_basis[0][*], then this could be done with fdot(x,x) and mag(x)
	// TODO: flip indices of primitive_basis
	ucell_params[0] = sqrt(primitive_basis[0][0]*primitive_basis[0][0] + primitive_basis[1][0]*primitive_basis[1][0] + primitive_basis[2][0]*primitive_basis[2][0]);
	ucell_params[1] = sqrt(primitive_basis[0][1]*primitive_basis[0][1] + primitive_basis[1][1]*primitive_basis[1][1] + primitive_basis[2][1]*primitive_basis[2][1]);
	ucell_params[2] = sqrt(primitive_basis[0][2]*primitive_basis[0][2] + primitive_basis[1][2]*primitive_basis[1][2] + primitive_basis[2][2]*primitive_basis[2][2]);
	// gamma - angle between basis0 and basis1 = arccos(fdot(basis0, basis1) / (mag(basis0) * mag(basis1))); from cos(theta) = fdot(a,b) / (mag(a)*mag(b))
	ucell_params[5] = primitive_basis[0][0]*primitive_basis[0][1] + primitive_basis[1][0]*primitive_basis[1][1] + primitive_basis[2][0]*primitive_basis[2][1];
	ucell_params[5] = ucell_params[5]/(ucell_params[0]*ucell_params[1]);
	ucell_params[5] = rad2deg*acos(ucell_params[5]);
	// beta - angle between basis0 and basis2
	ucell_params[4] = primitive_basis[0][0]*primitive_basis[0][2] + primitive_basis[1][0]*primitive_basis[1][2] + primitive_basis[2][0]*primitive_basis[2][2];
	ucell_params[4] = ucell_params[4]/(ucell_params[0]*ucell_params[2]);
	ucell_params[4] = rad2deg*acos(ucell_params[4]);
	// alpha - angle between basis1 and basis2
	ucell_params[3] = primitive_basis[0][2]*primitive_basis[0][1] + primitive_basis[1][2]*primitive_basis[1][1] + primitive_basis[2][2]*primitive_basis[2][1];
	ucell_params[3] = ucell_params[3]/(ucell_params[1]*ucell_params[2]);
	ucell_params[3] = rad2deg*acos(ucell_params[3]);
	// ENHANCE: double arithmetic leads to imprecise values
	return;
}

const double DEFAULT_EPSILON = 1e-3;

// checks if double is within range of an integer - returns 1 for true, 0 for false
int int_check(double fvalue, int ireference, double epsilon){
	return fabs(fvalue - (double) ireference) < epsilon;
}

void lattice2int(double fcoords[3], int coords[3], double epsilon){
	for (int dim_idx = 0; dim_idx < 3; dim_idx++){
		double x = fcoords[dim_idx];
		int comp = round(x);
		int res = int_check(x, comp, epsilon);
		assert(res == 1);
		coords[dim_idx] = comp;
	}
}

// convert a site from cartesian coordinates to lattice coordinates
void cartesian2lattice_site(double ccart[3], double invert_primitive_basis[3][3], int clattice[3]){
	double fclattice[3];
	vecmul(ccart, invert_primitive_basis, fclattice);
	lattice2int(fclattice, clattice, DEFAULT_EPSILON);
}


void cartesian2lattice(double ccart[3], double invert_primitive_basis[3][3], double clattice[3]){
	vecmul(ccart, invert_primitive_basis, clattice);
}

// convert lattice coordinates to cartesian coordinates
void lattice2cartesian(int clattice[3], double primitive_basis[3][3], double ccart[3]){
	for (int dim_idx = 0; dim_idx < 3; dim_idx++){ // vecmul
		ccart[dim_idx] = 
			primitive_basis[dim_idx][0] * (double)clattice[0]
			+ primitive_basis[dim_idx][1] * (double)clattice[1]
			+ primitive_basis[dim_idx][2] * (double)clattice[2];
	}
}

// round floating-point val towards nearest integer in direction of target
int round_towards(double val, int target)
{
	// ENHANCE: do using math.h rounding modes
	if (target >= val)
		return (int) ceil(val);
	else
		return (int) floor(val);
}
