#include "stdafx.h"
#include "Defs.h"
#include "Geometry.h"
#include "Vector.h"
#include "Random.h"
#include "Atoms.h"
#include "Simulation_Aux.h"
#include "Simulation.h"

// [ ]: atoms, zones, orientation

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

Atom_Color atom_color[10];

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

void create_default_atom(int n)
{
	int i,j;
	char errorstring[256];

	atom[n] = (Atom *)malloc(sizeof(Atom));

	if (atom[n] == NULL)
	{
		printf("ERROR! Not enough memory to allocate atom %d\n", n);
		return;
	}

	strcpy(atom[n]->name, DEFAULT_ATOM_NAME); 			
	atom[n]->type = 1;
		
	for (i=0;i<3;++i)
		{
			atom[n]->coord[i] = 0.0;
			atom[n]->lattice[i] = 0.0;
		}

	atom[n]->bsradius = DEFAULT_BS_RADIUS; //set or optional
	//atom[n]->sfradius = DEFAULT_SF_RADIUS; //set or optional
		
	//atom[n]->visible = true; //can remove
	//atom[n]->selected = false; //can remove

	//atom[n]->style = DEFAULT_ATOM_STYLE; //can remove

	/*atom[n]->color[0] = DEFAULT_ATOM_COLOR_R; //can remove
	atom[n]->color[1] = DEFAULT_ATOM_COLOR_G; //can remove
	atom[n]->color[2] = DEFAULT_ATOM_COLOR_B; //can remove*/

	for (i=0;i<MAXIMUM_NUMBER_OF_NEIGHBORS + DISSOLUTION;++i)
		atom[n]->position_on_transition_list[i] = -1;

	for (i=0;i<MAXIMUM_NUMBER_OF_NEIGHBORS;++i)
		atom[n]->occupied_neighbor_sites[i] = -1;

	// linked list structure

	atom[n]->next_atom = -1;
	atom[n]->previous_atom = -1;
		
	// bonding
	/*for (i=0;i<MAXIMUM_NUMBER_OF_COSMETIC_BONDS;++i) //can remove this part?
		atom[n]->bond[i] = -1;								// links to the bond drawing list
	atom[n]->nob = 0;*/										// number_of_bonds

	// things specific to x-ray structures

	/*atom[n]->biso = 0.0; //can remove all of this
	for (i=0;i<3;++i)
		for (j=0;j<3;++j)
			atom[n]->ecos[i][j] = 0.0;
	for (i=0;i<3;++i)
		atom[n]->erms[i] = 0.0;*/

	return;
}

/*******************************************************************************
*******************************************************************************/

int add_atom(double x, double y, double z, int type, int special)
{
	/*if (x > 60)
		printf("made it in!\n");*/
	int i, j, k, m, n1;

	int xzone, yzone, zzone;
	double checkx, checky, checkz;

	int pos, ct, n2;

	double sp[3], spo[3];

	if (atom_at(x,y,z) >= 0)
	{
		int num_overlapping = 0;
		for (i=0;i<nat;++i)
		{
			if ((atom[i]->lattice[0] == x)&&(atom[i]->lattice[1] == y)&&(atom[i]->lattice[2] == z))
				++num_overlapping;
		}

		printf("ERROR! Unable to add atom %d; %d other atoms found at (%lf, %lf, %lf)\n", nat, num_overlapping, x, y, z);
		simulation_should_kill_itself = true;
		return nat;
	}

	// allocate memory pointed to by the last element of the atom list
	/*if (x > 60)
		printf("atom of type %d being added at %lf %lf %lf\n", type, x, y, z);*/
	pos = nat;
	create_default_atom(nat);

	/*if (x > 60)
		printf("made it past making a default: nat = %d\n", nat+1);*/
	// set atom color by type
		
	/*atom[nat]->color[0] = atom_color[type].r;
	atom[nat]->color[1] = atom_color[type].g;
	atom[nat]->color[2] = atom_color[type].b;*/

	++nat;

	findzone(&xzone, &yzone, &zzone, x, y, z);

	/*if (x > 60)
		printf("found zone\n");*/
	// xzone, yzone, zzone now have a position open at the end of the zone
	// pos points to this location.  mark the spot and increment the number of atoms
	// in the zone.
	
	// update the zone.  Increment the number of elements.  If the zone was
	// empty, create a link to the first element in that zone

	if (zone[xzone][yzone][zzone].offset == -1)		// first atom in linked list
	{
		zone[xzone][yzone][zzone].offset = pos;

		atom[pos]->next_atom = -1;							// no valid link
		atom[pos]->previous_atom = -1;
	}
	else
	{
		// link this atom to the others in the zone linked list

		j = zone[xzone][yzone][zzone].offset;				// first element of list

		while (atom[j]->next_atom != -1)
			j = atom[j]->next_atom;

		// j points to the previous last atom in the zone linked list and points to nothing

		atom[j]->next_atom = pos;
		atom[pos]->previous_atom = j;
		atom[pos]->next_atom = -1;
	}

	/*if (x > 60)
		printf("other zone logic done\n");*/

	atom[pos]->lattice[0] = x;
	atom[pos]->lattice[1] = y;
	atom[pos]->lattice[2] = z;

	atom[pos]->type = type;
	strcpy(atom[nat-1]->name, atom_names[type-1]);

	/*if (x > 60)
		printf("copied the name: atom is type %s\n", atom[nat-1]->name);*/

	// find (or set) the occupied neighbor sites

	// saturate all the bonds

	for (i=0;i<number_of_possible_neighbors;++i)
	{
		// mark that this atom cannot yet jump in direction i
		/*if (x > 60)
			printf("i am testing neighbor %d\n", i);*/

		atom[pos]->position_on_transition_list[i] = -1;

		if ((ssz == 1)&&(jump_offset[i].dz != 0))
		{
			/*if (x > 60)
				printf("met the corner case\n");*/
			continue;
		}

		// checkx, checky, checkz point to the neighboring site
		// so update occupied neighbor site according to atom_at(checkx, checky, checkz);

		checkx = x + jump_offset[i].dx;
		checky = y + jump_offset[i].dy;
		checkz = z + jump_offset[i].dz;
		/*if (x > 60)
			printf("before pbc xyz %lf %lf %lf\n", checkx, checky, checkz);*/
		adjust_pbc(&checkx, &checky, &checkz);
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
				atom[pos]->occupied_neighbor_sites[i] = atom_at(checkx, checky, checkz);

				/*if (x>60)
					printf("pos = %d, i = %d, atom[pos]->occupied[i] = %d\n", pos, i, atom[pos]->occupied_neighbor_sites[i]);*/

				if (atom[pos]->occupied_neighbor_sites[i] >= 0 ) {
					/*if (x>60)
						printf("reverse: index %d, opposite offset %d, pos %d\n", atom[pos]->occupied_neighbor_sites[i], opposite_offset[i], pos);*/
					atom[atom[pos]->occupied_neighbor_sites[i]]->occupied_neighbor_sites[opposite_offset[i]] = pos;
					/*if (x>60)
						printf("uno reverse didn't break me\n");*/
				}

				/*if (x > 60)
					printf("About to leave the normal case\n");*/
				break;

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
					
				vecmul(sp, latmat, spo);
				unit(spo, spo);

				if (dot(spo, a_pn) > 0.5)
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
					
				vecmul(sp, latmat, spo);
				unit(spo, spo);

				if (dot(spo, a_pn) > 0.5)
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
				//printf("I'm in specified\n");
				atom[pos]->occupied_neighbor_sites[i] = temp_atom.occupied_neighbor_sites[i];
				break;
			default:
				//printf("I made it to the default????\n");
				break;
		}
	}

	//printf("bond saturated\n");
	// can't evaporate either

	atom[pos]->position_on_transition_list[number_of_possible_neighbors] = -1;

	//printf("did you cause a problem\n");
	// set the hopping rates for this atom

	if (special != NORMAL)		// rates will be refreshed soon
		return pos;

	// now set the transition rates

	refresh_transitions(pos);

	//printf("my transition refreshed\n");
	// cycle through the near neighbors, refresh their transitions [or bury as necessary]

	for (i=0;i<number_of_possible_neighbors;++i)
	{
		//printf("trying to refresh neighbor %d\n", i);
		j = atom[pos]->occupied_neighbor_sites[i];

		if (j >= 0)
		{
			k = refresh_transitions(j);					// refresh transitions of neighbor
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

void move_atom(int ia, int fa)
	{
		int n;
		int n2;
		int xzone, yzone, zzone;

		atom[fa]->lattice[0] = atom[ia]->lattice[0];
		atom[fa]->lattice[1] = atom[ia]->lattice[1];
		atom[fa]->lattice[2] = atom[ia]->lattice[2];

		/*atom[fa]->color[0] = atom[ia]->color[0];
		atom[fa]->color[1] = atom[ia]->color[1];
		atom[fa]->color[2] = atom[ia]->color[2];*/

		atom[fa]->bsradius = atom[ia]->bsradius;

		atom[fa]->type = atom[ia]->type;

		//atom[fa]->biso = atom[ia]->biso;

		atom[fa]->next_atom = atom[ia]->next_atom;
		atom[fa]->previous_atom = atom[ia]->previous_atom;

		if (atom[fa]->next_atom >= 0)							// update valid link
      		atom[atom[fa]->next_atom]->previous_atom = fa;

		if (atom[fa]->previous_atom >= 0)
			atom[atom[fa]->previous_atom]->next_atom = fa;
		else													// fa is first element of a zone
			{
				findzone(&xzone, &yzone, &zzone, atom[fa]->lattice[0], atom[fa]->lattice[1], atom[fa]->lattice[2]);
				zone[xzone][yzone][zzone].offset = fa;
			}

		for (n=0;n < number_of_possible_neighbors;++n)
			{
				atom[fa]->occupied_neighbor_sites[n] =	atom[ia]->occupied_neighbor_sites[n];

				n2 = atom[ia]->occupied_neighbor_sites[n];
				if (n2 >= 0) atom[n2]->occupied_neighbor_sites[opposite_offset[n]] = fa;

				atom[fa]->position_on_transition_list[n] = atom[ia]->position_on_transition_list[n];

				if (atom[fa]->position_on_transition_list[n] >= 0)
            	transition_list[atom[fa]->position_on_transition_list[n]]->number_in_list = fa;
			}

		atom[fa]->position_on_transition_list[number_of_possible_neighbors]
      		= atom[ia]->position_on_transition_list[number_of_possible_neighbors];

		if (atom[fa]->position_on_transition_list[number_of_possible_neighbors] >= 0)
      		transition_list[atom[fa]->position_on_transition_list[number_of_possible_neighbors]]->number_in_list = fa;

		return;
	}

/******************************************************************************/
/******************************************************************************/

void remove_atom(int at)
{

	int i,j;
	int type;
	int xzone, yzone, zzone;

	int number_of_new_atoms, number_of_new_random_atoms;
	double x, y, z;
	int vc;

	double subv;

	struct {
			double x;
			double y;
			double z;
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

  	type = atom[at]->type;

	for (i=0;i<number_of_possible_neighbors;++i)
	{
		j = atom[at]->occupied_neighbor_sites[i];

		switch(j) //might be irrelevant if burial removed
		{
			case -2:
				// re-incarnate the buried atom.  this atom will be of type "type"
				// we'll have to add an atom at this point, but we'll do this
				// only after we remove the existence of the current atom
							

				new_atom[number_of_new_atoms].x = atom[at]->lattice[0] + jump_offset[i].dx;
				new_atom[number_of_new_atoms].y = atom[at]->lattice[1] + jump_offset[i].dy;
				new_atom[number_of_new_atoms].z = atom[at]->lattice[2] + jump_offset[i].dz;

				new_atom[number_of_new_atoms].vc = opposite_offset[i];	// only allowed direction

				adjust_pbc(&new_atom[number_of_new_atoms].x,
							&new_atom[number_of_new_atoms].y,
							&new_atom[number_of_new_atoms].z);

				++number_of_new_atoms;
				break;

			case -3:
				// incarnate a random atom

				new_random_atom[number_of_new_random_atoms].x = atom[at]->lattice[0] + jump_offset[i].dx;
				new_random_atom[number_of_new_random_atoms].y = atom[at]->lattice[1] + jump_offset[i].dy;
				new_random_atom[number_of_new_random_atoms].z = atom[at]->lattice[2] + jump_offset[i].dz;

				new_random_atom[number_of_new_random_atoms].vc = opposite_offset[i];	// only allowed direction

				adjust_pbc(&new_random_atom[number_of_new_random_atoms].x,
                     		&new_random_atom[number_of_new_random_atoms].y,
							&new_random_atom[number_of_new_random_atoms].z);

				++number_of_new_random_atoms;
				break;

			case -1:
				if (atom[at]->position_on_transition_list[i] != -1)		// i.e., this is a spot to jump to 
					take_off_transition_list(at, i);
				break;

			default:	// make other atom see this spot as empty
				atom[j]->occupied_neighbor_sites[opposite_offset[i]] = -1;
				nt[nnt] = j;
				++nnt;
				break;
		}
	}

	take_off_transition_list(at, number_of_possible_neighbors);			// dissolution

	// now get rid of the atom.  This is almost equivalent to burying it.
	// find out what zone we're in

	// remove the atom from the atom list

	i = atom[at]->next_atom;
	j = atom[at]->previous_atom;

	if (j == -1)
    {
        // this is the first atom on this list, so make the zone point to
		// the next element in the list.  Note that if the zone had only
		// one element, i should be -1, which will alert the offset that
		// the zone is empty

		findzone(&xzone, &yzone, &zzone, atom[at]->lattice[0], atom[at]->lattice[1], atom[at]->lattice[2]);
		zone[xzone][yzone][zzone].offset = i;

		if (i != -1)
            atom[i]->previous_atom = -1;
	}
	else
    {
        if (i == -1)
        {
		    // this is the last element on this list,
			atom[j]->next_atom = -1;
		}
		else
        {
			// atom is embedded in the list, nothing special needs be done
			atom[i]->previous_atom = j;
			atom[j]->next_atom = i;
		}
	}

	// now move atom from the end of the atom list to this spot

	if (at != (nat-1))
	{
		//copy_atom(nat-1, at);
        move_atom((nat-1), at);

	    for (i=0;i<nnt;++i)
			if (nt[i] == (nat-1))
				nt[i] = at;
	}

	free(atom[nat-1]);
	--nat;

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

		subv = drandj(&seed);

        /*if (subv < substrate_percent_a)
			nr[nnr] = random_reincarnate_atom(x,y,z,1,vc);
		else if (subv < (substrate_percent_a + IMPURITY_CONCENTRATION))
            nr[nnr] = random_reincarnate_atom(x,y,z,3,vc);
		else
			nr[nnr]= random_reincarnate_atom(x,y,z,2,vc);
		*/

		if (subv < substrate_percent_a)
	        nr[nnr] = random_reincarnate_atom(x,y,z,1,vc);
		else if (subv < (substrate_percent_a + substrate_percent_b))
            nr[nnr] = random_reincarnate_atom(x,y,z,2,vc);
		else
			nr[nnr]= random_reincarnate_atom(x,y,z,3,vc);

        ++nnr;
        }

	// all atoms have now been incarnated, so refresh transitions of new atoms

	for (i=0;i<nnt;++i)
		refresh_transitions(nt[i]);

	for (i=0;i<nnb;++i)
		refresh_transitions(nb[i]);

	for (i=0;i<nnr;++i)
		refresh_transitions(nr[i]);

	return;
}

/******************************************************************************/
/******************************************************************************/

// checks if there is an atom at point (cx, cy, cz).
// If so, it returns the index to that atom.  If not, return -1.

int atom_at(double cx, double cy, double cz)
{
	int i;
	int zx, zy, zz;

	// find proper zone
	findzone(&zx, &zy, &zz, cx, cy, cz);

	// cycle through the zone linked list
	i = zone[zx][zy][zz].offset;

	while (i != -1)
	{
		if ((atom[i]->lattice[0] == cx)&&
			(atom[i]->lattice[1] == cy)&&
			(atom[i]->lattice[2] == cz))
			return i;
		else i = atom[i]->next_atom;
	}

	return -1;			// no atom
}

int reincarnate(double x, double y, double z, int type, int vc, int buried) {
	int i,j;

	int xzone, yzone, zzone;
	double checkx, checky, checkz;

	create_default_atom(nat);

	// first, find out what zone we're in

	findzone(&xzone, &yzone, &zzone, x,y,z);

	if (zone[xzone][yzone][zzone].offset == -1)		// first atom in linked list
	{
        zone[xzone][yzone][zzone].offset = nat;
	    atom[nat]->next_atom = -1;							// no valid link
		atom[nat]->previous_atom = -1;
	}
	else		// link this atom to the others in the zone linked list
	{

		j = zone[xzone][yzone][zzone].offset;				// first element of list

		while (atom[j]->next_atom != -1)
      		j = atom[j]->next_atom;

		// j points to the previous last atom in the zone linked list
		// and points to nothing

		atom[j]->next_atom = nat;
		atom[nat]->previous_atom = j;
		atom[nat]->next_atom = -1;
	}

	atom[nat]->lattice[0] = x;
	atom[nat]->lattice[1] = y;
	atom[nat]->lattice[2] = z;

	atom[nat]->type = type;
	strcpy(atom[nat]->name, atom_names[type-1]);

	/*atom[nat]->color[0] = atom_color[type].r; //TODO: might not need this part
	atom[nat]->color[1] = atom_color[type].g;
	atom[nat]->color[2] = atom_color[type].b;*/


	// find (or set) the occupied neighbor sites
	// except for the vc direction, all other spots should be occupied

	for (i=0;i<number_of_possible_neighbors;++i)
	{
		checkx = x + jump_offset[i].dx;
		checky = y + jump_offset[i].dy;
		checkz = z + jump_offset[i].dz;

		adjust_pbc(&checkx, &checky, &checkz);

		atom[nat]->occupied_neighbor_sites[i] = atom_at(checkx, checky, checkz);

		if ((atom[nat]->occupied_neighbor_sites[i] == -1)&&(i != vc))
			atom[nat]->occupied_neighbor_sites[i] = buried; //-2 for non-random, -3 for random
	}

	// we'll have to add this atom to the transition lists, and also
	// update the configuration of the atoms contained in n[i]

	for (i=0;i<number_of_possible_neighbors;++i)
	{
		if (atom[nat]->occupied_neighbor_sites[i] >= 0)
			atom[atom[nat]->occupied_neighbor_sites[i]]->occupied_neighbor_sites[opposite_offset[i]] = nat;

		atom[nat]->position_on_transition_list[i] = -1;		// initialization
	}

	atom[nat]->position_on_transition_list[number_of_possible_neighbors] = -1;

	refresh_transitions(nat);

	return nat++; //or do in 2 lines if this doesn't work
}


/* Maybe want to combine the following 2 methods to just the 1 reincarnate */
int reincarnate_atom(double x, double y, double z, int type, int vc)
{
	return reincarnate(x, y, z, type, vc, -2);
}

int random_reincarnate_atom(double x, double y, double z, int type, int vc)
{
	return reincarnate(x, y, z, type, vc, -3);
}

/******************************************************************************/
/******************************************************************************/

// Remove the given atom from the list of atoms

void kill_atom(int atom_number)
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

	i = atom[atom_number]->next_atom;
	j = atom[atom_number]->previous_atom;

	if (j == -1)
    {
        // this is the first atom on this list, so make the zone point to
		// the next element in the list.  Note that if the zone had only
		// one element, i should be -1, which will alert the offset that
		// the zone is empty

		findzone(&xzone, &yzone, &zzone, atom[atom_number]->lattice[0], atom[atom_number]->lattice[1], atom[atom_number]->lattice[2]);
		zone[xzone][yzone][zzone].offset = i;

		if (i != -1)
            atom[i]->previous_atom = -1;
	}
	else
    {
        if (i == -1)
        {
		    // this is the last element on this list,
			atom[j]->next_atom = -1;
		}
		else
        {
			// atom is embedded in the list, nothing special needs be done

			atom[i]->previous_atom = j;
			atom[j]->next_atom = i;
		}
	}

	// make all atoms whose neighbor this was see an empty spot

	for (i=0;i<number_of_possible_neighbors;++i)
	{
		if (atom[atom_number]->occupied_neighbor_sites[i] >= 0)
		{
			j = atom[atom_number]->occupied_neighbor_sites[i];
			atom[j]->occupied_neighbor_sites[opposite_offset[i]] = -1;
		}
	}

	// now move atom from the end of the atom list to this spot

	if (atom_number != (nat-1))
		move_atom((nat-1), atom_number);

	copy_atom(atom_number, nat-1); //TODO: does this realistically need to happen?

	free(atom[nat-1]);
	--nat;

	return;
}

/******************************************************************************/
/******************************************************************************/

// Copy all aspects of atom j into atom i

void copy_atom(int i, int j)
{
	int m, bn;

	strncpy(atom[i]->name, atom[j]->name, 24); //limited to 24 bc buffer size
	atom[i]->type = atom[j]->type;

	for (m=0;m<3;++m)
      	{
	        atom[i]->coord[m] = atom[j]->coord[m];
			atom[i]->lattice[m] = atom[j]->lattice[m];
		}

	atom[i]->bsradius = atom[j]->bsradius;
	//atom[i]->sfradius = atom[j]->sfradius;

	//atom[i]->visible = atom[j]->visible;
	//atom[i]->selected = atom[j]->selected;

	//atom[i]->style = atom[j]->style;

	/*for (m=0;m<3;++m)
		atom[i]->color[m] = atom[j]->color[m];*/

	for (m=0;m<MAXIMUM_NUMBER_OF_NEIGHBORS+DISSOLUTION;++m)
		atom[i]->position_on_transition_list[m] = atom[j]->position_on_transition_list[m];

	// should there be something like transition_list[] = new atom number?

	//copy_atom, unlike move_atom, does not preserve the simulation linked list

	atom[i]->next_atom = atom[j]->next_atom;
	atom[i]->previous_atom = atom[j]->previous_atom;

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

void organize(Atom* a[], int na)
{
	// like copy_xyz_to_coord, but adjusts center of gravity, too

	orthomol(a, na, latmat);					// use the cell dimensions to orthogonalize
	//centerg(a, na); //don't do this!
}

// Orthogonalize all the coordinates according to the cell orthogonalization matrix
void orthomol(Atom* atm[], int na, double com[3][3])
{
	int k;

	for (k=0;k<na;++k)
    	vecmul(atm[k]->lattice, com, atm[k]->coord);

	return;
}

// Translate the coordinates of the atoms so that their center of gravity is on the origin

void centerg(Atom* atm[], int na)
{
	int i,j;

	for (i=0;i<3;++i)
		cg[i] = (double)0.0;

	if (na == 0)
      	return;

	for(j=0;j<na;++j)
		for(i=0;i<3;++i)
			cg[i]=cg[i]+atm[j]->coord[i];

	for(i=0;i<3;++i)
		cg[i]=cg[i]/(double)na;

	for(j=0;j<na;++j)
		for(i=0;i<3;++i)
			atm[j]->coord[i]=atm[j]->coord[i]-cg[i];

	dax -= cg[0];
	day -= cg[1];
	daz -= cg[2];

	return;
}

// [ ]: what does this do?
void latmat_to_cell(double com[3][3], double celld[6])
{
	double pir = 180.0/PI;

	celld[0] = sqrt(com[0][0]*com[0][0] + com[1][0]*com[1][0] + com[2][0]*com[2][0]);
	celld[1] = sqrt(com[0][1]*com[0][1] + com[1][1]*com[1][1] + com[2][1]*com[2][1]);
	celld[2] = sqrt(com[0][2]*com[0][2] + com[1][2]*com[1][2] + com[2][2]*com[2][2]);
		
	celld[5] = com[0][0]*com[0][1] + com[1][0]*com[1][1] + com[2][0]*com[2][1];
	celld[5] = celld[5]/(celld[0]*celld[1]);
	celld[5] = pir*acos(celld[5]);

	celld[4] = com[0][0]*com[0][2] + com[1][0]*com[1][2] + com[2][0]*com[2][2];
	celld[4] = celld[4]/(celld[0]*celld[2]);
	celld[4] = pir*acos(celld[4]);

	celld[3] = com[0][2]*com[0][1] + com[1][2]*com[1][1] + com[2][2]*com[2][1];
	celld[3] = celld[3]/(celld[1]*celld[2]);
	celld[3] = pir*acos(celld[3]);

	return;
}

