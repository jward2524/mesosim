#include "Atoms.h"
#include "ErrorM.h"
#include "Simulation.h"
#include "Utils.h"
#include "Vector.h"
#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// updates atom_arr // [ ]: but doesn't update atom_cnt?
void create_default_atom(long int atom_idx, Atom **atom_arr, struct SimulationEnv *se)
{
    atom_arr[atom_idx] = (Atom *)malloc(sizeof(Atom));

    if (atom_arr[atom_idx] == NULL) {
        // TODO: free mallocs before exiting
        fprintf(stderr, "ERROR! Not enough memory to allocate atom %ld\n", atom_idx);
        fprintf(stderr, "Couldn't allocate memory for atom %ld: %s", atom_idx, strerror(errno));
        clean_and_error(errno);
    }

    atom_arr[atom_idx]->type = 1;

    for (int i = 0; i < 3; ++i) {
        atom_arr[atom_idx]->cartesian[i] = 0.0;
        atom_arr[atom_idx]->lattice[i] = 0.0;
    }

    for (int i = 0; i < se->num_transition_vectors + se->dissolution; ++i)
        atom_arr[atom_idx]->transition_indices[i] = -1;

    for (int i = 0; i < se->num_transition_vectors; ++i)
        atom_arr[atom_idx]->neighbor_atom_idxs[i] = -1;

    // linked list structure, for faster finding of atoms by position

    atom_arr[atom_idx]->next_atom = -1;
    atom_arr[atom_idx]->previous_atom = -1;
    atom_arr[atom_idx]->bsradius = -1;
    atom_arr[atom_idx]->energy = 0;

    return;
}

/*******************************************************************************
*******************************************************************************/
// updates [atom_arr], atom_arr[i]->lattice, [atom_arr[n]->cart_coord], atom_cnt,
// zone_arr[xzone][yzone][zzone].offset, atom_arr[pos]->next_atom/prev_atom,
// atom_arr[pos]->transition_indices; returns index in atom_arr
// energy is updated in refresh_transitions
long int add_atom(int u, int v, int w, unsigned char type, int special, struct SimulationState *ss,
                  struct SimulationEnv *se)
{
    // XXX: special isn't really used
    // lattice coordinates xyz, atom type, special atom conditions (unused)

    long int atom_idx;

    int zone_u, zone_v, zone_w;

    long int pos;

    // [ ]: this is a sanity check? iterating over atom list instead of zone (like atom_at)
    if (atom_at(u, v, w, ss->atom_arr, ss->zone_arr, se) >= 0) {
        int num_overlapping = 0;
        for (int i = 0; i < ss->atom_cnt; ++i) {
            if ((ss->atom_arr[i]->lattice[0] == u) && (ss->atom_arr[i]->lattice[1] == v) &&
                (ss->atom_arr[i]->lattice[2] == w))
                ++num_overlapping;
        }

        fprintf(stderr, "ERROR! Unable to add atom %ld; %d other atoms found at (%d, %d, %d)\n",
                ss->atom_cnt, num_overlapping, u, v, w);
        // ss->simulation_should_kill_itself = true;
        clean_and_error(EXIT_FAILURE);
        return ss->atom_cnt;
    }

    // allocate memory pointed to by the last element of the atom list
    pos = ss->atom_cnt; // position in atom array, presumably
    if (pos > se->max_atoms) {
        fprintf(stderr, "More atoms (%ld) than allocated in atom array (%ld)\n", pos,
                se->max_atoms);
        clean_and_error(EXIT_FAILURE);
    }
    create_default_atom(ss->atom_cnt, ss->atom_arr, se);
    ++ss->atom_cnt;

    if (ss->atom_cnt > se->max_atoms) {
        fprintf(stderr, "Number of atoms (%ld) is exceeding set maximum (%ld)\n", ss->atom_cnt,
                se->max_atoms);
        clean_and_error(errno);
    }

    // TODO: this is already done in atom_at - why repeat it
    findzone(&zone_u, &zone_v, &zone_w, u, v, w, se);

    // xzone, yzone, zzone now have a position open at the end of the zone
    // pos points to this location.  mark the spot and increment the number of atoms
    // in the zone.

    // update the zone.  Increment the number of elements.  If the zone was
    // empty, create a link to the first element in that zone

    if (ss->zone_arr[zone_u][zone_v][zone_w].offset == -1) // first atom in [zone?] linked list
    {
        ss->zone_arr[zone_u][zone_v][zone_w].offset = pos;

        ss->atom_arr[pos]->next_atom = -1; // no valid link
        ss->atom_arr[pos]->previous_atom = -1;
    } else {
        // link this atom to the others in the zone linked list

        // first element of list
        atom_idx = ss->zone_arr[zone_u][zone_v][zone_w].offset;

        while (ss->atom_arr[atom_idx]->next_atom != -1)
            atom_idx = ss->atom_arr[atom_idx]->next_atom;

        // j points to the previous last atom in the zone linked list and points to nothing

        ss->atom_arr[atom_idx]->next_atom = pos;
        ss->atom_arr[pos]->previous_atom = atom_idx;
        ss->atom_arr[pos]->next_atom = -1;
    }

    ss->atom_arr[pos]->lattice[0] = u;
    ss->atom_arr[pos]->lattice[1] = v;
    ss->atom_arr[pos]->lattice[2] = w;

    lattice2cartesian(ss->atom_arr[pos]->lattice, se->primitive_basis,
                      ss->atom_arr[pos]->cartesian);

    ss->atom_arr[pos]->type = type;
    // TODO: use pos instead of atom_cnt-1
    // TODO: use snprintf instead of strcpy

    // update neighbors
    // update neighbor's neighbor (only updates index pointing to current atom)
    // find (or set) the occupied neighbor sites
    // identify neighbors and notify them of presence
    for (int i = 0; i < se->num_transition_vectors; ++i) {
        // mark that this atom cannot yet jump in direction i

        ss->atom_arr[pos]->transition_indices[i] = -1;
        // [ ]: if system size is still 1 (in w direction?) and jump isn't zero, skip it?
        // if ((system_size_z == 1)&&(se->transition_vectors[i].dz != 0))
        // {
        // 	continue;
        // }

        // query_u, query_v, query_w point to the neighboring site
        // so update occupied neighbor site according to atom_at(query_u, query_v, query_w);

        switch (special) {
        case NORMAL:
            // normal bonding considerations
            // case NORMAL_NOGO:
            // //remainder of cases become irrelevant when removing burial
            //  set occupied_neighbor_site[i] to the atom at that site.
            //  If there really is an atom there, cross-link it to our new atom.
            ss->atom_arr[pos]->neighbor_atom_idxs[i] =
                atom_at_offset(u, v, w, i, ss->atom_arr, ss->zone_arr, se);

            // if atom is present at potential jump site, fill position in neighbor_atom_idxs of
            // this atom and the found neighbor atom
            if (ss->atom_arr[pos]->neighbor_atom_idxs[i] >= 0) {
                ss->atom_arr[ss->atom_arr[pos]->neighbor_atom_idxs[i]]
                    ->neighbor_atom_idxs[se->opposite_tvectors[i]] = pos;
            }

            break;

        case SPECIFIED: // only used when reading kmc files?
            printf("I'm in specified\n");
            // ss->atom_arr[pos]->neighbor_atom_idxs[i] = temp_atom.neighbor_atom_idxs[i];
            break;
        default:
            // printf("I made it to the default????\n");
            break;
        }
    }

    if (special != NORMAL) // rates will be refreshed soon
        return pos;

    return pos;
}

// copies atom from initial_idx of atom_arr to final_idx
// atom does not change position within simulation, only position in atom_arr
// used only within the simulation routines, things like bonds, etc. are not copied.
void move_atom(long int initial_idx, long int final_idx, Atom **atom_arr, Zone ***zone_arr,
               Transition **transition_arr, struct SimulationEnv *se)
{
    int zone_u, zone_v, zone_w;

    atom_arr[final_idx]->lattice[0] = atom_arr[initial_idx]->lattice[0];
    atom_arr[final_idx]->lattice[1] = atom_arr[initial_idx]->lattice[1];
    atom_arr[final_idx]->lattice[2] = atom_arr[initial_idx]->lattice[2];

    lattice2cartesian(atom_arr[final_idx]->lattice, se->primitive_basis,
                      atom_arr[final_idx]->cartesian);

    /*atom[fa]->color[0] = atom[ia]->color[0];
    atom[fa]->color[1] = atom[ia]->color[1];
    atom[fa]->color[2] = atom[ia]->color[2];*/

    atom_arr[final_idx]->bsradius = atom_arr[initial_idx]->bsradius;
    atom_arr[final_idx]->energy = atom_arr[initial_idx]->energy;

    atom_arr[final_idx]->type = atom_arr[initial_idx]->type;

    // atom[fa]->biso = atom[ia]->biso;

    atom_arr[final_idx]->next_atom = atom_arr[initial_idx]->next_atom;
    atom_arr[final_idx]->previous_atom = atom_arr[initial_idx]->previous_atom;

    if (atom_arr[final_idx]->next_atom >= 0) {
        // update valid link
        atom_arr[atom_arr[final_idx]->next_atom]->previous_atom = final_idx;
    }

    if (atom_arr[final_idx]->previous_atom >= 0) {
        atom_arr[atom_arr[final_idx]->previous_atom]->next_atom = final_idx;
    } else {
        // fa is first element of a zone
        findzone(&zone_u, &zone_v, &zone_w, atom_arr[final_idx]->lattice[0],
                 atom_arr[final_idx]->lattice[1], atom_arr[final_idx]->lattice[2], se);
        zone_arr[zone_u][zone_v][zone_w].offset = final_idx;
    }

    // copy diffusion transitions
    long n2;
    for (int n = 0; n < se->num_transition_vectors; ++n) {
        atom_arr[final_idx]->neighbor_atom_idxs[n] = atom_arr[initial_idx]->neighbor_atom_idxs[n];

        n2 = atom_arr[initial_idx]->neighbor_atom_idxs[n];
        if (n2 >= 0) {
            atom_arr[n2]->neighbor_atom_idxs[se->opposite_tvectors[n]] = final_idx;
        }

        atom_arr[final_idx]->transition_indices[n] = atom_arr[initial_idx]->transition_indices[n];

        if (atom_arr[final_idx]->transition_indices[n] >= 0) {
            transition_arr[atom_arr[final_idx]->transition_indices[n]]->atom_idx = final_idx;
        }
    }

    // copy evaportation transition
    if (se->is_soluble[atom_arr[initial_idx]->type]) {
        atom_arr[final_idx]->transition_indices[se->num_transition_vectors] =
            atom_arr[initial_idx]->transition_indices[se->num_transition_vectors];

        if (atom_arr[final_idx]->transition_indices[se->num_transition_vectors] >= 0) {
            long transition_idx =
                atom_arr[final_idx]->transition_indices[se->num_transition_vectors];
            transition_arr[transition_idx]->atom_idx = final_idx;
        }
    }

    return;
}

void remove_atom(long int atom_idx, struct SimulationState *ss, struct SimulationEnv *se)
{

    long neighbor_idx;
    int xzone, yzone, zzone;

    int number_of_new_atoms, number_of_new_random_atoms;

    double subv;

    struct {
        int u;
        int v;
        int w;
        int vc;
    } new_atom[MAXIMUM_NUMBER_OF_NEIGHBORS], new_random_atom[MAXIMUM_NUMBER_OF_NEIGHBORS];

    long occupied_neighbors[MAXIMUM_NUMBER_OF_NEIGHBORS];
    int occupied_neighbors_cnt = 0;

    int buried_neighbors[MAXIMUM_NUMBER_OF_NEIGHBORS];
    int buried_neighbors_cnt = 0;

    int reincar_neighbors[MAXIMUM_NUMBER_OF_NEIGHBORS];
    int reincar_neighbors_cnt = 0;

    number_of_new_atoms = 0;
    number_of_new_random_atoms = 0;

    // destroy cross-references to neighboring atoms.  Mark this spot as empty
    // [ ]: isn't this done in refresh_transitions too
    // A: yes, but can't call refresh transitions on an atom without an index

    for (int i = 0; i < se->num_transition_vectors; ++i) {
        // [?]: real use of neighbor_atoms_idxs
        neighbor_idx = ss->atom_arr[atom_idx]->neighbor_atom_idxs[i];

        switch (neighbor_idx) // might be irrelevant if burial removed // [ ]: burried
        {
        // -2 and -1 are not used - reincarnation and incarnation?
        case -2:
            // re-incarnate the buried atom.  this atom will be of type "type"
            // we'll have to add an atom at this point, but we'll do this
            // only after we remove the existence of the current atom

            new_atom[number_of_new_atoms].u =
                ss->atom_arr[atom_idx]->lattice[0] + se->transition_vectors[i].dx;
            new_atom[number_of_new_atoms].v =
                ss->atom_arr[atom_idx]->lattice[1] + se->transition_vectors[i].dy;
            new_atom[number_of_new_atoms].w =
                ss->atom_arr[atom_idx]->lattice[2] + se->transition_vectors[i].dz;

            new_atom[number_of_new_atoms].vc = se->opposite_tvectors[i]; // only allowed direction

            adjust_pbc(&new_atom[number_of_new_atoms].u, &new_atom[number_of_new_atoms].v,
                       &new_atom[number_of_new_atoms].w, se);

            ++number_of_new_atoms;
            break;

        case -3:
            // incarnate a random atom

            new_random_atom[number_of_new_random_atoms].u =
                ss->atom_arr[atom_idx]->lattice[0] + se->transition_vectors[i].dx;
            new_random_atom[number_of_new_random_atoms].v =
                ss->atom_arr[atom_idx]->lattice[1] + se->transition_vectors[i].dy;
            new_random_atom[number_of_new_random_atoms].w =
                ss->atom_arr[atom_idx]->lattice[2] + se->transition_vectors[i].dz;

            // only allowed direction
            new_random_atom[number_of_new_random_atoms].vc = se->opposite_tvectors[i];

            adjust_pbc(&new_random_atom[number_of_new_random_atoms].u,
                       &new_random_atom[number_of_new_random_atoms].v,
                       &new_random_atom[number_of_new_random_atoms].w, se);

            ++number_of_new_random_atoms;
            break;

        case -1:
            if (ss->atom_arr[atom_idx]->transition_indices[i] != -1)
                // an open site to jump to
                remove_from_transition_list(atom_idx, i, ss);
            break;

        default: // make other atom see this spot as empty
            ss->atom_arr[neighbor_idx]->neighbor_atom_idxs[se->opposite_tvectors[i]] = -1;
            occupied_neighbors[occupied_neighbors_cnt] = neighbor_idx;
            ++occupied_neighbors_cnt;
            break;
        }
    }

    // dissolution - if not soluble, won't have dissolution transition
    if (se->is_soluble[ss->atom_arr[atom_idx]->type]) {
        if (ss->atom_arr[atom_idx]->transition_indices[se->num_transition_vectors] != -1)
            remove_from_transition_list(atom_idx, se->num_transition_vectors, ss);
    }

    ss->total_internal_energy -= ss->atom_arr[atom_idx]->energy;

    // now get rid of the atom.  This is almost equivalent to burying it.
    // find out what zone we're in

    // remove the atom from the atom list

    long next_atom_idx = ss->atom_arr[atom_idx]->next_atom;
    neighbor_idx = ss->atom_arr[atom_idx]->previous_atom;

    if (neighbor_idx == -1) {
        // this is the first atom on this list, so make the zone point to
        // the next element in the list.  Note that if the zone had only
        // one element, i should be -1, which will alert the offset that
        // the zone is empty

        findzone(&xzone, &yzone, &zzone, ss->atom_arr[atom_idx]->lattice[0],
                 ss->atom_arr[atom_idx]->lattice[1], ss->atom_arr[atom_idx]->lattice[2], se);
        ss->zone_arr[xzone][yzone][zzone].offset = next_atom_idx;

        if (next_atom_idx != -1)
            ss->atom_arr[next_atom_idx]->previous_atom = -1;
    } else {
        if (next_atom_idx == -1) {
            // this is the last element on this list,
            ss->atom_arr[neighbor_idx]->next_atom = -1;
        } else {
            // atom is embedded in the list, nothing special needs be done
            ss->atom_arr[next_atom_idx]->previous_atom = neighbor_idx;
            ss->atom_arr[neighbor_idx]->next_atom = next_atom_idx;
        }
    }

    // now move atom from the end of the atom list to this spot
    if (atom_idx != (ss->atom_cnt - 1)) {
        // copy_atom(nat-1, at);
        move_atom((ss->atom_cnt - 1), atom_idx, ss->atom_arr, ss->zone_arr, ss->transition_arr, se);

        for (int i = 0; i < occupied_neighbors_cnt; ++i)
            if (occupied_neighbors[i] == (ss->atom_cnt - 1))
                occupied_neighbors[i] = atom_idx;
    }

    free(ss->atom_arr[ss->atom_cnt - 1]);
    ss->atom_arr[ss->atom_cnt - 1] = NULL;
    --ss->atom_cnt;

    // re-incarnate the buried atoms of the same type

    // buried_neighbors_cnt = 0;

    for (int i = 0; i < number_of_new_atoms; ++i) {
        // ++buried_neighbors_cnt;
    }

    // reincar_neighbors_cnt = 0;

    for (int i = 0; i < number_of_new_random_atoms; ++i) {
        subv = drand();

        double bar = 0;
        int type = -1;
        do {
            type++;
            bar += se->substrate_composition[type];
        } while (bar < subv);
        // ++reincar_neighbors_cnt;
    }

    return;
}

/******************************************************************************/
/******************************************************************************/

// checks if there is an atom at point (u, v, w) in lattice coordinates.
// If so, it returns the index to that atom.  If not, return -1.
long atom_at(int u, int v, int w, Atom **atom_arr, Zone ***zone_arr, struct SimulationEnv *se)
{
    // lattice coordinates uvw
    long i;
    int zone_u, zone_v, zone_w;

    adjust_pbc(&u, &v, &w, se);

    // find proper zone
    findzone(&zone_u, &zone_v, &zone_w, u, v, w, se);

    // cycle through the zone linked list
    i = zone_arr[zone_u][zone_v][zone_w].offset;
    while (i != -1) {
        if ((atom_arr[i]->lattice[0] == u) && (atom_arr[i]->lattice[1] == v) &&
            (atom_arr[i]->lattice[2] == w))
            return i;
        else
            i = atom_arr[i]->next_atom;
    }

    return -1; // no atom
}

long atom_at_offset(int u, int v, int w, int offset, Atom **atom_arr, Zone ***zone_arr,
                    struct SimulationEnv *se)
{
    int neighbor_x, neighbor_y, neighbor_z;

    neighbor_x = u + se->transition_vectors[offset].dx;
    neighbor_y = v + se->transition_vectors[offset].dy;
    neighbor_z = w + se->transition_vectors[offset].dz;

    adjust_pbc(&neighbor_x, &neighbor_y, &neighbor_z, se);

    long atom_idx = atom_at(neighbor_x, neighbor_y, neighbor_z, atom_arr, zone_arr, se);
    return atom_idx;
}

// fills initial_config with type of neighbors to atom[at], before jump offset_idx
int get_initial_configuration(long atom_idx, int num_transition_vectors, Atom **atom_arr,
                              int initial_config[]) // atom_idx is position in atom list, offset_idx
                                                    // is index in se->transition_vectors
{
    int offset_idx;
    long neighbor_idx;
    int nn_count = 0; // nearest-neighbors

    for (offset_idx = 0; offset_idx < num_transition_vectors; ++offset_idx) {
        // [?]: real use of neighbor_atoms_idxs
        neighbor_idx = atom_arr[atom_idx]->neighbor_atom_idxs[offset_idx];

        if (neighbor_idx == -1)
            initial_config[offset_idx] = -1; // site is empty
        else {
            // site is occupied by some atom
            ++nn_count; // increment number of near neighbors
            initial_config[offset_idx] = atom_arr[neighbor_idx]->type;
        }
    }

    return nn_count;
}

// fills initial_config with type of neighbors to atom[at], after jump in direction
// se->transition_vectors[offset_idx]
int get_final_configuration(long at, int offset_idx, struct SimulationState *ss,
                            struct SimulationEnv *se,
                            int final_config[]) // offset_idx is position in offset list
{
    long atom_idx;
    int new_x, new_y, new_z;

    int nn_cnt = 0; // nearest-neighbors
    // atom position after jump offset_idx
    new_x = ss->atom_arr[at]->lattice[0] + se->transition_vectors[offset_idx].dx;
    new_y = ss->atom_arr[at]->lattice[1] + se->transition_vectors[offset_idx].dy;
    new_z = ss->atom_arr[at]->lattice[2] + se->transition_vectors[offset_idx].dz;

    adjust_pbc(&new_x, &new_y, &new_z, se);

    for (int i = 0; i < se->num_transition_vectors; ++i) {
        if (i == se->opposite_tvectors[offset_idx]) {
            // if direction is where the jump came from, set as empty
            final_config[i] = -1; // hardcode this? // [ ]: is there a case where it won't be empty?
            continue;
        }
        // location of neighbor
        atom_idx = atom_at_offset(new_x, new_y, new_z, i, ss->atom_arr, ss->zone_arr, se);

        if (atom_idx != -1) {
            // if there is an atom present, 'return' its type
            final_config[i] = ss->atom_arr[at]->type;
            ++nn_cnt;
        } else
            final_config[i] = -1;
    }

    return nn_cnt;
}

int get_coordination(long int atom_idx, struct SimulationState *ss, struct SimulationEnv *se)
{
    int coordination = 0;

    for (int i = 0; i < se->num_transition_vectors; ++i) {
        // [?]: real use of neighbor_atoms_idxs
        // long neighbor_idx = ss->atom_arr[atom_idx]->neighbor_atom_idxs[i];
        long neighbor_idx =
            atom_at_offset(ss->atom_arr[atom_idx]->lattice[0], ss->atom_arr[atom_idx]->lattice[1],
                           ss->atom_arr[atom_idx]->lattice[2], i, ss->atom_arr, ss->zone_arr, se);
        if (neighbor_idx >= 0) {
            ++coordination;
        }
    }

    return coordination;
}
