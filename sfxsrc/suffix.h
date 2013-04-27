// Patrick Brodie
// CS 471
// <A.Kalyanaraman>
// 3/23/13
// PA2 - Suffix Trees / McCreight's Algorithm


#ifndef SUFFIX_H_
#define SUFFIX_H_


// ============================================================================
// suffix.h declares the interface for building a suffix tree, along with the
// structures and global variables used during execution.
// ============================================================================


#include <stdio.h>
#include <stdlib.h>
#include <string.h>


// Insertion Types ===============

#define IA		0 	// Known link, Not root
#define IB		1	// Known link, root
#define IIA		2	// Unknown link, Not root
#define IIB		3	// Unknown link, root

// ===============================


// Suffix tree node structure ====

// Tree uses LEFT CHILD / RIGHT SIBLING structure

struct node {
	int id;			// Unique identifier
	int sfxnum;		// Index number 
	int strdepth;	// Path length in characters from root to this node.
	int starti;		// First index in slice of input_string stored by this node
	int endi;		// Last+1 index in slice of input_string stored by this node
	int array_start;			// First index of leaf array range that marks this node's leaves.
	int array_end;				// Last index of leaf array range that marks this node's leaves.
	struct node *sfxlink;		// Pointer to this node's suffix link
	struct node *leftchild;		// Pointer to first child (sorted alphabetically)
	struct node *rightsib;		// Pointer to next sibling (sorted alphabetically)
	struct node *parent;		// Pointer to parent node.
};

// ================================


// Global Variables ===============

struct node *root;	// Stores the Tree for a given session.
struct node *deepest; 		// Stored for reporting the longest exact matching sequence.
char *input_string;	// Stores one instance of the input string to be referenced by all nodes.
int idCnt, slen;	// Unique ID counter, length of input string
int numleaves, numints; 	// For counting leaves and internal nodes.

// ================================

// Build a suffix tree from the given string over the given alphabet.
struct node *build_tree (char*, char*);
// Free the memory allocated to a suffix tree.
void free_tree (struct node*);
// Print the children of the given node.
void print_children (struct node*);
// Print a depth-first traversal of the given tree.
void print_DFS (struct node*);

#endif
