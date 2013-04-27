

#include "mapread.h"


// ============================================================================
// ============================================================================



// ============================================================================
// Prepare Tree Sequence 
// ============================================================================


int *init_leafarray (int length)
// Allocate an array size of length and initialize it to hold -1s.
{
	int *leafarray, i;

	// Allocate array.
	leafarray = (int*) malloc (sizeof (int) * (length));
	if (!leafarray) {
		perror ("Unable to allocate leaf array");
		exit (1);
	}

	// Initialize array to hold all -1s
	for (i = 0; i < length; ++i) {
		leafarray[i] = -1;
	}

	return leafarray;
}	


void prepare_tree_DFS (struct node *tree, int *leafarray)
// Perform a Depth-First traversal of the given tree, populating
// the leaf array and marking each internal node's leaf list as the
// interval [leafarray[node->array_start], leafarray[node->array_end]).
{
	struct node *curr, *rightchild;
	if (tree) {
		if (tree -> sfxnum != -1) { 	// leaf case
			leafarray[nextindex] = tree -> sfxnum;
			if (tree -> strdepth >= LAMBDA) {
				tree -> array_start = nextindex;
				tree -> array_end = nextindex;
			}
			++nextindex;
		} else {						// internal node case
			// Recursively visit all children from left to right.
			curr = tree -> leftchild;
			while (curr) {
				prepare_tree_DFS (curr, leafarray);
				if (!(curr -> rightsib)) {	// Remember rightmost child.
					rightchild = curr;
				}
				curr = curr -> rightsib;
			}

			// Set internal node's leaf list to be the start of its leftmost 
			// child and the end of its rightmost child.
			if (tree -> strdepth >= LAMBDA) {
				tree -> array_start = tree -> leftchild -> array_start;
				tree -> array_end = rightchild -> array_end;
			}
		}
	}
}	


void prepare_tree (struct node *tree)
// Prepare the suffix tree for the read mapping.
{
	int *leafarray, i;

	// Allocate an array the length of the input genome
	leafarray = init_leafarray (slen + 1);

	// Perform a depth-first traversal of the tree recording the leaf list
	// of each node visited, and marking each leaf in the leafarray.
	prepare_tree_DFS (tree, leafarray);
}


// End of Prepare Tree Sequence. ++++++++++++++++++++++++++++++++++++++++++++++

// ============================================================================
// Map Reads Sequence 
// ============================================================================


void map_reads (const char *readfile)
// Map the reads one-by-one onto the genome.
{
	char read[READ_LENGTH], readname[NAME_LENGTH];

	get_next_read (read, readname, fp);

}


// End of Map Reads Sequence. +++++++++++++++++++++++++++++++++++++++++++++++++


void exec_mapread (const char *genomefile, const char *readfile, const char *alphabetfile)
// Execute the read mapping sequence.
//	1.  Build ST
//	2.	Prepare ST
//	3. 	Map Reads
//	4.	Output
{
	char *alphabet, *genome, *name;
	struct node *tree;

	// 0. Read in genome and alphabet.
	read_alphabet (&alphabet, alphabetfile);
	read_fasta (&name, &genome, alphabet, genomefile);

	// 1. Build the suffix tree.
	tree = build_tree (genome, alphabet);

	//	2. Prepare the tree and record leaf lists.
	nextindex = 0;
	prepare_tree (tree);

	// 	3.  Map Reads onto Genome
	map_reads (readfile);


}


void print_usage_and_exit ()
// Advise the user of usage and exit.
{
	printf ("USAGE: <map read exe> <FASTA genome> <FASTA reads> <alphabet file>\n");
	exit (1);
}





