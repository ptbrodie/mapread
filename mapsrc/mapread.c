// Author: Patrick Brodie

#include "mapread.h"


// ============================================================================
// ============================================================================



// ============================================================================
// Prepare Tree Sequence 
// ============================================================================


int *init_leafarray (int length)
// Allocate an array size of length and initialize it to hold -1s.
{
	int i, *leaves;

	// Allocate array.
	leaves = (int*) malloc (sizeof (int) * (length));
	if (!leaves) {
		perror ("Unable to allocate leaf array");
		exit (1);
	}

	// Initialize array to hold all -1s
	for (i = 0; i < length; ++i) {
		leaves[i] = -1;
	}

	return leaves;
}	


void prepare_tree_DFS (struct node *tree)
// Perform a Depth-First traversal of the given tree, populating
// the leaf array and marking each internal node's leaf list as the
// interval [leafarray[node->array_start], leafarray[node->array_end]).
{
	struct node *curr, *rightchild;
	if (tree) {
		if (tree -> sfxnum != -1) { 	// leaf case
			leafarray[nextindex] = tree -> sfxnum;
			tree -> array_start = nextindex;
			tree -> array_end = nextindex;
			++nextindex;
		} else {						// internal node case

			// Recursively visit all children from left to right.
			curr = tree -> leftchild;
			while (curr) {
				prepare_tree_DFS (curr);
				if (!(curr -> rightsib)) {	// Remember rightmost child.
					rightchild = curr;
				}
				curr = curr -> rightsib;
			}

			// Set internal node's leaf list to be the start of its leftmost 
			// child and the end of its rightmost child.
			tree -> array_start = tree -> leftchild -> array_start;
			tree -> array_end = rightchild -> array_end;
		}
	}
}	


void prepare_tree (struct node *tree)
// Prepare the suffix tree for the read mapping.
{
	// Allocate an array the length of the input genome
	leafarray = init_leafarray (slen + 1);

	// Perform a depth-first traversal of the tree recording the leaf list
	// of each node visited, and marking each leaf in the leafarray.
	prepare_tree_DFS (tree);
}


// End of Prepare Tree Sequence. ++++++++++++++++++++++++++++++++++++++++++++++

// ============================================================================
// Map Reads Sequence 
// ============================================================================



struct node *find_loc_BF (struct node *tree, char *read)
// Find the location of the longest common substring between an input read and the genome
// represented by the given suffix tree.
// TODO: Optimize this section by traversing suffix links to find match locations.
{
	struct node *curr, *parent, *deepest;
	int matches, readi, readlen, maxmatches, i;

	readlen = strlen (read) - LAMBDA + 1;	// No need to continue once strlen < LAMBDA
	curr = tree;
	parent = tree;
	deepest = tree;
	matches = 0;
	maxmatches = 0;

	// Iterate over the read matching its suffices against the tree.
	while (*read && readlen) {
		readi = 0;
		curr = get_branch_by_match (read[readi], tree);
		if (curr) {
			// start matching
			i = curr -> starti;
			while (input_string[i] == read[readi]) {
				++matches;
				if (i + 1 == curr -> endi) {
					parent = curr;
					if ((curr = get_branch_by_match (read[++readi], curr)) == NULL) break;
					i = curr -> starti;
				} else {
					++i; ++readi;
				}
			}
			// Exiting loop means a mismatch was seen.
			if (matches > LAMBDA && matches > maxmatches) {
				maxmatches = matches;
				if (!curr) deepest = parent;
				else deepest = curr;
				printf ("Saw %d matches, assigning deepest = %d\n", maxmatches, deepest -> id);
			}
		}
		++read; --readlen; matches = 0;
	}
	return deepest;
}



struct node *find_loc (struct node *tree, char *read)
// Find the location of the longest common substring between an input read and the genome
// represented by the given suffix tree.
{
	struct node *deepest, *curr, *parent;
	int readi = 0, i, r, e, mismatch = 1;
	int readlen = strlen (read) - LAMBDA + 1;

	parent = tree;
	deepest = tree;
	curr = tree;

	while (read[readi] && readi < readlen) {
		parent = curr;
		curr = get_branch_by_match (read[readi], parent);
		if (curr) {
			r = 0;
			i = curr -> starti;
			while (read[readi] == input_string[i]) {
				if (i + 1 == curr -> endi) {
					++readi; mismatch = 0; break;
				} else {
					++readi; ++i; ++r;
				}
			}
			if (mismatch) {
				readi -= r;
				if (parent -> strdepth > deepest -> strdepth)
					deepest = parent;
				curr = parent -> sfxlink;
			} else { 	
				mismatch = 1;
			}
		} else {
			curr = parent -> sfxlink;
		}
	}
	return deepest;
}


void retrieve_substring (char *substr, int start, int end) 
// Retrieve the substring of the input genome[start: end]
{
	int i;
	char *cp;
	printf ("start = %d, end = %d, len = %d\n", start, end, end-start+1);
	if (start < 0) start = 0;
	if (end > slen) end = slen;
	cp = substr;
	for (i = start; i <= end; ++i) {
		*cp++ = input_string[i];
	}
	*cp = 0;
}


void map_reads (struct node *tree, const char *readfile)
// Map the reads one-by-one onto the genome.
{
	char read[READ_LENGTH], readname[NAME_LENGTH], gslice[READ_LENGTH];
	int i = 1, j, readlen, score, comp, matchalign[2];
	struct node *bfres;
	FILE *fp;
	

	fp = open_file_read (readfile);
	fp = get_next_read (read, readname, fp);
	while (fp && i < 20) {
		readlen = strlen (read);						/// !!!!! HACK ALERT !!!!! Remove i < 20 condition! (inserted for testing speed)
		bfres = find_loc_BF (tree, read);
		//printf ("deepest node id = %d\n", bfres -> id);
		//printf ("leaf start = %d, leaf end = %d\n", bfres -> array_start, bfres -> array_end);
		// Perform an alignment between each location in 
		for (j = bfres -> array_start; j <= bfres -> array_end; ++j) {
			retrieve_substring (gslice, j - readlen, j + readlen);
			//printf ("LEAF = %d\n", j);
			//printf ("read = %s\n", read);
			//printf ("slice = %s\n", gslice);
			score = align_loc (gslice, read, matchalign);
			printf ("local align len = %d, #matches = %d\n", matchalign[0], matchalign[1]);
			comp = bf_align (gslice, read);
			printf ("BF found %d exact matching repeat\n", comp);
		}
		printf ("%d\n", i);
		++i;

		fp = get_next_read (read, readname, fp);
	}
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
	int i;
	struct node *tree;


	// 0. Read in genome and alphabet.
	printf ("Reading files ....\n");
	read_alphabet (&alphabet, alphabetfile);
	read_fasta (&name, &genome, alphabet, genomefile);

	// 1. Build the suffix tree.
	printf ("Building suffix tree ....\n");
	tree = build_tree (genome, alphabet);

	//	2. Prepare the tree and record leaf lists.
	printf ("Preparing suffix tree ....\n");
	nextindex = 0;
	prepare_tree (tree);
	/*for (i = 0; i < slen + 1; ++i) {
		printf ("%d ", leafarray[i]);
	}
	printf ("\n");
	getchar ();*/
	// 	3.  Map Reads onto Genome
	printf ("Mapping reads ....\n");
	map_reads (tree, readfile);


	// Clean up
	free_tree (tree);
	free (leafarray);

}


void print_usage_and_exit ()
// Advise the user of usage and exit.
{
	printf ("USAGE: <map read exe> <FASTA genome> <FASTA reads> <alphabet file>\n");
	exit (1);
}





