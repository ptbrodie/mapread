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
	int i;

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


int *prepare_tree (struct node *tree)
// Prepare the suffix tree for the read mapping.
{
	int i;

	// Allocate an array the length of the input genome
	leafarray = init_leafarray (slen + 1);

	// Perform a depth-first traversal of the tree recording the leaf list
	// of each node visited, and marking each leaf in the leafarray.
	prepare_tree_DFS (tree, leafarray);
	return leafarray;
}


// End of Prepare Tree Sequence. ++++++++++++++++++++++++++++++++++++++++++++++

// ============================================================================
// Map Reads Sequence 
// ============================================================================

/**************** Algorithm for FindLoc function (Step 3b) begin ***********************

// note: throughout this procedure you will never modify anything in the suffix tree. In other words, 
the suffix tree will be used as read-only.

0. Initializations: 
        i) struct node *T = the root of the suffix tree.         // "tree pointer"
        ii) int  read_ptr = 1;                                   // "read pointer" (again, use 0 in your code). 
                                                                 // Update this pointer as you match up each 
                                                                 // consecutive character along the read successfully.
                                                                 // Don't increment it if you see a mismatch.

1. Starting at T, find a path below the node T in the tree that spells out as many remaining characters 
of r starting from read_ptr. This would be a simple call to the FindPath() routine starting at the root (T) 
and the starting index for comparison on the read (read_ptr). 

2. Let say, after some number of successful character comparisons, a mismatch is observed and so the matching 
path terminates. There are two subcases here. The matching path could either terminate at the end of an edge 
(case A), or it could terminate be in the middle of an edge (case B). Either way, let u be the internal node 
last visited along the matching path. In case A, u will be the node at which this edge ends. For case B, u will 
be the node from which this edge spawns off.  If case B, then decrease the value of read_ptr by the number of 
characters that successfully matched below u along the edge before you hit a mismatch - i.e., if you matched r 
characters successfully along an edge below u before you hit a mismatch along that edge, then 
read_ptr = read_ptr - r. This will effectively reset the read_ptr back to where it was right after u. 
(Note, for case A, you don't need to do this since the mismatch occurred right after visiting node u.)

3. If the string-depth(u) ≥ l and if the string-depth is the longest seen so far for this read, then store a 
pointer, called "deepestNode" to the node u.  We will update this pointer in future iterations if need be.

4. Now, simply take the suffix link pointer from u to v. Set T=v, and then go back to  step 1, and repeat the process until you 
find the next mismatching point and so on.

5. At some point you will exhaust comparing all characters in your read. That signifies the end of 
iterations. Exit out of the while/for loop (from steps 1-4).

6. Upon exiting the loop, go to the node pointed to by the most up-to-date deepestNode pointer. The suffix 
ids in the interval A[deepestNode->start_index] to A[deepestNode->end_index] is to be returned as the candidate 
list of genomic positions Li for the read.  (Now, there is a possibility that this node's path-label doesn't 
really correspond to the "longest" common substring between the read and genome, but if that happens it will 
only be slightly shy of the length in practice. So ignore this slight approximation in algorithm and use this 
algorithm.)



		// Match along this branch. 
			r = 1; 
			i = curr -> starti;
			e = curr -> endi - curr -> starti;

			while (read[readi] == input_string[i]) {
				//printf ("Match    %c    %c\n", read[readi], input_string[i]);
				if (!(--e)) {
					parent = curr;
					r = 1;
					if ((curr = get_branch_by_match (read[++readi], curr)) == NULL) break;
					i = curr -> starti;
					e = curr -> endi - curr -> starti;
				} else {
					++readi; ++i; ++r;
				}
			}	// Exiting this loop means we have found a mismatch
			//printf ("MISMATCH %c    %c\n", read[readi], input_string[i]);
			readi -= r;
			if (readi > LAMBDA && parent -> strdepth > deepest -> strdepth) deepest = parent;
			curr = parent -> sfxlink;
	


**************** Algorithm for FindLoc function (Step 3b) end ***********************/


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
				deepest = parent;
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


char *retrieve_substring (int start, int end) 
// Retrieve the substring of the input genome[start: end]
{
	int i;
	char *substring, *cp;
	printf ("start = %d, end = %d, len = %d\n", start, end, end-start+1);
	if (start < 0) start = 0;
	if (end > slen) end = slen;
	substring = (char*) malloc (sizeof (char) * (end - start + 1));
	cp = substring;
	for (i = start; i <= end; ++i) {
		*cp++ = input_string[i];
	}
	*cp = 0;
	return substring;

}


void map_reads (struct node *tree, const char *readfile)
// Map the reads one-by-one onto the genome.
{
	char read[READ_LENGTH], readname[NAME_LENGTH], *gslice;
	struct node *bfres;
	FILE *fp;
	int i = 0, j, readlen;

	fp = open_file_read (readfile);
	fp = get_next_read (read, readname, fp);
	while (fp && i < 20) {
		readlen = strlen (read);						/// !!!!! HACK ALERT !!!!! Remove i < 20 condition! (inserted for testing speed)
		bfres = find_loc_BF (tree, read);
		// Perform an alignment between each location in 
		for (j = bfres -> array_start; j <= bfres -> array_end; ++j) {
			gslice = retrieve_substring (j - readlen, j + readlen);
			printf ("read = %s\n", read);
			printf ("slice = %s\n", gslice);
			getchar ();
			align_loc (gslice, read);
		}
		// align

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
	int *leafarray;
	struct node *tree;
	// TIMER VARIABLES =================
	struct timeval starttime, endtime;
    double elapsed;
    // =================================

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
	leafarray = prepare_tree (tree);

	// BEGIN TIMER ==========================================
	gettimeofday(&starttime, NULL);

	// 	3.  Map Reads onto Genome
	printf ("Mapping reads ....\n");
	map_reads (tree, readfile);

	// END TIMER ============================================
	gettimeofday(&endtime, NULL);

	// Compute and print elapsed time in milliseconds
    elapsed = (endtime.tv_sec - starttime.tv_sec) * 1000.0;      // sec to ms
    elapsed += (endtime.tv_usec - starttime.tv_usec) / 1000.0;   // us to ms

    printf ("Elapsed time for map reads: %lf ms\n", elapsed);
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





