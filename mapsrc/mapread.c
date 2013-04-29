// Author: Patrick Brodie

#include "mapread.h"


// ============================================================================
// Mapread performs a read-mapping algorithm on a set of reads against an input
// Genome.  
// ============================================================================

double X = 90.0, Y = 80.0;

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
			rightchild = curr;
			while (curr) {
				prepare_tree_DFS (curr);
				rightchild = curr;
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



struct node *find_loc_BF (int len, struct node *tree, char *read, int *maxmatches)
// Find the location of the longest common substring between an input read and the genome
// represented by the given suffix tree.
// NOTE: This is the brute force version of the find_loc algorithm.  Start at root for each
// suffix of the read and match it down the tree.
{
	struct node *curr, *parent, *deepest;
	int matches, readi, readlen, i;

	readlen = len - LAMBDA + 1;	// No need to continue once strlen < LAMBDA
	curr = tree;
	parent = tree;
	deepest = tree;
	matches = 0;
	*maxmatches = 0;

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
			if (matches > LAMBDA && matches > *maxmatches) {
				*maxmatches = matches;
				deepest = parent;
			}
		}
		++read; --readlen; matches = 0;
	}
	return deepest;
}



struct node *find_loc (int len, struct node *tree, char *read, int *maxmatches)
// Find the location of the longest common substring between an input read and the genome
// represented by the given suffix tree.
// NOTE: This is the optimized version of the find_loc algorithm.
{
	struct node *deepest, *curr, *parent;
	int readi = 0, i, r, e, mismatch = 1;
	int readlen = len - LAMBDA + 1;

	*maxmatches = 0;
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
				if (parent -> strdepth + r > *maxmatches) {
					deepest = parent;
					*maxmatches = parent -> strdepth + r;
				}
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


char *retrieve_substring (int *len, int start, int end) 
// Retrieve the substring of the input genome[start: end]
{
	int i;
	char *cp;
	
	if (start < 0) start = 0;
	if (end > slen) end = slen;
	*len = end - start;
	return &input_string[start];
}


void map_reads (struct node *tree, const char *readfile, const char *writefile)
// Map the reads one-by-one onto the genome.
{
	char read[READ_LENGTH], readname[NAME_LENGTH], *gslice;
	int i = 0, j, readlen, score, comp, matchalign[2], slicelen, avg;
	int hitstart, hitend, matches, hits = 0, nohits = 0, numleaves = 0;
	double identity, coverage, maxcoverage = 0.0;
	struct node *deepest;
	CELL **table;
	FILE *fp, *fpout;
	
	allocate_table (&table, READ_LENGTH*2, READ_LENGTH);
	printf ("Redirecting output to %s\n", writefile);

	// Open the read file and the output file.
	fpout = open_file_write (writefile);
	fp = open_file_read (readfile);
	fp = get_next_read (read, readname, fp);

	// For each read, find a viable location in the suffix tree and align it with the genome.
	while (fp /*&& i <= 3000*/) {  /// !!!!! HACK ALERT !!!!! Remove i condition! (inserted for testing speed)
		readlen = strlen (read);						
		deepest = find_loc_BF (readlen, tree, read, &matches);

		// Perform an alignment between each location
		if (matches > LAMBDA) { 
			// Loop over the range of values in the leaf array for alignment locales in the genome.
			numleaves += (deepest -> array_end - deepest -> array_start + 1);
			for (j = deepest -> array_start; j <= deepest -> array_end; ++j) {
				gslice = retrieve_substring (&slicelen, leafarray[j] - readlen, leafarray[j] + readlen);
				// Perform local align between the genome slice and the read.
				score = align_loc (gslice, slicelen, read, matchalign, &table);
				identity = ((double) matchalign[1] / (double) matchalign[0]) * 100.0;
				coverage = ((double) matchalign[0] / (double) readlen) * 100.0;

				// Check if the read was a hit.  If so, record it if it was the best so far.
				if (identity >= X && coverage >= Y) {
					if (coverage > maxcoverage) {
						maxcoverage = coverage;
						hitstart = leafarray[j] - readlen;
						hitend = leafarray[j] + readlen;
					}
				}
			}
		}

		// Output a hit if found.
		if (maxcoverage > 0.0) {
			hits++;
			fprintf (fpout, "%s %d %d\n", readname, hitstart, hitend);
		} else {
			nohits++;
			fprintf (fpout, "%s: No hit found.\n", readname);
		}
		maxcoverage = 0.0; ++i;

		// Get the next read from the read file.
		fp = get_next_read (read, readname, fp);
	}
	free_table (&table, READ_LENGTH*2, READ_LENGTH);

	// Print results of the read mapping.
	printf ("\n***************       RESULTS      ********************\n");
	printf ("Number of reads mapped:  %d\n", i);
	printf ("Genome length:           %d\n", slen);
	printf ("Number of HITS:          %d\n", hits);
	printf ("Number of MISSES:        %d\n", nohits);
	avg = numleaves / i;
	printf ("Average number of alignments per read = %d\n", avg);
	printf ("*******************************************************\n\n");
}


// End of Map Reads Sequence. +++++++++++++++++++++++++++++++++++++++++++++++++

// ============================================================================
// Execution of Algorithm
// ============================================================================

void exec_mapread (const char *genomefile, const char *readfile, const char *alphabetfile)
// Execute the read mapping sequence.
//	1.  Build ST
//	2.	Prepare ST
//	3. 	Map Reads
//	4.	Output
{
	char *alphabet, *genome, *name, writefile[256];
	int i;
	struct node *tree;

	// TIMER VARIABLES =================
	struct timeval startwhole, endwhole, startbuild, endbuild, 
			startprep, endprep, startread, endread;
    double elapsedwhole, elapsedbuild, elapsedprep, elapsedread;
    // =================================


	// 0. Read in genome and alphabet.
	printf ("\n0.  Reading files ....\n");
	read_alphabet (&alphabet, alphabetfile);
	read_fasta (&name, &genome, alphabet, genomefile);
	bzero (writefile, 256);
	strcat (writefile, "MappingResults_");
	strcat (writefile, readfile + 7);
	strcat (writefile, ".txt");

	// BEGIN TIMER WHOLE EXECUTION ==========================================
	gettimeofday(&startwhole, NULL);

		// BEGIN TIMER ST BUILD ==========================================
		gettimeofday(&startbuild, NULL);

		// 1. Build the suffix tree.
		printf ("1.  Building suffix tree ....\n");
		tree = build_tree (genome, alphabet);

		// END TIMER ST BUILD ============================================
		gettimeofday(&endbuild, NULL);

		// Compute and print elapsed time in milliseconds
		elapsedbuild = (endbuild.tv_sec - startbuild.tv_sec) * 1000.0;      // sec to ms
		elapsedbuild += (endbuild.tv_usec - startbuild.tv_usec) / 1000.0;   // us to ms

		printf ("      >Elapsed time (ST Build): %lf ms\n", elapsedbuild);

		// BEGIN TIMER PREPARATION ==========================================
		gettimeofday(&startprep, NULL);

		//	2. Prepare the tree and record leaf lists.
		printf ("2.  Preparing suffix tree ....\n");
		nextindex = 0;
		prepare_tree (tree);

		// END TIMER PREPARATION ============================================
		gettimeofday(&endprep, NULL);

		// Compute and print elapsed time in milliseconds
		elapsedprep = (endprep.tv_sec - startprep.tv_sec) * 1000.0;      // sec to ms
		elapsedprep += (endprep.tv_usec - startprep.tv_usec) / 1000.0;   // us to ms

		printf ("      >Elapsed time (Preparation): %lf ms\n", elapsedprep);

		// BEGIN TIMER READ MAPPING ==========================================
		gettimeofday(&startread, NULL);
		
		// 	3.  Map Reads onto Genome
		printf ("3.  Mapping reads ....\n");
		map_reads (tree, readfile, writefile);

		// END TIMER READ MAPPING ============================================
		gettimeofday(&endread, NULL);

		// Compute and print elapsed time in milliseconds
		elapsedread = (endread.tv_sec - startread.tv_sec) * 1000.0;      // sec to ms
		elapsedread += (endread.tv_usec - startread.tv_usec) / 1000.0;   // us to ms

		printf ("      >Elapsed time (Read Mapping): %lf ms\n", elapsedread);

	// END TIMER WHOLE EXECUTION ============================================
	gettimeofday(&endwhole, NULL);

	// Compute and print elapsed time in milliseconds
	elapsedwhole = (endwhole.tv_sec - startwhole.tv_sec) * 1000.0;      // sec to ms
	elapsedwhole += (endwhole.tv_usec - startwhole.tv_usec) / 1000.0;   // us to ms

	printf ("==== Elapsed time (Total): %lf ms\n", elapsedwhole);


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





