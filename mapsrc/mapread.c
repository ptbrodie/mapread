

#include "mapread.h"


// ============================================================================
// ============================================================================





void prepare_tree ()
// Prepare the suffix tree for the read mapping.
{


}



void exec_mapread (const char *genomefile, 
					const char *readfile, 
					const char *alphabetfile)
// Execute the read mapping sequence.
//	1.  Build ST
//	2.	Prepare ST
//	3. 	Map Reads
//	4.	Output
{
	char *alphabet, *genome, *name;

	// 0. Read in genome and alphabet.
	read_alphabet (&alphabet, alphabetfile);
	read_fasta (&name, &genome, alphabet, genomefile);
	// 1. 
	build_tree (genome, alphabet);
	//	2.
	prepare_tree ();



}


void print_usage_and_exit ()
// Advise the user of usage and exit.
{
	printf ("USAGE: <map read exe> <FASTA reference file> <FASTA read file> <alphabet file>\n");
	exit (1);
}





