#ifndef ALIGN_H_
#define ALIGN_H_


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>


// =====================================================================

#define INF (2 ^ 31 - 1)	// Infinity
#define S 0					// Substitution
#define I 1 				// Insertion
#define D 2 				// Deletion

extern int MATCH, MISMATCH, HGAP, GAP;

// Cells to make up the dynamic programming table.
typedef struct DP_cell {
	int sub;		// Max score from diagonal adjacent cell
	int ins;		// Max score from left adjacent cell
	int del;		// Max score from up adjacent cell
	int score;
} CELL;

// Structure to store data about a particular alignment.
typedef struct report_q {
	struct report_q *next;	// For linking in a list.
	char *align1;			// Alignment sequence for string 1
	char *align2;			// Alignment sequence for string 2
	char *connect;			// Connection sequence (matches and non-matches)
	int score;				// Optimal score
	int match;				// Number of matches
	int mismatch;			// Number of mismatches
	int gap;				// Number of gaps
	int hgap;				// Number of opening gaps
	int mini;				// Uppermost row bounding the alignment.
	int minj;				// Leftmost column bounding the alignment.
	int maxi;				// Lowermost row bounding the alignment.
	int maxj;				// Rightmost column bounding the alignment.
	double identity;		// Identity percentage
	double coverage;		// Percent coverage of the smaller of the two strings.
} REPORT;


// ============================================================================

int align_loc (char*,int,char*,int*,CELL***);



#endif