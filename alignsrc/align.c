// Author: Patrick Brodie


#include "align.h"

// ====================================================================
// Implementation  Smith-Waterman Local Alignment algorithm using the Affine
// Gap Penalty.
// 
// CONVENTION: i denotes ROW
//             j denotes COL
// ====================================================================

#define min(X, Y) ((X) < (Y)? (X) : (Y))
#define max(X, Y) ((X) > (Y)? (X) : (Y))
#define max4(W, X, Y, Z) (max((max(W, X)),(max(Y, Z))))
#define sub(a, b) ((a) == (b)? (MATCH) : (MISMATCH))


void allocate_table (CELL ***table, int cols, int rows)
// Allocate an x-cols x y-rows 2D-array in memory to be used as the alignment table.
{
	int i;
	*table = (CELL**) malloc (sizeof(CELL*) * rows);
	if (!*table) {
		printf ("Not enough memory.");
		exit (1);
	}
	for (i = 0; i < rows; ++i) {
		(*table)[i] = (CELL*) malloc (sizeof(CELL) * cols);
		if (!(*table)[i]) {
			printf ("Not enough memory.");
			exit (1);
		}
	}
}


void free_table (CELL ***table, int cols, int rows)
// Free the memory allocated to an x-cols x y-rows 2D-array.
{
	int i;
	for (i = 0; i < rows; ++i) {
		free ((*table)[i]);
	}
	free (*table);
}



void report_id (double identity) 
// Report the identity percentage of a local alignment.
{
	printf ("Identity Percentage: %.2lf\n", identity);
}


void print_table (CELL **table, int cols, int rows, char *s1, char *s2)
// Print an x-cols x y-rows 2D-array.
{
	int i, j, score;
	printf ("  ");
	for (i = 0; i < cols; ++i) {
		if (i == 0)
			printf (". %4c ", '-');
		else 
			printf (". %4c ", s1[i-1]);
	}
	printf ("\n");
	for (i = 0; i < rows; ++i) {
		if (i == 0)
			printf ("%c ", '-');
		else 
			printf ("%c ", s2[i-1]);
		for (j = 0; j < cols; ++j) {
			score = max4 (table[i][j].sub, table[i][j].ins, table[i][j].del, 0);
			printf ("| %4d ", score);
		}
		printf("|\n");
	}
}


void init_table (CELL ***table, int ilo, int jlo, int cols, int rows, char align_type)
// Initialize the dynamic programming table.
{
	int i, j;
	int penalty;
	if (align_type == 'g') penalty = GAP;
	else if (align_type == 'l') penalty = 0;
	else {
		printf ("Invalid alignment type %c\n", align_type);
		exit (1);
	}
	(*table)[ilo][jlo].sub = 0;
	(*table)[ilo][jlo].ins = 0;
	(*table)[ilo][jlo].del = 0;
	(*table)[ilo][jlo].score = 0;
	for (i = ilo + 1; i < rows; ++i) {
		(*table)[i][jlo].sub = -1 * INF;
		(*table)[i][jlo].ins = -1 * INF;
		if (penalty)		// Local does not init with a penalty
			(*table)[i][jlo].del = HGAP + i * penalty;
		else 
			(*table)[i][jlo].del = 0;
		(*table)[i][jlo].score = 0;
	}
	for (j = jlo + 1; j < cols; ++j) {
		(*table)[ilo][j].sub = -1 * INF;
		if (penalty)		// Local does not init with a penalty
			(*table)[ilo][j].ins = HGAP + j * penalty;
		else
			(*table)[ilo][j].ins = 0;
		(*table)[ilo][j].del = -1 * INF;
		(*table)[ilo][j].score = 0;
	}
	for (i = ilo + 1; i < rows; ++i) {
		(*table)[i][cols].sub = -1 * INF;
		(*table)[i][cols].ins = -1 * INF;
		if (penalty)		// Local does not init with a penalty
			(*table)[i][cols].del = HGAP + i * penalty;
		else 
			(*table)[i][cols].del = 0;
		(*table)[i][cols].score = 0;
	}
	for (j = jlo + 1; j < cols; ++j) {
		(*table)[rows][j].sub = -1 * INF;
		if (penalty)		// Local does not init with a penalty
			(*table)[rows][j].ins = HGAP + j * penalty;
		else
			(*table)[rows][j].ins = 0;
		(*table)[rows][j].del = -1 * INF;
		(*table)[rows][j].score = 0;
	}
}


int calc_t_loc (int *type, CELL **table, int i, int j)
// Calculate the local T value at ith row, jth col and store its type (S,I,or D)
{
	int sb = table[i][j].sub;
	int in = table[i][j].ins;
	int de = table[i][j].del;
	int maxm;
	maxm = max4 (sb, in, de, 0);
	if (maxm == 0) *type = -1;
	else if (maxm == de) *type = D;
	else if (maxm == in) *type = I;
	else if (maxm == sb) *type = S;
	else printf ("Error at i=%d, j=%d\n", i, j);
	return maxm;
}


int calculate_table_loc (CELL ***table, int *maxi, int *maxj, char *s1, char *s2, 
							int ilo, int jlo, int ihi, int jhi)
// Calculate the Dynamic Programming Table for local alignment using affine gap
// penalty.
{
	int i, j, type, maxm = 0, score = 0;
	struct timeval start, end;
	double elapsed;
	CELL *curr, *up, *diag, *left;

	for (i = ilo + 1; i < ihi; ++i) {
		for (j = jlo + 1; j < jhi; ++j) {
            curr = &((*table)[i][j]); up = &((*table)[i-1][j]);
            diag = &((*table)[i-1][j-1]); left = &((*table)[i][j-1]);
			curr -> sub = max (diag -> score + sub (s1[j-1], s2[i-1]), 0);
			curr -> ins = max4 (left -> ins + GAP, 
								left -> sub + HGAP + GAP,
								left -> del + HGAP + GAP, 0);
			curr -> del = max4 (up -> del + GAP, 
								up -> sub + HGAP + GAP,
								up -> ins + HGAP + GAP, 0);
			curr -> score = calc_t_loc (&type, *table, i, j);
			if (curr -> score > maxm) {
				maxm = curr -> score; *maxi = i; *maxj = j;
			}
		}
	}

	return maxm;
}


int traceback_loc (int *match, int *mismatch, int *gap, int *hgap,
					CELL **table, int maxi, int maxj, int *mini, int *minj, 
					int ilo_bnd, int jlo_bnd, char *s1, char *s2)
// Traceback along the optimal local alignment path and store it in align1 and align2.
{
	int score, tmp, type;
	int i = maxi, j = maxj;
	int maxlength = i + j;

	score = calc_t_loc (&type, table, i, j);

	// Traverse the path from the optimal score to where it began, collecting
	// symbols to represent it in a report and counting the penalty occurrences.
	while (type >= 0 && (i > ilo_bnd || j > jlo_bnd)) {
		switch (type) {
			case S:	// substitution
				score = calc_t_loc (&type, table, i-1, j-1);
				if (sub (s1[j-1], s2[i-1]) == MATCH) {
					++(*match);
				} else {
					++(*mismatch);
				}
				--i; --j;
				break;
			case D:	// deletion
				tmp = score;
				score = calc_t_loc (&type, table, i-1, j);
				if ((tmp - HGAP - GAP) == score) ++(*hgap);
				++(*gap);
				--i;
				break;
			case I:	// insertion
				tmp = score;
				score = calc_t_loc (&type, table, i, j-1);
				if ((tmp - HGAP - GAP) == score) ++(*hgap);
				++(*gap);
				--j;
				break;
			default:
				printf ("Error aligning at or near T(%d, %d).\n", i, j);
				exit (1);
		}
	}
	*mini = i; *minj = j;

	return score;
}


int align_loc (char *s1, int s1len, char *s2, int *matchalign, CELL ***table)
// Calculate the optimal local alignment for two strings s1 and s2.
{
	int i, ilo, jlo, ihi, jhi, n, m, opt_score, maxi, maxj, mini, minj;
	int match, mismatch, gap, hgap;
	int alignlen;
	struct timeval start, end;
	double identity, elapsed;
	match = mismatch = gap = hgap = 0;

	// Do not try to align null strings.
	if (s1 && s2) {
		n = s1len, m = strlen (s2);
		ilo = 0; jlo = 0; ihi = m; jhi = n;

		// Calculate the alignment between s1 and s2
		init_table (table, ilo, jlo, jhi + 1, ihi + 1, 'l');
			
		opt_score = calculate_table_loc (table, &maxi, &maxj, s1, s2, ilo, jlo, ihi + 1, jhi + 1);

		traceback_loc (&match, &mismatch, &gap, &hgap, *table, 
						maxi, maxj, &mini, &minj, ilo, jlo, s1, s2);

		alignlen = match + mismatch + gap + hgap;
		matchalign[0] = alignlen;
		matchalign[1] = match;

		return opt_score;
	}
	printf ("Cannot align null string\n");
	exit (1);
}


int bf_align (char *s1, char *s2) {
	int i, j, s1len, s2len;
	int jsave, isave, match, max;
	s1len = strlen (s1);
	s2len = strlen (s2);
	max = 0; match = 0;

	for (i = 0; i < s1len; ++i) {

		for (j = 0; j < s2len; ++j) {
			if (s1[i] == s2[j]) {
				jsave = j; isave = i;
				while (s1[i++] == s2[j++]) match++;
				if (match > max) max = match;
				j = jsave; i = isave; match = 0;
			}
		}
	}
	return max;
}

