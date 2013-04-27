// Patrick Brodie
// CptS 471
// <A.Kalyanaraman>
// 2/12/13
// PA1: Global and Local Alignment Using Affine Gap Penalty.


#include "align.h"

// ====================================================================
// Implementation of both Needleman-Wunsch Global Alignment algorithm
// and the Smith-Waterman Local Alignment algorithm using the Affine
// Gap Penalty.
// 
// Usage:  <exe name> <fasta file> <0:glb; 1:loc> <opt: parm config>
//
// CONVENTION: i denotes ROW
//             j denotes COL
// ====================================================================


REPORT *allocate_report(char *align1, char *align2, char *connect, 
					int score, int match, int mismatch, int gap, int hgap,
					int mini, int minj, int maxi, int maxj, 
					double identity, double coverage)
// Allocate a report for the given alignment stats.
{
	REPORT *report = NULL;
	report = (REPORT*) malloc (sizeof (REPORT));
	if (!report) { printf ("Error mallocing report.\n"); exit (1); }
	report -> align1 = (char *) malloc (strlen (align1)+1);
	if (!report -> align1) { printf ("Error mallocing report -> align1.\n"); exit (1); }
	strcpy (report -> align1, align1);
	report -> align1[strlen (align1)] = 0;
	report -> align2 = (char *) malloc (strlen (align2)+1);
	if (!report -> align2) { printf ("Error mallocing report -> align1.\n"); exit (1); }
	strcpy (report -> align2, align2);
	report -> align2[strlen (align2)] = 0;
	report -> connect = (char *) malloc (strlen (connect)+1);
	if (!report -> connect) { printf ("Error mallocing report -> connect.\n"); exit (1); }
	strcpy (report -> connect, connect);
	report -> connect[strlen (connect)] = 0;
	report -> score = score;
	report -> match = match;
	report -> mismatch = mismatch;
	report -> gap = gap;
	report -> hgap = hgap;
	report -> mini = mini;
	report -> minj = minj;
	report -> maxi = maxj;
	report -> maxj = maxj;
	report -> identity = identity;
	report -> coverage = coverage;
	report -> next = NULL;
	return report;
}


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


void free_reports (REPORT **r)
// free a list of reports
{
	REPORT *tmp;
	while (*r) {
		free ((*r) -> align1);
		free ((*r) -> align2);
		free ((*r) -> connect);
		tmp = *r;
		*r = (*r) -> next;
		free (tmp);
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


int min (int x, int y)
// Return the minimum integer from a set of two integers x and y.
{
	return (x < y)? x : y;
}


int max (int x, int y) 
// Return the maximum integer from a set of two integers x and y.
{
	return (x > y)? x : y;
}


int max3 (int x, int y, int z)
// Return the maximum integer from a set of three integers x, y, z.
{
	return max (x, max (y, z));
}


void reverse (char **str, int len)
// reverse a string in place.
{
	char *i, *j, tmp;
	if (*str && len > 0) {
		i = *str; j = &(*str)[len-1];
		while (i < j) {
			tmp = *i;
			*i = *j;
			*j = tmp;
			++i; --j;
		}
	}
}


void print_reportq (REPORT *r) 
// Print the linked list of reports.
{
    if (r) {
        printf ("->[%d,%lf]", r -> score, r -> coverage);
        print_reportq (r -> next);
    } else {
        printf ("-> NULL\n");
    }
}


void push (REPORT **rhead, REPORT *report)
// Push the given report onto the front of the given report list.
{
    report -> next = *rhead;
    *rhead = report;
}


void insert_report (REPORT **rhead, REPORT *report)
// Insert a report into the report queue sorted by score.
{
    REPORT *curr = *rhead;
    if (!*rhead || report -> score > (*rhead) -> score) {
        push (rhead, report);
    } else {
        while (curr -> next && report -> score <= curr -> next -> score) {
            curr = curr -> next;
        }
        push (&(curr -> next), report);
    }
}


void print_alignment (char *align1, char *align2, char *connect, int mini, int minj)
// print the alignment 60 chars per line, with the number of chars
// used before and after the print.
{
	int count = 0, s1ct = mini, s2ct = minj;
	if (!align1 || !align2 || !connect 
		|| ((strlen (connect) != strlen (align1)) 
		&& (strlen (connect) != strlen (align2)))) {
		printf ("Alignment is null or incorrect.\n");
		return;
	}
	while (*align1 && *align2 && *connect) {
		printf ("s1:%6d  ", s1ct+1);
		while (*align1 && count < 60) {
			printf ("%c", *align1); ++count;
			if (*align1++ != '-') ++s1ct;
		}
		printf ("  %d\n           ", s1ct);
		count = 0; 
		while (*connect && count < 60) {
			printf ("%c", *connect++); ++count;
		}
		printf ("\ns2:%6d  ", s2ct+1);
		count = 0;
		while (*align2 && count < 60) {
			printf ("%c", *align2); ++count;
			if (*align2++ != '-') ++s2ct;
		}
		printf ("  %d\n\n", s2ct); count = 0; 
	}
}


void report (char *align1, char *align2, char *connect, 
				 int opt, int match, int mismatch, int gap, int hgap, int mini, int minj)
// Print a report for the alignment, 60 chars per line.
// Report optimal alignment score, #matches, #mismatches, and #gaps.
{
	print_alignment (align1, align2, connect, mini, minj);
	printf ("Report: \n\n");
	printf ("Optimal alignment: %d\n\n", opt);
	printf ("Number of matches: %d\n", match);
	printf ("Number of mismatches: %d\n", mismatch);
	printf ("Number of gaps: %d\n", gap);
	printf ("Number of entry gaps: %d\n", hgap);
}


void report_id (double identity) 
// Report the identity percentage of a local alignment.
{
	printf ("Identity Percentage: %.2lf\n", identity);
}


void report_cov (double coverage) 
// Report the identity percentage of a local alignment.
{
	printf ("Coverage Percentage: %.2lf\n", coverage);
}


void starline ()
// Print a line of stars to separate print sections.
{
	printf ("\n******************************************************************************\n\n");
}


void print_reports (REPORT *r, int num_reports) 
// Print the queue of reports that contains disjoint local alignments
// starting with the alignment with the highest score.
{
	int i = 0; double coverage = 0.0;
	while (r && i < num_reports && coverage < 90.0) {
		printf ("========== Local Alignment %d\n\n", i + 1);
		report (r -> align1, r-> align2, r -> connect, r -> score, r -> match, r -> mismatch,
				r -> gap, r -> hgap, r -> mini, r -> minj);
		report_id (r -> identity);
		report_cov (r -> coverage);
		printf ("\n\n");
		++i; coverage += r -> coverage;
		r = r -> next;
	}
	starline ();
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
			score = max3 (table[i][j].sub, table[i][j].ins, table[i][j].del);
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
	for (i = ilo + 1; i < rows; ++i) {
		(*table)[i][jlo].sub = -1 * INF;
		(*table)[i][jlo].ins = -1 * INF;
		if (penalty)		// Local does not init with a penalty
			(*table)[i][jlo].del = HGAP + i * penalty;
		else 
			(*table)[i][jlo].del = 0;
	}
	for (j = jlo + 1; j < cols; ++j) {
		(*table)[ilo][j].sub = -1 * INF;
		if (penalty)		// Local does not init with a penalty
			(*table)[ilo][j].ins = HGAP + j * penalty;
		else
			(*table)[ilo][j].ins = 0;
		(*table)[ilo][j].del = -1 * INF;
	}
}


int sub (char ai, char bi)
// Return the substitution score for character ai and bi.
{
	if (ai == bi) 
		return MATCH;
	else
		return MISMATCH;
}


int calc_t (int *type, CELL **table, int i, int j) 
// Calculate the T value at the ith row, jth col, and store its type (S,I,or D)
{
	int sb = table[i][j].sub;
	int in = table[i][j].ins;
	int de = table[i][j].del;
	int maxm = max3 (sb,in,de);
	
	if (maxm == de) *type = D;
	else if (maxm == in) *type = I;
	else if (maxm == sb) *type = S;
	else printf ("Error at i=%d j=%d\n", i, j);
	return maxm;
}


int calc_t_loc (int *type, CELL **table, int i, int j)
// Calculate the local T value at ith row, jth col and store its type (S,I,or D)
{
	int sb = table[i][j].sub;
	int in = table[i][j].ins;
	int de = table[i][j].del;
	int maxm = max (max3 (sb,in,de), 0);
	if (maxm == 0) *type = -1;
	else if (maxm == de) *type = D;
	else if (maxm == in) *type = I;
	else if (maxm == sb) *type = S;
	else printf ("Error at i=%d, j=%d\n", i, j);
	return maxm;
}


int get_bounded_max_ij (int *maxi, int *maxj, int ilo_bnd, int jlo_bnd,
						int ihi_bnd, int jhi_bnd, CELL **table)
// Return the maximum local alignment score and store its i, j coords in maxi, maxj.
{
	int i, j, type, curr, max = 0;
	for (i = ihi_bnd; i > ilo_bnd; --i) {
		for (j = jhi_bnd; j > jlo_bnd; --j) {
			curr = calc_t_loc (&type, table, i, j);
			if (curr > max) {
				max = curr;
				*maxi = i; *maxj = j;
			}
		}
	}
	return max;
}


void calculate_table_glb (CELL ***table, char *s1, char *s2, int cols, int rows)
// Calculate the Dynamic Programming Table for global alignment using affine gap
// penalty.
{
	int i, j, type;
	for (i = 1; i < rows; ++i) {
		for (j = 1; j < cols; ++j) {
			(*table)[i][j].sub = calc_t (&type, *table, i-1, j-1) + sub (s1[j-1], s2[i-1]);
			(*table)[i][j].ins = max3 ((*table)[i][j-1].ins + GAP, 
										(*table)[i][j-1].sub + HGAP + GAP,
										(*table)[i][j-1].del + HGAP + GAP);
			(*table)[i][j].del = max3 ((*table)[i-1][j].del + GAP, 
										(*table)[i-1][j].sub + HGAP + GAP,
										(*table)[i-1][j].ins + HGAP + GAP);
		}
	}
}


void calculate_table_loc (CELL ***table, char *s1, char *s2, int ilo, int jlo, int ihi, int jhi)
// Calculate the Dynamic Programming Table for local alignment using affine gap
// penalty.
{
	int i, j, type;
	for (i = ilo + 1; i < ihi; ++i) {
		for (j = jlo + 1; j < jhi; ++j) {
			(*table)[i][j].sub = max (calc_t (&type, *table, i-1, j-1) + sub (s1[j-1], s2[i-1]), 0);
			(*table)[i][j].ins = max (max3 ((*table)[i][j-1].ins + GAP, 
										(*table)[i][j-1].sub + HGAP + GAP,
										(*table)[i][j-1].del + HGAP + GAP), 0);
			(*table)[i][j].del = max (max3 ((*table)[i-1][j].del + GAP, 
										(*table)[i-1][j].sub + HGAP + GAP,
										(*table)[i-1][j].ins + HGAP + GAP), 0);
		}
	}
}


void traceback_glb (char **align1, char **align2, char **connect, 
					int *match, int *mismatch, int *gap, int *hgap,
					CELL **table, char *s1, char *s2)
// Traceback along the optimal alignment path and store it in align1 and align2.
{
	int score, index = 0, tmp, type;
	int i = strlen (s2), j = strlen (s1);
	int maxlength = i + j;
	*align1 = (char*) malloc (sizeof (char) * maxlength);
	*align2 = (char*) malloc (sizeof (char) * maxlength);
	*connect = (char*) malloc (sizeof (char) * maxlength);

	score = calc_t (&type, table, i, j);
	while (i > 0 || j > 0) {
		if (j == 0) type = D;
		if (i == 0) type = I;
		switch (type) {
			case S:	// substitution
				score = calc_t (&type, table, i-1, j-1);
				if (sub (s1[j-1], s2[i-1]) == MATCH) {
					++(*match);
					(*connect)[index] = '|';
				} else {
					++(*mismatch);
					(*connect)[index] = ' ';
				}
				(*align1)[index] = s1[--j];
				(*align2)[index] = s2[--i];
				break;
			case D:	// deletion
				tmp = score;
				score = calc_t (&type, table, i-1, j);
				if ((tmp - HGAP - GAP) == score) ++(*hgap);
				++(*gap);
				(*align1)[index] = '-';
				(*align2)[index] = s2[--i];
				(*connect)[index] = ' ';
				break;
			case I:	// insertion
				tmp = score;
				score = calc_t (&type, table, i, j-1);
				if ((tmp - HGAP - GAP) == score) ++(*hgap);
				++(*gap);
				(*align1)[index] = s1[--j];
				(*align2)[index] = '-';
				(*connect)[index] = ' ';
				break;
			default:
				printf ("Error aligning at or near T(%d, %d).\n", i, j);
				exit (1);
		}
		
		++index;
	}
	(*align1)[index] = 0;
	(*align2)[index] = 0;
	(*connect)[index] = 0;

	// right the strings since we collected them in reverse order.
	reverse (align1, index);
	reverse (align2, index);
	reverse (connect, index);
}


int traceback_loc (char **align1, char **align2, char **connect, 
					int *match, int *mismatch, int *gap, int *hgap,
					CELL **table, int maxi, int maxj, int *mini, int *minj, 
					int ilo_bnd, int jlo_bnd, char *s1, char *s2)
// Traceback along the optimal local alignment path and store it in align1 and align2.
{
	int score, index = 0, tmp, type;
	int i = maxi, j = maxj;
	int maxlength = i + j;
	*align1 = (char*) malloc (sizeof (char) * maxlength);
	*align2 = (char*) malloc (sizeof (char) * maxlength);
	*connect = (char*) malloc (sizeof (char) * maxlength);

	score = calc_t_loc (&type, table, i, j);

	// Traverse the path from the optimal score to where it began, collecting
	// symbols to represent it in a report and counting the penalty occurrences.
	while (type >= 0 && (i > ilo_bnd || j > jlo_bnd)) {
		switch (type) {
			case S:	// substitution
				score = calc_t_loc (&type, table, i-1, j-1);
				if (sub (s1[j-1], s2[i-1]) == MATCH) {
					++(*match);
					(*connect)[index] = '|';
				} else {
					++(*mismatch);
					(*connect)[index] = ' ';
				}
				(*align1)[index] = s1[--j];
				(*align2)[index] = s2[--i];
				break;
			case D:	// deletion
				tmp = score;
				score = calc_t_loc (&type, table, i-1, j);
				if ((tmp - HGAP - GAP) == score) ++(*hgap);
				++(*gap);
				(*align1)[index] = '-';
				(*align2)[index] = s2[--i];
				(*connect)[index] = ' ';
				break;
			case I:	// insertion
				tmp = score;
				score = calc_t_loc (&type, table, i, j-1);
				if ((tmp - HGAP - GAP) == score) ++(*hgap);
				++(*gap);
				(*align1)[index] = s1[--j];
				(*align2)[index] = '-';
				(*connect)[index] = ' ';
				break;
			default:
				printf ("Error aligning at or near T(%d, %d).\n", i, j);
				exit (1);
		}
		++index;
	}
	(*align1)[index] = 0;
	(*align2)[index] = 0;
	(*connect)[index] = 0;
	*mini = i; *minj = j;

	// right the strings since we collected them in reverse order.
	reverse (align1, index);
	reverse (align2, index);
	reverse (connect, index);
	return score;
}


void get_disj_reports (REPORT **rhead, double coverage, int count, int ilo, int jlo, int ihi, int jhi,
						CELL **table, char *s1, char *s2)
// Collect reports for the first "count" levels of local alignments.
// An alignment is found within the subtable of table bounded by (ilo, jlo) and (ihi, jhi).
{
	int opt_score=0, match=0, mismatch=0, gap=0, hgap=0, maxi=0, maxj=0, mini=0, minj=0; 
	char *align1, *align2, *connect;
	double identity, loc_coverage;
	REPORT *report;

	// Return if no more alignments should be collected.
	if (ilo > ihi || jlo > jhi || coverage >= 90.0 || count <= 0) {
		return;
	} else {
		// Collect all stats for optimal local alignment in subtable (ilo, jlo) -> (ihi, jhi)
		init_table (&table, ilo, jlo, jhi + 1, ihi + 1, 'l');	
		calculate_table_loc (&table, s1, s2, ilo, jlo, ihi + 1, jhi + 1);
		opt_score = get_bounded_max_ij (&maxi, &maxj, ilo, jlo, ihi, jhi, table);
		traceback_loc (&align1, &align2, &connect, &match, &mismatch, &gap, &hgap, table, 
						maxi, maxj, &mini, &minj, ilo, jlo, s1, s2);
		identity = ((double) match / strlen (align1)) * 100;
		loc_coverage = ((double) (match + mismatch) / min (strlen (s1), strlen (s2))) * 100;
		coverage += loc_coverage;
		// create report of the current alignment
		report = allocate_report (align1, align2, connect, opt_score, match, 
									mismatch, gap, hgap, mini, minj, maxi, maxj, 
									identity, loc_coverage);
		
		// insert report into sorted report queue.
		insert_report (rhead, report);
		free (align1); free (align2); free (connect);

		// get report for low and hi alignments.
		get_disj_reports (rhead, coverage, count - 1, ilo, jlo, mini, minj, table, s1, s2);
		get_disj_reports (rhead, coverage, count - 1, maxi, maxj, ihi, jhi, table, s1, s2);
	}
}


int align_loc (char *s1, char *s2)
// Calculate the optimal local alignment for two strings s1 and s2.
{
	CELL **table;
	REPORT *rq = NULL;
	int i, ilo, jlo, ihi, jhi, n, m;

	// Cannot align null strings.
	if (s1 && s2) {
		
		n = strlen (s1), m = strlen (s2);
		ilo = 0; jlo = 0; ihi = m; jhi = n;

		// Allocate an m+1 by n+1 table.
		allocate_table(&table, n + 1, m + 1);

		// Calculate progressively smaller disjoint alignments and collect
		// Reports describing them.
		get_disj_reports (&rq, 0.0, 3, 0, 0, m, n, table, s1, s2);

		// Output the reports collected in get_disj_reports.
		print_reports (rq, 4);

		// Clean up.
		free_table (&table, n + 1, m + 1);
		free_reports (&rq);
		return 0;
	}
	printf ("Cannot align null string\n");
	exit (1);
}


int align_glb (char *s1, char *s2) 
// Calculate the optimal alignment for two strings s1 and s2.
{
	CELL **table;
	char *align1, *align2, *connect;
	int opt_score=0, type=0, match=0, mismatch=0, gap=0, hgap=0, n, m;
	double identity;
	
	if (s1 && s2) {		// cannot align null string.

		n = strlen (s1), m = strlen (s2);

		// Allocate an m+1 by n+1 table.
		allocate_table(&table, n + 1, m + 1);

		// Initialize the table by lining the top row and left column with gap penalties.
		init_table (&table, 0, 0, n + 1, m + 1, 'g');

		// Calculate T(i, j) for all i, j in the table.
		calculate_table_glb (&table, s1, s2, n + 1, m + 1);

		// Calculate the optimal score in the table.
		opt_score = calc_t (&type, table, m, n);

		// Traceback along the path taken to get to the optimal score.
		traceback_glb (&align1, &align2, &connect, &match, &mismatch, &gap, &hgap, table, s1, s2);
		identity = ((double) match / strlen (align1)) * 100;

		// Print a report for the global alignment.
		report (align1, align2, connect, opt_score, match, mismatch, gap, hgap, 0, 0);
		report_id (identity);
		starline ();
		// Clean up.
		free_table (&table, n + 1, m + 1);
		free (align1); free (align2); free (connect);
		return opt_score;
	}
	printf ("Cannot align null string\n");
	exit (1);
}


void read_fasta_a (char **n1, char **n2, char **s1, char **s2, char *filename)
// Read a fasta file containing two sequences into memory to be aligned.
{
	FILE *fp = NULL;
	struct stat st;
	int c; char *curr;

	if (!(fp = fopen (filename, "r"))) {
		printf ("Cannot open file %s\n.", filename);
		exit(1);
	}
	if (stat (filename, &st) < 0) {
		printf ("Cannot stat file %s\n.", filename);
		fclose (fp);
		exit (1);
	}
	
	// Allocate enough room for the input sequences
	*n1 = (char*) malloc (sizeof (char) * 128);
	*n2 = (char*) malloc (sizeof (char) * 128);
	*s1 = (char*) malloc (st.st_size);
	*s2 = (char*) malloc (st.st_size);
	
	if (!*n1 || !*n2 || !*s1 || !*s2) {
		printf ("Malloc failed while reading fasta.\n");
		fclose (fp);
		exit (1);
	}

	// Read first sequence name.
	curr = *n1;
	c = getc (fp);
	while (c != EOF && c != ' ' && c != '\t' && c != '\n') {
		if (c != '>') *curr++ = c;
		c = getc (fp);
	}
	// Rest of first line is don't-care.
	while (c != EOF && c != '\n') { 
		c = getc (fp); 
	}
	// Read first sequence
	curr = *s1;
	while (c != EOF && c != '>') {
		if (c == 'a' || c == 'g' || c == 'c' || c == 't' 
			|| c == 'A' || c == 'G' || c == 'C' || c == 'T') *curr++ = c;
		c = getc (fp);
	}
	// Read second sequence name
	curr = *n2;
	c = getc (fp);
	while (c != EOF && c != ' ' && c != '\t' && c != '\n') {
		if (c != '>') *curr++ = c;
		c = getc (fp);
	}
	// Rest of sequence line is don't-care.
	while (c != EOF && c != '\n') { c = getc (fp); }
	// Read second sequence
	curr = *s2;
	while (c != EOF && c != '>') {
		if (c == 'a' || c == 'g' || c == 'c' || c == 't' 
			|| c == 'A' || c == 'G' || c == 'C' || c == 'T') *curr++ = c;
		c = getc (fp);
	}
	fclose (fp);
}


print_meta (char *n1, char *n2, char *s1, char *s2, char *type)
// Print a the preliminary data for the given alignment.
{
	starline ();
	printf ("Alignment Type: ");
	if (strcmp (type, "0") == 0) printf ("Global Affine Gap\n\n");
	else if (strcmp (type, "1") == 0) printf ("Local Affine Gap\n\n");
	else {
		printf ("Invalid.\n"); exit (1);
	}
	printf ("Scores:   Match: %3d, Mismatch: %3d, H: %3d, G: %3d\n\n", MATCH, MISMATCH, HGAP, GAP);
	printf ("Sequence 1: \"%s\", length = %d characters\n", n1, (int)strlen (s1));
	printf ("Sequence 2: \"%s\", length = %d characters\n\n", n2, (int)strlen (s2));
}









