// Author: Patrick Brodie

// MAIN BLOCK

#include "mapread.h"

int main (int argc, char *argv[])
// Get it!
{

	// TIMER VARIABLES =================
	struct timeval starttime, endtime;
    double elapsed;
    // =================================

	if (argc != 4) {
		print_usage_and_exit ();
	} else {
		read_parms ("INPUTS/parameters.config");
		
		// BEGIN TIMER ==========================================
		gettimeofday(&starttime, NULL);
		
		exec_mapread (argv[1], argv[2], argv[3]);

			// END TIMER ============================================
		gettimeofday(&endtime, NULL);

		// Compute and print elapsed time in milliseconds
    	elapsed = (endtime.tv_sec - starttime.tv_sec) * 1000.0;      // sec to ms
    	elapsed += (endtime.tv_usec - starttime.tv_usec) / 1000.0;   // us to ms

    	printf ("Elapsed time: %lf ms\n", elapsed);
	}
	return 0;
}