// Author: Patrick Brodie

// MAIN BLOCK

int main (int argc, char *argv[])
// Get it!
{
	if (argc != 4) {
		print_usage_and_exit ();
	} else {
		exec_mapread (argv[1], argv[2], argv[3]);
	}
	return 0;
}