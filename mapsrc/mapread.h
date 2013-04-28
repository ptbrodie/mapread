// Author: Patrick Brodie

#ifndef MAPREAD_H_
#define MAPREAD_H_


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include "../sfxsrc/suffix.h"
#include "../alignsrc/align.h"
#include "../iosrc/fileio.h"


// MAX LENGTH OF READ is assumed to be 512 here.  In the future this parameter should be discovered by
// A size file or by counting the number of characters on a line in the read file.

#define LAMBDA				25


// References the next index to insert into during the recursive
// preparation of the tree.
int nextindex;

int *leafarray;







#endif




