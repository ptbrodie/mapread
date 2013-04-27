#ifndef FILEIO_H_
#define FILEIO_H_


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>



int MATCH, MISMATCH, HGAP, GAP;


// INTERFACE PROTOTYPES



void read_parms (const char*);
void read_fasta (char**, char**, char*, const char*);
void read_alphabet (char**, const char*);



#endif
