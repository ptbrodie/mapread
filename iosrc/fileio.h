// Author: Patrick Brodie

#ifndef FILEIO_H_
#define FILEIO_H_


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>


#define 	READ_LENGTH		512
#define		NAME_LENGTH		256


int MATCH, MISMATCH, HGAP, GAP;


// INTERFACE PROTOTYPES



void read_parms (const char*);
void read_fasta (char**, char**, char*, const char*);
void read_alphabet (char**, const char*);
FILE *open_file_read (const char*);
FILE *get_next_read (char*,char*,FILE*);



#endif
