/*
* Header file for file_openclose.c
* v1.0 Nov 2019
* Author: J. Janek
*/

#ifndef HEADER_FILE_OC
#define HEADER_FILE_OC
#include <cstdio>

int my_fopen_r(FILE **fr, char *filename, char *mode);
int my_fclose(FILE **f, char *filename);
int my_fopen_w(FILE **fw, char *filename, char *mode);
int my_copy(char *srcname, char *destname);
// #define DEBUG // for debugging only...

#endif
