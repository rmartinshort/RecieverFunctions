/*
   *******************************************************************
   * (c) Copyright 2002 The Regents of the University of California. *
   *     All Rights Reserved.                                        *
   *                                                                 *
   *     This work was produced under the sponsorship of the         *
   *     United States Department of Energy.                         *
   *     The Government retains certain rights therein.              *
   *******************************************************************

 author: Peter Goldstein
 date:   august 2002
 where:  LLNL
 Updates: Brian Savage <savage13@dtm.ciw.edu>
          - Added Error checking for reads and file length
	  - Remove Hard coded numbers for readability

 *********************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

#define PROGRAM "sacswap"

#define NFLOATS    70
#define NINTS      40
#define NCHARS     192

#define FLOAT_SIZE 4
#define INT_SIZE   4
#define CHAR_SIZE  1
#define DATA_SIZE  4


#define LONG_NPTS  9 /* Position of number of points in longs */
#define LONG_NVHDR 6 /* Position of Header Value in longs */
#define HEADER_MIN 1 /* Minimum Header Value */
#define HEADER_MAX 6 /* Maximum Header Value */

#define MAXBYTES 280 /* 70*4 for floats */

#define HEADER_SIZE (NFLOATS * FLOAT_SIZE) + (NINTS * INT_SIZE) + (NCHARS * CHAR_SIZE)

int file_size(char *file);
float float_swap(char cbuf[]);
long  long_swap(char cbuf[]);
long  long_join(char cbuf[]);


int 
main(int ac, char **av) {
  
  int i, j, npts, nvhdr, nnvhdr;
  float *data;
  float sacfloat[NFLOATS];
  long saclong[NINTS];
  char cbuf[MAXBYTES], sacfile[2048], sacchars[NCHARS];
  FILE *fp;
  char swap[2048];
  int size;

  /* Check usage */
  if (ac < 2) { 
    fprintf(stderr, "Usage: %s sacfile(s)\n", PROGRAM);
    fprintf(stderr, "\tconvert sacfiles between byte orders\n");
    fprintf(stderr, "\toutput files are appended with .swap\n");
    fprintf(stderr, "\twill read either byte order\n");
    exit(-1);
  }
  
  /* Loop over files */
  for(i = 1; i < ac; i++) {

    /* Determine file existance and size */
    if((size = file_size(av[i])) == -1) {
      continue;
    }

    /* Open file or skip to next file */
    if ( (fp = fopen(av[i], "rb")) == NULL) {
      fprintf(stderr, "%s Error Opening file %s\n", PROGRAM, av[i]);
      continue;
    }
    
    /* Read the Header */
    /* floats */
    if(fread(cbuf, FLOAT_SIZE * NFLOATS, 1, fp) != 1) {
      fprintf(stderr, "%s: bad read in header: %s\n", PROGRAM, av[i]);
      continue;
    }
    /* Swap the floats */
    for(j = 0; j< NFLOATS; j++) {
      sacfloat[j] = float_swap(cbuf + FLOAT_SIZE * j);
    }

    /* ints */
    if(fread(cbuf, INT_SIZE * NINTS, 1, fp) != 1) {
      fprintf(stderr, "%s: bad read in header: %s\n", PROGRAM, av[i]);
      continue;
    }
    /* Swap the ints */
    for(j = 0; j < NINTS; j++) {
      saclong[j] = long_swap(cbuf + INT_SIZE * j);
    }
    
    /* Figure out which way we're going: To or from native format? */
    /* store swapped and unswapped header values */
    nvhdr  = (int) saclong[LONG_NVHDR];
    nnvhdr = (int) long_join(cbuf + INT_SIZE * LONG_NVHDR);

    if ( nvhdr >= HEADER_MIN && nvhdr <= HEADER_MAX ) {
      npts = (int) saclong[LONG_NPTS];
      sprintf(swap, "non-native => native");
    } else if ( nnvhdr >= HEADER_MIN && nnvhdr <= HEADER_MAX ) {
      npts = (int) long_join(cbuf + INT_SIZE * LONG_NPTS);
      sprintf(swap, "native => non-native");
    } else {
      fprintf(stderr, "%s: %s has bad header (native,nonnative): %d, %d\n", 
	      PROGRAM, av[i], nvhdr, nnvhdr);
      continue;
    }
    /* Check the length of the file */
    if(HEADER_SIZE + DATA_SIZE * npts != size) {
      fprintf(stderr, "%s: Number of points in file is incorrect: %s\n", 
	      PROGRAM, av[i]);
      continue;
    }

    /* chars */
    if(fread(cbuf, NCHARS * CHAR_SIZE, 1, fp) != 1) {
      fprintf(stderr, "%s: bad read in header: %s\n", PROGRAM, av[i]);
      continue;
    }
    for(j = 0; j < NCHARS; j++) {
      sacchars[j] = cbuf[j];
    }
    
    /* Read the Data */
    data = (float *) malloc(npts * DATA_SIZE);
    for (j = 0; j < npts; j++) {
      if(fread(cbuf, DATA_SIZE, 1, fp) != 1) {
	fprintf(stderr, "%s: bad read in data: %s\n", PROGRAM, av[i]);
	free(data);
	fprintf(stderr, "%s: exiting\n", PROGRAM);
	exit(-1);
      }
      data[j] = (float) float_swap(cbuf);
    }
    fclose(fp);
    
    sprintf(sacfile, "%s.swap", av[i]); 
    if ( (fp = fopen(sacfile, "w")) == NULL ) {
      fprintf(stdout, "%s: Error in opening output file: %s\n", 
	      PROGRAM, sacfile);
    }

    fprintf(stderr, "%s: Writing %s.swap with npts = %d %s\n", 
	    PROGRAM, av[i], npts, swap);
    
    /* Write the header */
    fwrite(&sacfloat, NFLOATS * FLOAT_SIZE, 1, fp);
    fwrite(&saclong,  NINTS   * INT_SIZE,   1, fp);
    fwrite(&sacchars, NCHARS  * CHAR_SIZE,  1, fp);
    
    /* Write the data */
    fwrite(&data[0], npts * DATA_SIZE, 1, fp);
    
    fclose(fp);
    free(data);
    
  }
  /* end loop over files*/
  return 0;
}

int 
file_size(char *file) {
  struct stat s;
  char *c;
  if(stat(file, &s) != 0) {
    c = (char *) malloc(sizeof(PROGRAM) + sizeof(file) + 3);
    sprintf(c, "%s: %s", PROGRAM, file);
    perror(c);
    free(c);
    return(-1);
  }
  return((int)s.st_size);
}

long 
long_swap(char cbuf[]) {
  union {
    char cval[4];
    long lval;
  } l_union;
  
  l_union.cval[3] = cbuf[0];
  l_union.cval[2] = cbuf[1];
  l_union.cval[1] = cbuf[2];
  l_union.cval[0] = cbuf[3];
  return(l_union.lval);
}

float 
float_swap(char cbuf[]) {
  union {
    char cval[4];
    float fval;
  } f_union;
  
  f_union.cval[3] = cbuf[0];
  f_union.cval[2] = cbuf[1];
  f_union.cval[1] = cbuf[2];
  f_union.cval[0] = cbuf[3];
  return(f_union.fval);
}

long 
long_join(char cbuf[]) {
  union {
    char cval[4];
    long lval;
  } l_union;
  
  l_union.cval[3] = cbuf[3];
  l_union.cval[2] = cbuf[2];
  l_union.cval[1] = cbuf[1];
  l_union.cval[0] = cbuf[0];
  return(l_union.lval);
}

