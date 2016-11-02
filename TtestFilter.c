#include <stdio.h>  //fprintf
#include <stdlib.h> //free
#include <zlib.h>
#include <inttypes.h>
#include <math.h> // pow()
#include <stdint.h>

#define VERSION "0.0.1"

int main(int argc, char *argv[])
{
  //char *counts_file;
  //int k_length = 31;

  int c;
  while ((c = getopt(argc, argv, "k:")) >= 0) {
    switch (c) {
      //case 'k': k_length = atoi(optarg); break;
    }
  }

  if (optind == argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   ttestFilter [options] <counts.tsv>\n\n");
		//fprintf(stderr, "Options: -k INT    length of k-mers (max_value: 32) [%d]\n", k_length);
		fprintf(stderr, "\n");
		return 1;
	}

  return 0;
}
