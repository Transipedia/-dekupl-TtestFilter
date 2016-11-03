#include <stdio.h>  //fprintf
#include <stdlib.h> //free
#include <zlib.h>
#include <inttypes.h>
#include <math.h> // pow()
#include <stdint.h>

// C++
#include <boost/math/distributions/students_t.hpp>
#include <iostream>
#include <iomanip>

// klib (H.li)
#include "kseq.h"
#include "kvec.h"
#include "kstring.h"

KSEQ_INIT(gzFile, gzread)

#define VERSION "0.0.1"

// Simplify usage:
using namespace boost::math;
using namespace std;

double log2(double n) {
  // log(n)/log(2) is log2.
  return log(n) / log(2);
}

double mean(double *counts, int n) {
  double sum = 0;
  for (int i = 0; i < n; i++) {
    sum += counts[i];
  }
  return sum / n;
}

double sd(double *counts, int n, double mean) {
  double sum = 0;
  for (int i = 0; i < n; i++) {
    sum += pow(counts[i] - mean, 2);
  }
  return sqrt(sum / (n-1));
}

int main(int argc, char *argv[])
{
  char *counts_file, *conditions_file, *conditionA, *conditionB;
  //int k_length = 31;
  double pvalue_threshold = 0.05;
  double log2fc_threshold = 1;

  int c;
  while ((c = getopt(argc, argv, "p:f:")) >= 0) {
    switch (c) {
      case 'p': pvalue_threshold = atof(optarg); break;
      case 'f': log2fc_threshold = atof(optarg); break;
    }
  }

  if ((optind + 3) >= argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   ttestFilter [options] <counts.tsv> <conditions.tsv> <conditionA> <conditionB>\n\n");
		fprintf(stderr, "Options: -p FLOAT  max pvalue [%.2f]\n", pvalue_threshold);
    fprintf(stderr, "         -f FLOAT  min log2 fold change (absolute) [%.2f]\n", log2fc_threshold);
		fprintf(stderr, "\n");
		return 1;
	}

  counts_file     = argv[optind++];
  conditions_file = argv[optind++];
  conditionA      = argv[optind++];
  conditionB      = argv[optind++];

  gzFile fp;
	kstream_t *ks;
	kstring_t *str;
  kvec_t(char*) samples;
  kvec_t(int) conditionA_indicies, conditionB_indicies;
  int dret = 0;

  str = (kstring_t*)calloc(1, sizeof(kstring_t));
  kv_init(samples);
  kv_init(conditionA_indicies);
  kv_init(conditionB_indicies);

  // 1. Get samples names from counts_file

  fp = gzopen(counts_file, "r");
  if(!fp) { fprintf(stderr, "Failed to open %s\n", counts_file); exit(EXIT_FAILURE); }
  ks = ks_init(fp);

  ks_getuntil(ks, KS_SEP_SPACE, str, &dret);
  while(ks_getuntil(ks, KS_SEP_SPACE, str, &dret) >= 0) {
    kv_push(char*, samples, ks_release(str));
    if (dret == '\n') break;
  }
  ks_destroy(ks);
  gzclose(fp);

  // 1. Get samples conditions indicies and Normalization factors

  double *normalization_factors = (double*)calloc(kv_size(samples), sizeof(double));

  fp = gzopen(conditions_file, "r");
  if(!fp) { fprintf(stderr, "Failed to open %s\n", conditions_file); exit(EXIT_FAILURE); }

  ks = ks_init(fp);

  // Skip header line
  ks_getuntil(ks, KS_SEP_LINE, str, &dret);
  while(ks_getuntil(ks, KS_SEP_SPACE, str, &dret) >= 0) {
    char *sample_name = ks_release(str);
    if(dret != '\n' && ks_getuntil(ks, 0, str, &dret) > 0) {
      char *condition = ks_release(str);
      if(dret != '\n' && ks_getuntil(ks, 0, str, &dret) > 0) {
        double norm_fact = atof(str->s);
        size_t j = 0;
        while(j < kv_size(samples) && strcmp(kv_A(samples,j),sample_name) != 0) { j++; }
        if(j < kv_size(samples)) {
          if(strcmp(condition,conditionA) == 0) {
            kv_push(int, conditionA_indicies, j);
          } else if(strcmp(condition, conditionB) == 0) {
            kv_push(int, conditionB_indicies, j);
          }
          normalization_factors[j] = norm_fact;
        }
      }
      free(condition);
    }
    free(sample_name);
    // skip the rest of the line
		if (dret != '\n') while ((dret = ks_getc(ks)) > 0 && dret != '\n');
  }
  ks_destroy(ks);
  gzclose(fp);

  fprintf(stderr, "CONDITION_A samples: ");
  for(size_t j = 0; j < kv_size(conditionA_indicies); j++) {
    if(j > 0) fprintf(stderr, ", ");
    fprintf(stderr, "%d => %s", kv_A(conditionA_indicies, j), kv_A(samples, kv_A(conditionA_indicies, j)));
  }
  fprintf(stderr, "\n");

  fprintf(stderr, "CONDITION_B samples: ");
  for(size_t j = 0; j < kv_size(conditionB_indicies); j++) {
    if(j > 0) fprintf(stderr, ", ");
    fprintf(stderr, "%d => %s", kv_A(conditionB_indicies, j), kv_A(samples, kv_A(conditionB_indicies, j)));
  }
  fprintf(stderr, "\n");

  double n1 = kv_size(conditionA_indicies), n2 = kv_size(conditionB_indicies);
  double *counts  = (double*)malloc(kv_size(samples) * sizeof(double));
  double *a       = (double*)malloc(n1 * sizeof(double));
  double *b       = (double*)malloc(n2 * sizeof(double));

  // Print output header line
  fprintf(stdout, "tag\tpvalue\tmeanA\tmeanB\tlog2FC");
  for(size_t j = 0; j < kv_size(samples); j++) { fprintf(stdout, "\t%s", kv_A(samples, j));}
  fprintf(stdout, "\n");

  fp = gzopen(counts_file, "r");
  if(!fp) { fprintf(stderr, "Failed to open %s\n", counts_file); exit(EXIT_FAILURE); }
  ks = ks_init(fp);

  // skip header line
  ks_getuntil(ks, KS_SEP_LINE, str, &dret);
  while(ks_getuntil(ks, KS_SEP_SPACE, str, &dret) >= 0) {
    char *kmer = ks_release(str);

    // load counts
    size_t j = 0;
    while(ks_getuntil(ks, KS_SEP_SPACE, str, &dret) >= 0 && j < kv_size(samples)) {
      counts[j] = (double) atoi(str->s) / normalization_factors[j];
      j++;
    }

    if (dret != '\n') {
      fprintf(stderr, "inconsistent number of column\n");
      exit(EXIT_FAILURE);
    }

    for(int i = 0; i < n1; i++) {
      int k = kv_A(conditionA_indicies, i);
      a[i] =  counts[k];
    }
    for(int i = 0; i < n2; i++) {
      int k = kv_A(conditionB_indicies, i);
      b[i]  = counts[k];
    }

    // Compute means
    double m1 = mean(a, n1), m2 = mean(b, n2);

    // Compute standard deviation
    double sd1 = sd(a, n1, m1), sd2 = sd(b, n2, m2);

    double df, t_stat, pvalue;

    if(sd1 == 0 && sd2 == 0) {
      df = 1;
      t_stat = -1;
      pvalue = 0.5;
    } else {
      // Compute freedom degree
      df = sd1 * sd1 / n1 + sd2 * sd2 / n2;
      df *= df;
      double t1 = sd1 * sd1 / n1;
      t1 *= t1;
      t1 /=  (n1 - 1);
      double t2 = sd2 * sd2 / n2;
      t2 *= t2;
      t2 /= (n2 - 1);
      df /= (t1 + t2);

      // Compute t-statistic:
      t_stat = (m1 - m2) / sqrt(sd1 * sd1 / n1 + sd2 * sd2 / n2);

      // Compute p-value
      students_t dist(df);
      pvalue = 2 * cdf(complement(dist, fabs(t_stat)));
    }

    // Compute log2FC
    double log2FC = log2(m2/m1);

    if(pvalue <= pvalue_threshold || fabs(log2FC)) {
      fprintf(stdout, "%s\t%f\t%f\t%f\t%f", kmer, pvalue, m1, m2, log2FC);
      for(size_t i = 0; i < kv_size(samples); i++) { fprintf(stdout, "\t%.2f", counts[i]); }
      fprintf(stdout, "\n");
    }

    // Free k-mer
    free(kmer);
  }
  ks_destroy(ks);
  gzclose(fp);

  return 0;
}
