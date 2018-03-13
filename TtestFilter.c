#include <stdio.h>  //fprintf
#include <stdlib.h> //free
#include <iostream> //cout
#include <zlib.h>
#include <inttypes.h>
#include <math.h> // pow()
#include <stdint.h>

// C++
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/poisson.hpp>

// klib (H.li)
#include "kseq.h"
#include "kvec.h"
#include "kstring.h"

KSEQ_INIT(gzFile, gzread)

#define VERSION "0.0.2"
#define TMP_FILENAME ".TtestFilter.tmp"

typedef struct {
  char *kmer;
  size_t n;
  double *counts;
  double log2FC, meanA, meanB, sumA, sumB;
} kmer_test_t;

void kmer_test_write(kmer_test_t *k, gzFile fp) {
  size_t kmer_l = strlen(k->kmer);
  gzwrite(fp, &kmer_l,    sizeof(kmer_l));
  gzwrite(fp, k->kmer,   kmer_l);
  gzwrite(fp, &k->n,      sizeof(size_t));
  gzwrite(fp, k->counts, sizeof(double) * k->n);
  //gzwrite(fp, &k->pvalue, sizeof(double));
  gzwrite(fp, &k->log2FC, sizeof(double));
  gzwrite(fp, &k->meanA,  sizeof(double));
  gzwrite(fp, &k->meanB,  sizeof(double));
}

int kmer_test_read(kmer_test_t *k, gzFile fp) {
  size_t kmer_l=0;
  gzread(fp, &kmer_l,     sizeof(size_t));
//  printf("%zu\n",kmer_l);
  if(k->kmer)
    k->kmer = (char*)realloc(k->kmer, kmer_l + 1);
  else
    k->kmer = (char*)malloc(kmer_l + 1);
  k->kmer[kmer_l] = '\0';
  gzread(fp, k->kmer,    kmer_l);
  gzread(fp, &k->n,       sizeof(size_t));
  size_t counts_size = sizeof(double) * k->n;
  if(k->counts)
    k->counts = (double*)realloc(k->counts, counts_size);
  else
    k->counts = (double*)malloc(counts_size);
  gzread(fp, k->counts,  counts_size);
  gzread(fp, &k->log2FC,  sizeof(double));
  gzread(fp, &k->meanA,   sizeof(double));
  return gzread(fp, &k->meanB,   sizeof(double));
}

void kmer_test_destroy(kmer_test_t *k) {
  if(k->kmer)
    free(k->kmer);
  if(k->counts)
    free(k->counts);
}

double log2(double n) {
  // log(n)/log(2) is log2.
  return log(n) / log(2);
}

double sum(double *counts, int n) {
  double sum = 0;
  int i = 0;
  for (; i < n; i++) {
    sum += counts[i];
  }
  return sum;
}


double mean(double *counts, int n) {
  return sum(counts,n) / n;
}

double sd(double *counts, int n, double mean) {
  double sum = 0;
  int i = 0;
  for (; i < n; i++) {
    sum += pow(counts[i] - mean, 2);
  }
  return sqrt(sum / (n-1));
}

int cmp_reverse_double_pointers(const void * a, const void * b) {
  const double aa = **(const double **)a;
  const double bb = **(const double **)b;
  if (aa > bb) {
    return -1;
  }
  if (bb > aa) {
    return  1;
  }
  return 0;
}

double dmin(double a, double b) {
  if (a > b) {
    return b;
  } else {
    return a;
  }
}

double compute_t_test(double m1, double m2, int n1, int n2, double sd1, double sd2)
{
    // Compute freedom degrees t-stat and pvalue
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
        boost::math::students_t dist(df);
        pvalue = 2 * cdf(complement(dist, fabs(t_stat)));
    }
    return pvalue;
}

double logpoisson(double lambda, double k)
{
    boost::math::poisson_distribution <double> dist(lambda);
    return boost::math::log1p(pdf(dist,k)-1);
    
}


/* Poisson model
 * similar to HAWK
 * https://www.biorxiv.org/content/biorxiv/early/2017/05/23/141267.full.pdf
 * K1 and K2 are sums of counts of a kmer between conditions
 * N1 and N2 are total kmer count sums between conditions
 * both as per HAWK notations
 */
double compute_poisson_likelihood_ratio(double K1, double K2, int64_t N1, int64_t N2)
{
    double t_stat = 0;
    double ll1 = logpoisson(K1, K1) + logpoisson(K2, K2);
    double theta = (K1 + K2) / (N1 + N2); // TODO this ratio should be related to normalization factors, check it sometimes
    double ll2 = logpoisson(theta * N1, K1) + logpoisson(theta * N2, K2);
    double llr = ll1 - ll2;
    t_stat = 2* llr;
    int df=1;
    boost::math::chi_squared dist(df);
    double pvalue = 1.0-cdf(dist, t_stat);
    return pvalue;
}


int main(int argc, char *argv[])
{
  char *counts_file, *conditions_file, *conditionA, *conditionB;
  //int k_length = 31;
  double pvalue_threshold = 0.05;
  double log2fc_threshold = 0;
  char* raw_pval_file = NULL;
  bool use_t_test = true;

  int c;
  while ((c = getopt(argc, argv, "p:f:r:l")) >= 0) {
    switch (c) {
      case 'p': pvalue_threshold = atof(optarg); break;
      case 'f': log2fc_threshold = atof(optarg); break;
      case 'r': raw_pval_file    = optarg; break;
      case 'l': use_t_test       = false; break;
    }
  }

  if ((optind + 3) >= argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   ttestFilter [options] <counts.tsv> <conditions.tsv> <conditionA> <conditionB>\n\n");
		fprintf(stderr, "Options: -p FLOAT  max pvalue [%.2f]\n", pvalue_threshold);
    fprintf(stderr, "         -f FLOAT  min log2 fold change (absolute) [%.2f]\n", log2fc_threshold);
    fprintf(stderr, "         -r STR    file were raw pvalues will be printed\n");
    fprintf(stderr, "         -l        disregard the program name, do a Poisson likelihood ratio (HAWK model)\n");
		fprintf(stderr, "\n");
		return 1;
	}

  counts_file     = argv[optind++];
  conditions_file = argv[optind++];
  conditionA      = argv[optind++];
  conditionB      = argv[optind++];

  gzFile fp, tmp_fp, pval_fp = NULL;
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


  double n1 = kv_size(conditionA_indicies), n2 = kv_size(conditionB_indicies);
  double *a       = (double*)malloc(n1 * sizeof(double));
  double *b       = (double*)malloc(n2 * sizeof(double));
  double *log_a   = (double*)malloc(n1 * sizeof(double));
  double *log_b   = (double*)malloc(n2 * sizeof(double));
  kmer_test_t kmer_test;
  kvec_t(double) pvalues;
  size_t j, i, line = 1;

  kmer_test.n      = kv_size(samples);
  kmer_test.counts = (double*)malloc(kv_size(samples) * sizeof(double));
  kv_init(pvalues);

  // Friendly print
  fprintf(stderr, "Condition A: %s\n", conditionA);
  fprintf(stderr, "Condition B: %s\n", conditionB);

  fprintf(stderr, "Group A: ");
  for(j = 0; j < n1; j++) {
    if(j > 0) fprintf(stderr, ", ");
    fprintf(stderr, "%s", kv_A(samples, kv_A(conditionA_indicies, j)));
  }
  fprintf(stderr, "\n");

  fprintf(stderr, "Group B: ");
  for(j = 0; j < n2; j++) {
    if(j > 0) fprintf(stderr, ", ");
    fprintf(stderr, "%s", kv_A(samples, kv_A(conditionB_indicies, j)));
  }
  fprintf(stderr, "\n");

  fprintf(stderr, "Computing t-test for all k-mers\n");

  // Print output header line
  fprintf(stdout, "tag\tpvalue\tmeanA\tmeanB\tlog2FC");
  for(j = 0; j < kv_size(samples); j++) { fprintf(stdout, "\t%s", kv_A(samples, j));}
  fprintf(stdout, "\n");

  // Open counts file
  fp = gzopen(counts_file, "r");
  if(!fp) { fprintf(stderr, "Failed to open %s\n", counts_file); exit(EXIT_FAILURE); }
  ks = ks_init(fp);

  // Open tmp file to put computed values before correction of pvalues
  tmp_fp = gzopen(TMP_FILENAME, "w");
  if(!tmp_fp) { fprintf(stderr, "Failed to open %s\n", TMP_FILENAME); exit(EXIT_FAILURE); }

  //file with raw pvals
  if(raw_pval_file) {
    pval_fp = gzopen(raw_pval_file, "w");
    if(!pval_fp) { fprintf(stderr, "Failed to open %s\n", raw_pval_file); exit(EXIT_FAILURE); }
  }

  // do a first pass over all counts to get total kmer counts
  uint64_t N1=0, N2=0; // total kmer counts as per Hawk notation
  ks_getuntil(ks, KS_SEP_LINE, str, &dret);
  while(ks_getuntil(ks, KS_SEP_SPACE, str, &dret) >= 0) {
    char *kmer = ks_release(str);
    kmer_test.kmer = kmer;

    // load counts
    j = 0;
    while(j < kmer_test.n && ks_getuntil(ks, KS_SEP_SPACE, str, &dret) >= 0) {
      kmer_test.counts[j] = (double) atoi(str->s);
      j++;
    }

    for(i = 0; i < n1; i++) {
      int k = kv_A(conditionA_indicies, i);
      N1 += kmer_test.counts[k];
    }
    for(i = 0; i < n2; i++) {
      int k = kv_A(conditionB_indicies, i);
      N2 += kmer_test.counts[k] ;
    }
  }

  // ks_rewind(ks);  // not enough
  // rewind counts file
  ks_destroy(ks);
  gzclose(fp); 
  fp = gzopen(counts_file, "r");
  ks = ks_init(fp);
  
  // skip header line
  ks_getuntil(ks, KS_SEP_LINE, str, &dret);
  while(ks_getuntil(ks, KS_SEP_SPACE, str, &dret) >= 0) {
    char *kmer = ks_release(str);
    kmer_test.kmer = kmer;

    // load counts
    j = 0;
    while(j < kmer_test.n && ks_getuntil(ks, KS_SEP_SPACE, str, &dret) >= 0) {
      kmer_test.counts[j] = (double) atoi(str->s); 
      if (use_t_test)
        kmer_test.counts[j] /= normalization_factors[j]; 
      j++;
    }

    if (dret != '\n') {
      fprintf(stderr, "inconsistent number of column (line %zu)\n", line);
      exit(EXIT_FAILURE);
    }

    // Extract count for each conditions
    for(i = 0; i < n1; i++) {
      int k = kv_A(conditionA_indicies, i);
      a[i] =  kmer_test.counts[k] + 1; // FIXME that +1 fix isn't totally satisfactory, it artifically increases the counts (rayan)
      log_a[i] = log(kmer_test.counts[k] + 1);
    }
    for(i = 0; i < n2; i++) {
      int k = kv_A(conditionB_indicies, i);
      b[i]  = kmer_test.counts[k] + 1;
      log_b[i] = log(kmer_test.counts[k] + 1);
    }

    // Compute base means
    kmer_test.meanA = mean(a, n1);
    kmer_test.meanB = mean(b, n2);

    // Compute log means
    double m1 = mean(log_a, n1), m2 = mean(log_b, n2);

    // Compute log standard deviation
    double sd1 = sd(log_a, n1, m1), sd2 = sd(log_b, n2, m2);

    // Compute pvalues (a big chunk of computation is being done inside those functions)
    double pvalue;
    if (use_t_test)
    {
        pvalue = compute_t_test(m1,m2,n1,n2,sd1,sd2);
    }
    else
    {   // poisson
        kmer_test.sumA = mean(a, n1);
        kmer_test.sumB = mean(b, n2);
        pvalue = compute_poisson_likelihood_ratio(kmer_test.sumA, kmer_test.sumB, N1, N2);
    }

    if(raw_pval_file) {
      char buffer[100];
      sprintf(buffer,"%s\t%f\n",kmer_test.kmer,pvalue);
      gzwrite(pval_fp, &buffer, strlen(buffer));
    }

    // Compute log2FC
    kmer_test.log2FC = log2(kmer_test.meanB/kmer_test.meanA);

    // Write kmer to tmp file
    kmer_test_write(&kmer_test, tmp_fp);

    // Place pvalue into array for multi-test correction
    kv_push(double, pvalues, pvalue);

    // Free k-mer
    free(kmer);
    line++;
  }
  if(raw_pval_file)
    gzclose(pval_fp);
  ks_destroy(ks);
  gzclose(fp);
  gzclose(tmp_fp);

  free(a);
  free(b);
  free(log_a);
  free(log_b);

  fprintf(stderr, "Adjusting pvalues with Benjamini–Hochberg\n");

  // Correction of pvalues using Benjamini–Hochberg
  double **pvalues_indicies = (double**)malloc(sizeof(double*) * kv_size(pvalues));
  for(i = 0; i < kv_size(pvalues); i ++) {
    pvalues_indicies[i] = &pvalues.a[i];
  }
  qsort(pvalues_indicies, kv_size(pvalues), sizeof(double*), cmp_reverse_double_pointers);
  double min = 1;
  for(i = 0; i < kv_size(pvalues); i++) {
    double new_pvalue = (double)kv_size(pvalues) / ((double)kv_size(pvalues) - (double)i) * (*pvalues_indicies[i]);

    if(new_pvalue > min) {
      new_pvalue = min;
    } else {
      min = new_pvalue;
    }

    *pvalues_indicies[i] = dmin(1, new_pvalue);
  }
  free(pvalues_indicies);

  kmer_test.kmer   = NULL;

  fprintf(stderr, "Print the final output\n");

  // Open tmp file and print final ouput with corrected pvalues
  tmp_fp = gzopen(TMP_FILENAME, "r");
  if(!tmp_fp) { fprintf(stderr, "Failed to open %s\n", TMP_FILENAME); exit(EXIT_FAILURE); }
  i = 0;
  while(kmer_test_read(&kmer_test, tmp_fp)) {
    double pvalue = kv_A(pvalues, i);
    //std::cout << "pvalue: "<< pvalue << " log2FC: " << kmer_test.log2FC << std::endl;
    if(pvalue <= pvalue_threshold && fabs(kmer_test.log2FC) >= log2fc_threshold) {
      fprintf(stdout, "%s\t%f\t%f\t%f\t%f", kmer_test.kmer, pvalue, kmer_test.meanA, kmer_test.meanB, kmer_test.log2FC);
      for(j = 0; j < kmer_test.n; j++) { fprintf(stdout, "\t%.2f", kmer_test.counts[j]); }
      fprintf(stdout, "\n");
    }
    i++;
  }
  gzclose(tmp_fp);
  remove(TMP_FILENAME);
  kmer_test_destroy(&kmer_test);

  return 0;
}
