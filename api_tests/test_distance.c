/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2025--, scikit-bio development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#include "scikit-bio-binaries/distance.h"
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

int global_failed = false;
int failed = false;

#define ASSERT(x) {if(!(x)) \
    { fprintf(stderr, "failed assert [%s:%i]\n", __FILE__, __LINE__); \
      failed = true; global_failed = true;}}

void test_permanova_ties() {
    failed = false;

    // Same as tests in src/tests
    const double matrix_fp64[] = { 
      0., 1., 1., 4.,
      1., 0., 3., 2.,
      1., 3., 0., 3.,
      4., 2., 3., 0.};
    const float matrix_fp32[] = { 
      0., 1., 1., 4.,
      1., 0., 3., 2.,
      1., 3., 0., 3.,
      4., 2., 3., 0.};

    const unsigned int grouping_equal[] = { 0, 0, 1, 1};
    const unsigned int grouping_equal_swapped[] = { 1, 1, 0, 0};

    const unsigned int  n_samples = 4;


    // value from skbio
    const float exp_stat = 2.0;
    // using a different random than skbio, so result different
    const float exp_pvalue = 0.66;

    double stat_fp64, pvalue_fp64;
    float stat_fp32, pvalue_fp32;

    skbb_permanova_fp64(n_samples, matrix_fp64,
                  grouping_equal, 999,
                  &stat_fp64, &pvalue_fp64);
    ASSERT(fabs(stat_fp64 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp64 - exp_pvalue) < 0.05);

    skbb_permanova_fp32(n_samples, matrix_fp32,
                  grouping_equal, 999,
                  &stat_fp32, &pvalue_fp32);
    ASSERT(fabs(stat_fp64 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp64 - exp_pvalue) < 0.05);

    skbb_permanova_fp32(n_samples, matrix_fp32,
                  grouping_equal_swapped, 999,
                  &stat_fp32, &pvalue_fp32);
    ASSERT(fabs(stat_fp32 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp32 - exp_pvalue) < 0.05);

    skbb_permanova_fp64(n_samples, matrix_fp64,
                  grouping_equal_swapped, 999,
                  &stat_fp64, &pvalue_fp64);
    ASSERT(fabs(stat_fp64 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp64 - exp_pvalue) < 0.05);

    if (failed) {
      printf("ERROR: Permanova test failed\n");
    } else {
      printf("INFO: Permanova test succeeded\n");
    }
}

int main(int argc, char** argv) {
    global_failed = false;
    test_permanova_ties();

    printf("\n");
    if (global_failed) {
      printf("ERROR: Some distance tests failed\n");
    } else {
      printf("INFO: All distance tests succeeded\n");
    }
    
    return failed ? 1 : 0;
}
