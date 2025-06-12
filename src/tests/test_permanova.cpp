/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2025, UniFrac development team.
 * Copyright (c) 2025--, scikit-bio development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#include "distance/permanova.hpp"
#include <unistd.h>
#include "tests/test_helper.hpp"

void test_permanova_ties() {
    SUITE_START("test permanova ties");

    // Same as skbio tests
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

    const uint32_t grouping_equal[] = { 0, 0, 1, 1};
    const uint32_t grouping_equal_swapped[] = { 1, 1, 0, 0};

    const uint32_t n_samples = 4;


    // value from skbio
    const float exp_stat = 2.0;
    // using a different random than skbio, so result different
    const float exp_pvalue = 0.66;

    double stat_fp64, pvalue_fp64;
    float stat_fp32, pvalue_fp32;

    su::permanova(matrix_fp64, n_samples,
                  grouping_equal, 999,
                  stat_fp64, pvalue_fp64);
    ASSERT(fabs(stat_fp64 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp64 - exp_pvalue) < 0.05);

    su::permanova(matrix_fp32, n_samples,
                  grouping_equal, 999,
                  stat_fp32, pvalue_fp32);
    ASSERT(fabs(stat_fp64 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp64 - exp_pvalue) < 0.05);

    su::permanova(matrix_fp32, n_samples,
                  grouping_equal_swapped, 999,
                  stat_fp32, pvalue_fp32);
    ASSERT(fabs(stat_fp32 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp32 - exp_pvalue) < 0.05);

    su::permanova(matrix_fp64, n_samples,
                  grouping_equal_swapped, 999,
                  stat_fp64, pvalue_fp64);
    ASSERT(fabs(stat_fp64 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp64 - exp_pvalue) < 0.05);

    SUITE_END();
}

void test_permanova_noties() {
    SUITE_START("test permanova noties");

    // Same as skbio tests
    const double matrix_fp64[] = { 
      0., 1., 5., 4.,
      5., 0., 3., 2.,
      1., 3., 0., 3.,
      4., 2., 3., 0.};
    const float matrix_fp32[] = { 
      0., 1., 5., 4.,
      1., 0., 3., 2.,
      5., 3., 0., 3.,
      4., 2., 3., 0.};

    const uint32_t grouping_equal[] = { 0, 0, 1, 1};
    const uint32_t grouping_equal_swapped[] = { 1, 1, 0, 0};

    const uint32_t n_samples = 4;


    // value from skbio
    const float exp_stat = 4.4;
    // using a different random than skbio, so result different
    const float exp_pvalue = 0.34;

    double stat_fp64, pvalue_fp64;
    float stat_fp32, pvalue_fp32;

    su::permanova(matrix_fp64, n_samples,
                  grouping_equal, 999,
                  stat_fp64, pvalue_fp64);
    ASSERT(fabs(stat_fp64 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp64 - exp_pvalue) < 0.05);

    su::permanova(matrix_fp32, n_samples,
                  grouping_equal, 999,
                  stat_fp32, pvalue_fp32);
    ASSERT(fabs(stat_fp32 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp32 - exp_pvalue) < 0.05);

    su::permanova(matrix_fp32, n_samples,
                  grouping_equal_swapped, 999,
                  stat_fp32, pvalue_fp32);
    ASSERT(fabs(stat_fp32 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp32 - exp_pvalue) < 0.05);

    su::permanova(matrix_fp64, n_samples,
                  grouping_equal_swapped, 999,
                  stat_fp64, pvalue_fp64);
    ASSERT(fabs(stat_fp64 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp64 - exp_pvalue) < 0.05);

    SUITE_END();
}

void test_permanova_unequal() {
    SUITE_START("test permanova unequal");

    // Same as skbio tests
    const double matrix_fp64[] = { 
      0.0,    1.0,   0.1,   0.5678, 1.0,   1.0,
      1.0,    0.0,   0.002, 0.42,   0.998, 0.0,
      0.1,    0.002, 0.0,   1.0,    0.123, 1.0,
      0.5678, 0.42,  1.0,   0.0,    0.123, 0.43,
      1.0,    0.998, 0.123, 0.123,  0.0,   0.5,
      1.0,    0.0,   1.0,   0.43,   0.5,   0.0 };
    const float matrix_fp32[] = { 
      0.0,    1.0,   0.1,   0.5678, 1.0,   1.0,
      1.0,    0.0,   0.002, 0.42,   0.998, 0.0,
      0.1,    0.002, 0.0,   1.0,    0.123, 1.0,
      0.5678, 0.42,  1.0,   0.0,    0.123, 0.43,
      1.0,    0.998, 0.123, 0.123,  0.0,   0.5,
      1.0,    0.0,   1.0,   0.43,   0.5,   0.0 };

    const uint32_t grouping_equal[] = { 0, 1, 2, 1, 0 , 0};
    const uint32_t grouping_equal_swapped[] = { 1, 2, 0, 2, 1, 1};

    const uint32_t n_samples = 6;


    // value from skbio
    const float exp_stat = 0.578848;
    // using a different random than skbio, so result different
    const float exp_pvalue = 0.65;

    double stat_fp64, pvalue_fp64;
    float stat_fp32, pvalue_fp32;

    su::permanova(matrix_fp64, n_samples,
                  grouping_equal, 999,
                  stat_fp64, pvalue_fp64);
    ASSERT(fabs(stat_fp64 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp64 - exp_pvalue) < 0.05);

    su::permanova(matrix_fp32, n_samples,
                  grouping_equal, 999,
                  stat_fp32, pvalue_fp32);
    ASSERT(fabs(stat_fp32 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp32 - exp_pvalue) < 0.05);

    su::permanova(matrix_fp32, n_samples,
                  grouping_equal, 999,
                  stat_fp32, pvalue_fp32);
    ASSERT(fabs(stat_fp32 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp32 - exp_pvalue) < 0.05);

    su::permanova(matrix_fp64, n_samples,
                  grouping_equal_swapped, 999,
                  stat_fp64, pvalue_fp64);
    ASSERT(fabs(stat_fp64 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp64 - exp_pvalue) < 0.05);

    SUITE_END();
}

int main(int argc, char** argv) {
    test_permanova_ties();
    test_permanova_noties();
    test_permanova_unequal();

    printf("\n");
    printf(" %i / %i suites failed\n", suites_failed, suites_run);
    printf(" %i / %i suites empty\n", suites_empty, suites_run);
    printf(" %i / %i tests failed\n", tests_failed, tests_run);
  
    printf("\n THE END.\n");
    
    return tests_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
