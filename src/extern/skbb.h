/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2025--, scikit-bio development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#ifndef SKBB_H
#define SKBB_H

#ifdef __cplusplus
#include <vector>
#define EXTERN extern "C"
#else
#define EXTERN
#endif

/*
 *  Test for significant differences between groups using PERMANOVA.
 *
 *  Permutational Multivariate Analysis of Variance (PERMANOVA) is a
 *  non-parametric method that tests whether two or more groups of objects
 *  (e.g., samples) are significantly different based on a categorical factor.
 *  It is conceptually similar to ANOVA except that it operates on a distance
 *  matrix, which allows for multivariate analysis. PERMANOVA computes a
 *  pseudo-F statistic.
 *
 *  Statistical significance is assessed via a permutation test. The assignment
 *  of objects to groups (`grouping`) is randomly permuted a number of times
 *  (controlled via `n_perm`). A pseudo-F statistic is computed for each
 *  permutation and the p-value is the proportion of permuted pseudo-F
 *  statisics that are equal to or greater than the original (unpermuted)
 *  pseudo-F statistic.
 *
 *  Input parameters:
 *   n_dims    - size of the matrix
 *   mat       - distance matrix (n_dims x n_dims)
 *   grouping  - grouping array of size n_dims 
 *   n_perm    - number of permutations
 *
 *  Output parameters, as pointers:
 *   fstat_ptr  - pseudo-F statistics
 *   pvalue_ptr - p-value
 */

EXTERN void skbb_permanova_fp64(unsigned int n_dims, const double * mat, const unsigned int *grouping, unsigned int n_perm, double *fstat_ptr, double *pvalue_ptr);

EXTERN void skbb_permanova_fp32(unsigned int n_dims, const float  * mat, const unsigned int *grouping, unsigned int n_perm, float  *fstat_ptr, float  *pvalue_ptr);


#endif
