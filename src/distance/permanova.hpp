/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2025, UniFrac development team.
 * Copyright (c) 2025--, scikit-bio development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 * Classes, methods and functions for computing permanova
 */

#ifndef PERMANOVA_HPP
#define PERMANOVA_HPP

#include <stdint.h>

namespace skbb {

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
 *   seed      - Optional random seed, if non-negative. Use system random seed if <0
 *
 *  Output parameters:
 *   fstat_out  - pseudo-F statistics
 *   pvalue_out - p-value
 */

void permanova(unsigned int n_dims, const double mat[], const uint32_t grouping[], unsigned int n_perm, int seed, double &fstat_out, double &pvalue_out);
void permanova(unsigned int n_dims, const float  mat[], const uint32_t grouping[], unsigned int n_perm, int seed, float  &fstat_out, float  &pvalue_out);

}

#endif
