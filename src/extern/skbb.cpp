/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2025--, scikit-bio development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#include "skbb.h"

#include "distance/permanova.hpp"

void skbb_permanova_fp64(unsigned int n_dims, const double mat[], const unsigned int grouping[], unsigned int n_perm, double *fstat_ptr, double *pvalue_ptr) {
  skbb::permanova(n_dims, mat, grouping, n_perm, *fstat_ptr, *pvalue_ptr);
}

void skbb_permanova_fp32(unsigned int n_dims, const float mat[], const unsigned int grouping[], unsigned int n_perm, float *fstat_ptr, float *pvalue_ptr) {
  skbb::permanova(n_dims, mat, grouping, n_perm, *fstat_ptr, *pvalue_ptr);
}

