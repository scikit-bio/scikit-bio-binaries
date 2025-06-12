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

namespace su {

// Compute Permanova
void permanova(const double * mat, unsigned int n_dims, const uint32_t *grouping, unsigned int n_perm, double &fstat_out, double &pvalue_out);
void permanova(const float  * mat, unsigned int n_dims, const uint32_t *grouping, unsigned int n_perm, float  &fstat_out, float  &pvalue_out);

}

#endif
