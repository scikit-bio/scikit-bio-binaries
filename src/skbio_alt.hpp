/*
 * Classes, methods and unction that provide skbio-like unctionality
 */

#ifndef UNIFRAC_SKBIO_ALT_H
#define UNIFRAC_SKBIO_ALT_H

#include <stdint.h>

namespace su {

// Compute Permanova
void permanova(const double * mat, unsigned int n_dims, const uint32_t *grouping, unsigned int n_perm, double &fstat_out, double &pvalue_out);
void permanova(const float  * mat, unsigned int n_dims, const uint32_t *grouping, unsigned int n_perm, float  &fstat_out, float  &pvalue_out);

}

#endif
