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
 * Classes, methods and functions for PCoA
 */

#ifndef SKBB_PCOA_HPP
#define SKBB_PCOA_HPP

#include <stdint.h>

namespace skbb {

// Note: Helper function
// Center the matrix
// mat and center must be nxn and symmetric
// centered must be pre-allocated and same size as mat...will work even if centered==mat
void mat_to_centered(const double mat[], const uint32_t n_samples, double centered[]);
void mat_to_centered(const float  mat[], const uint32_t n_samples, float  centered[]);
void mat_to_centered(const double mat[], const uint32_t n_samples, float  centered[]);

// Note: Helper function
// Find eigen values and vectors
// Based on N. Halko, P.G. Martinsson, Y. Shkolnisky, and M. Tygert.
//     Original Paper: https://arxiv.org/abs/1007.5510
// centered == n x n, must be symmetric, Note: will be used in-place as temp buffer
void find_eigens_fast(const uint32_t n_samples, const uint32_t n_dims, double centered[], double * &eigenvalues, double * &eigenvectors);
void find_eigens_fast(const uint32_t n_samples, const uint32_t n_dims, float  centered[], float  * &eigenvalues, float  * &eigenvectors);

/*
 *  Perform Principal Coordinate Analysis.
 *
 *  Principal Coordinate Analysis (PCoA) is a method similar
 *  to Principal Components Analysis (PCA) with the difference that PCoA
 *  operates on distance matrices, typically with non-euclidian and thus
 *  ecologically meaningful distances like UniFrac in microbiome research.
 *
 *  In ecology, the euclidean distance preserved by Principal
 *  Component Analysis (PCA) is often not a good choice because it
 *  deals poorly with double zeros (Species have unimodal
 *  distributions along environmental gradients, so if a species is
 *  absent from two sites at the same site, it can't be known if an
 *  environmental variable is too high in one of them and too low in
 *  the other, or too low in both, etc. On the other hand, if an
 *  species is present in two sites, that means that the sites are
 *  similar.).
 *
 *  Implemented using FSVD method.
 *
 *  Note that the returned eigenvectors are not normalized to unit length.
 *
 *  Input parameters:
 *   n_dims    - Size of the matrix
 *   mat       - Distance matrix (n_dims x n_dims)
 *   n_eighs   - Number of eigenvalues to return
 *
 *  Output parameters: (returns a new buffer, must be released by caller)
 *   eigenvalues          - Array of size n_eighs
 *   samples              - The position of the samples in the ordination space,
 *                          row-indexed by the sample id.
 *                          Matrix of size (n_eighs x n_dims)
 *   proportion_explained - Array of size n_eighs.
 *                          The index corresponds to the ordination axis labels.
*/

// Perform Principal Coordinate Analysis using FSVD method
// n_samples - in, size of the matrix (n x n)
// mat       - in, 
// n_dims    - in, Dimensions to reduce the distance matrix to. This number determines how many eigenvectors and eigenvalues will be returned.
// eigenvalues - out, alocated buffer of size n_dims
// samples     - out, alocated buffer of size n_dims x n_samples
// proportion_explained - out, allocated buffer of size n_dims
void pcoa_fsvd(const uint32_t n_dims, const double mat[], const uint32_t n_eighs, double * &eigenvalues, double * &samples, double * &proportion_explained);
void pcoa_fsvd(const uint32_t n_dims, const float  mat[], const uint32_t n_eighs, float  * &eigenvalues, float  * &samples, float  * &proportion_explained);
void pcoa_fsvd(const uint32_t n_dims, const double mat[], const uint32_t n_eighs, float  * &eigenvalues, float  * &samples, float  * &proportion_explained);

// in-place version, will use mat as temp buffer internally
void pcoa_fsvd_inplace(const uint32_t n_dims, double mat[], const uint32_t n_eighs, double * &eigenvalues, double * &samples, double * &proportion_explained);
void pcoa_fsvd_inplace(const uint32_t n_dims, float  mat[], const uint32_t n_eighs, float  * &eigenvalues, float  * &samples, float  * &proportion_explained);

}

#endif
