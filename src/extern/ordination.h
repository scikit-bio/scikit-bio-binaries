/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2025--, scikit-bio development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#ifndef SKBB_EXTERN_ORDINATION_H
#define SKBB_EXTERN_ORDINATION_H

#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN
#endif

/*
 *  Centers a distance matrix.
 *
 *  Following computation as per
 *  Numerical Ecology (Legendre & Legendre 1998).
 *
 *  Input parameters:
 *   n_dims    - Size of the matrix
 *   mat       - Distance matrix (n_dims x n_dims)
 *
 *  Output parameters:
 *   centered  - Resulting matrix (n_dims x n_dims), pre-allocated
 *
 *   Note: mat and centered are allowed to point to the same memory area (e.g. in-place)
 */

EXTERN void skbb_center_distance_matrix_fp64(unsigned int n_dims, const double mat[], double centered[]);
EXTERN void skbb_center_distance_matrix_fp32(unsigned int n_dims, const float  mat[], float  centered[]);
EXTERN void skbb_center_distance_matrix_fp64_to_fp32(unsigned int n_dims, const double mat[], float  centered[]);

/*
 * Performs singular value decomposition, or more specifically in
 * this case eigendecomposition, using fast heuristic algorithm
 * nicknamed "FSVD" (FastSVD), adapted and optimized from the algorithm
 * described by Halko et al (2011).
 *   N. Halko, P.G. Martinsson, Y. Shkolnisky, and M. Tygert.
 *  Original Paper: https://arxiv.org/abs/1007.5510
 *
 *  Input parameters:
 *   n_dims    - Size of the matrix
 *   centered  - Centered distance matrix (n_dims x n_dims), will be overwritten during compute
 *   n_eighs   - Number of eigenvalues to return
 *
 *  Output parameters:
 *   eigenvalues          - Array of size n_eighs, pre-allocated
 *   eigenvectors         - Matrix of size (n_dims x n_eighs), pre-allocated
 */
EXTERN void skbb_fsvd_inplace_fp64(unsigned int n_dims, double centered[], unsigned int n_eighs, double eigenvalues[], double eigenvectors[]);
EXTERN void skbb_fsvd_inplace_fp32(unsigned int n_dims, float  centered[], unsigned int n_eighs, float  eigenvalues[], float  eigenvectors[]);

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
 *  Output parameters:
 *   eigenvalues          - Array of size n_eighs, pre-allocated
 *   samples              - The position of the samples in the ordination space,
 *                          row-indexed by the sample id.
 *                          Matrix of size (n_eighs x n_dims), pre-allocated
 *   proportion_explained - Array of size n_eighs, pre-allocated
 *                          The index corresponds to the ordination axis labels.
*/

EXTERN void skbb_pcoa_fsvd_fp64(unsigned int n_dims, const double mat[], unsigned int n_eighs, double eigenvalues[], double samples[], double proportion_explained[]);
EXTERN void skbb_pcoa_fsvd_fp32(unsigned int n_dims, const float  mat[], unsigned int n_eighs, float  eigenvalues[], float  samples[], float  proportion_explained[]);
EXTERN void skbb_pcoa_fsvd_fp64_to_fp32(unsigned int n_dims, const double mat[], unsigned int n_eighs, float  eigenvalues[], float  samples[], float  proportion_explained[]);

// in-place version, will use mat as temp buffer internally
EXTERN void skbb_pcoa_fsvd_inplace_fp64(unsigned int n_dims, double mat[], unsigned int n_eighs, double eigenvalues[], double samples[], double proportion_explained[]);
EXTERN void skbb_pcoa_fsvd_inplace_fp32(unsigned int n_dims, float  mat[], unsigned int n_eighs, float  eigenvalues[], float  samples[], float  proportion_explained[]);

#endif
