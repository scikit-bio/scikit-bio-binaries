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
 *
 * This file is used to create the necessary interfaces
 * by means of
 *   generate_permanova_dyn.py
 *
 * Anything ending in _T will get a acc-specific function wrapper.
 *
 */

#include "distance/permanova_dyn.hpp"
#include <cstdlib>
#include <algorithm>

#if !(defined(_OPENACC) || defined(OMPGPU))

#include <omp.h>

#elif defined(_OPENACC)

#include <openacc.h>

#endif

static inline int pmn_get_max_parallelism_T() {
#if !(defined(_OPENACC) || defined(OMPGPU))
  // No good reason to do more than max threads
  // (but use 2x to reduce thread spawning overhead)
  // but we do use 4x blocking, so account for that, too
  return 2*omp_get_max_threads()*4;
#else
  // 1k is enough for consumer-grade GPUs
  // 4k needed for larger GPUs
  // Keeping it at 4k not slowing down consumer GPUs
  // Getting much higher than 4k seems to slow down the compute
  return 4000; // should likely by dynamic, but there are no portable functions avaialble
#endif
}

#if !(defined(_OPENACC) || defined(OMPGPU))
// CPU version, uses tilling
// no internal parallelism, all parallelism handled outside

// Compute PERMANOVA pseudo-F partial statistic
// mat is symmetric matrix of size n_dims x in_n
// grouping is an array of size in_n
// inv_group_sizes is an array of size maxel(grouping)
template<class TFloat>
static inline TFloat pmn_f_stat_sW_one(
		const uint32_t n_dims,
		const TFloat * mat,
		const uint32_t *grouping,
		const TFloat *inv_group_sizes) {
  // Use full precision for intermediate compute, to minimize accumulation errors
  double s_W = 0.0;

  constexpr uint32_t TILE = 128;  // 128 is big enough to speed up the code, without cache trashing

  for (uint32_t trow=0; trow < (n_dims-1); trow+=TILE) {   // no columns in last row
    for (uint32_t tcol=trow+1; tcol < n_dims; tcol+=TILE) { // diagonal is always zero
      const uint32_t max_row = std::min(trow+TILE,n_dims-1);
      const uint32_t max_col = std::min(tcol+TILE,n_dims);

      for (uint32_t row=trow; row < max_row; row++) {
        const uint32_t min_col = std::max(tcol,row+1);
        uint32_t group_idx = grouping[row];

        // Use full precision for intermediate compute, to minimize accumulation errors
        double local_s_W = 0.0;
        const TFloat * mat_row = mat + uint64_t(row)*uint64_t(n_dims);
        for (uint32_t col=min_col; col < max_col; col++) {
            if (grouping[col] == group_idx) {
                TFloat val = mat_row[col];  // mat[row,col];
                local_s_W += val * val;
            }
        }
        s_W += local_s_W*inv_group_sizes[group_idx];
      }
    }
  }

  return s_W;
}


// 4-at-a-time version, to minimize memory reads while still fitting in registers
template<class TFloat>
static inline void pmn_f_stat_sW_four(
		const uint32_t n_dims,
		const TFloat * mat,
		const uint32_t *grouping0,
		const uint32_t *grouping1,
		const uint32_t *grouping2,
		const uint32_t *grouping3,
		const TFloat *inv_group_sizes,
		TFloat& out_sW0, TFloat& out_sW1, TFloat& out_sW2, TFloat& out_sW3) {
  // Use full precision for intermediate compute, to minimize accumulation errors
  double s_W0 = 0.0;
  double s_W1 = 0.0;
  double s_W2 = 0.0;
  double s_W3 = 0.0;

  constexpr uint32_t TILE = 128;  // 128 is big enough to speed up the code, without cache trashing

  for (uint32_t trow=0; trow < (n_dims-1); trow+=TILE) {   // no columns in last row
    for (uint32_t tcol=trow+1; tcol < n_dims; tcol+=TILE) { // diagonal is always zero
      const uint32_t max_row = std::min(trow+TILE,n_dims-1);
      const uint32_t max_col = std::min(tcol+TILE,n_dims);

      for (uint32_t row=trow; row < max_row; row++) {
        const uint32_t min_col = std::max(tcol,row+1);
        uint32_t group_idx0 = grouping0[row];
        uint32_t group_idx1 = grouping1[row];
        uint32_t group_idx2 = grouping2[row];
        uint32_t group_idx3 = grouping3[row];

        // Use full precision for intermediate compute, to minimize accumulation errors
        double local_s_W0 = 0.0;
        double local_s_W1 = 0.0;
        double local_s_W2 = 0.0;
        double local_s_W3 = 0.0;
        const TFloat * mat_row = mat + uint64_t(row)*uint64_t(n_dims);
        for (uint32_t col=min_col; col < max_col; col++) {
	    // speculatively read, we will likely use it at least in one of the ifs
	    // small penalty if we mis-predicted, no effect on semantics
            TFloat val = mat_row[col];  // mat[row,col];
	    val = val*val;
            if (grouping0[col] == group_idx0) {
                local_s_W0 += val;
            }
            if (grouping1[col] == group_idx1) {
                local_s_W1 += val;
            }
            if (grouping2[col] == group_idx2) {
                local_s_W2 += val;
            }
            if (grouping3[col] == group_idx3) {
                local_s_W3 += val;
            }
        }
        s_W0 += local_s_W0*inv_group_sizes[group_idx0];
        s_W1 += local_s_W1*inv_group_sizes[group_idx1];
        s_W2 += local_s_W2*inv_group_sizes[group_idx2];
        s_W3 += local_s_W3*inv_group_sizes[group_idx3];
      }
    }
  }

  out_sW0 = s_W0;
  out_sW1 = s_W1;
  out_sW2 = s_W2;
  out_sW3 = s_W3;
}
#endif

// Compute PERMANOVA pseudo-F partial statistic
// mat is symmetric matrix of size n_dims x n_dims
// groupings is a matrix of size n_dims x n_grouping_dims
// inv_group_sizes is an array of size maxel(groupings)
// Results in group_sWs, and array of size n_grouping_dims

#if !(defined(_OPENACC) || defined(OMPGPU))

template<class TFloat>
static inline void pmn_f_stat_sW_cpu(
		const uint32_t n_dims,
		const TFloat * mat,
		const uint32_t n_grouping_dims,
		const uint32_t *groupings,
		const TFloat *inv_group_sizes,
		TFloat *group_sWs) {
// CPU version, call function
#pragma omp parallel for
 for (uint32_t gblock=0; gblock < n_grouping_dims; gblock+=4) {
    if ((gblock+3)<n_grouping_dims) {
      // can process multiple permutations per pass
      const uint32_t grouping_el = gblock;
      const uint32_t *grouping0 = groupings + uint64_t(grouping_el)*uint64_t(n_dims);
      const uint32_t *grouping1 = groupings + uint64_t(grouping_el+1)*uint64_t(n_dims);
      const uint32_t *grouping2 = groupings + uint64_t(grouping_el+2)*uint64_t(n_dims);
      const uint32_t *grouping3 = groupings + uint64_t(grouping_el+3)*uint64_t(n_dims);
      pmn_f_stat_sW_four(n_dims,mat,
		    grouping0, grouping1, grouping2, grouping3,
		    inv_group_sizes,
		    group_sWs[grouping_el], group_sWs[grouping_el+1], group_sWs[grouping_el+2], group_sWs[grouping_el+3]);
    } else {
      // just do one at a time
      for (uint32_t grouping_el=gblock; grouping_el < n_grouping_dims; grouping_el++) {
        const uint32_t *grouping = groupings + uint64_t(grouping_el)*uint64_t(n_dims);
        group_sWs[grouping_el] = pmn_f_stat_sW_one<TFloat>(n_dims,mat,grouping,inv_group_sizes);
      }
    }
 } 
}

#else

template<class TFloat>
static inline void pmn_f_stat_sW_gpu(
		const uint32_t n_dims,
		const TFloat * mat,
		const uint32_t n_grouping_dims,
		const uint32_t *groupings,
		const TFloat *inv_group_sizes,
		TFloat *group_sWs) {
// GPU version, just put all the code in here
 const uint64_t groupings_size = uint64_t(n_dims)*uint64_t(n_grouping_dims);
#ifdef OMPGPU
#pragma omp target teams distribute map(to:groupings[0:groupings_size]) map(from:group_sWs[0:n_grouping_dims])
#else
#pragma acc parallel loop gang copyin(groupings[0:groupings_size]) copyout(group_sWs[0:n_grouping_dims]) default(present)
#endif
 for (uint32_t grouping_el=0; grouping_el < n_grouping_dims; grouping_el++) {
    const uint32_t *grouping = groupings + uint64_t(grouping_el)*uint64_t(n_dims);
    // Use full precision for intermediate compute, to minimize accumulation errors
    double s_W = 0.0;
    for (uint32_t row=0; row < (n_dims-1); row++) {   // no columns in last row
      uint32_t group_idx = grouping[row];
#ifdef OMPGPU
#pragma omp parallel for simd reduction(+:s_W)
#else
#pragma acc loop vector reduction(+:s_W)
#endif
      for (uint32_t col=row+1; col < n_dims; col++) { // diagonal is always zero
        if (grouping[col] == group_idx) {
            const TFloat * mat_row = mat + uint64_t(row)*uint64_t(n_dims);
            TFloat val = mat_row[col];  // mat[row,col];
            s_W += val * val * inv_group_sizes[group_idx];
        }
      }
    }
    group_sWs[grouping_el] = s_W;
 } 
}
#endif

template<class TFloat>
static inline void pmn_f_stat_sW_T(
		const uint32_t n_dims,
		const TFloat * mat,
		const uint32_t n_grouping_dims,
		const uint32_t *groupings,
		const TFloat *inv_group_sizes,
		TFloat *group_sWs) {
#if !(defined(_OPENACC) || defined(OMPGPU))
  pmn_f_stat_sW_cpu(n_dims, mat, n_grouping_dims, groupings, inv_group_sizes, group_sWs);
#else
  pmn_f_stat_sW_gpu(n_dims, mat, n_grouping_dims, groupings, inv_group_sizes, group_sWs);
#endif
}
