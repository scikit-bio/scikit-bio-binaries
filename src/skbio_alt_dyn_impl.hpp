/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2025, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#include "unifrac_accapi.hpp"
#include <cstdlib>
#include <algorithm>

#if !(defined(_OPENACC) || defined(OMPGPU))

#include <omp.h>

#elif defined(_OPENACC)

#include <openacc.h>

#endif

static inline int pmn_get_max_parallelism_T() {
#if !(defined(_OPENACC) || defined(OMPGPU))
  // no good reason to do more than max threads
  return omp_get_max_threads();
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
		const TFloat * mat,
		const uint32_t n_dims,
		const uint32_t *grouping,
		const TFloat *inv_group_sizes) {
  // Use full precision for intermediate compute, to minimize accumulation errors
  double s_W = 0.0;

  constexpr uint32_t TILE = 1024;  // 1k grouping els fit nicely in L1 cache (4kB total)

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
#endif

// Compute PERMANOVA pseudo-F partial statistic
// mat is symmetric matrix of size n_dims x n_dims
// groupings is a matrix of size n_dims x n_grouping_dims
// inv_group_sizes is an array of size maxel(groupings)
// Results in group_sWs, and array of size n_grouping_dims
// Note: Best results when n_grouping_dims fits in L1 cache
template<class TFloat>
static inline void pmn_f_stat_sW_T(
		const TFloat * mat,
		const uint32_t n_dims,
		const uint32_t *groupings,
		const uint32_t n_grouping_dims,
		const TFloat *inv_group_sizes,
		TFloat *group_sWs) {
#if !(defined(_OPENACC) || defined(OMPGPU))
// CPU version, call function
#pragma omp parallel for
 for (uint32_t grouping_el=0; grouping_el < n_grouping_dims; grouping_el++) {
    const uint32_t *grouping = groupings + uint64_t(grouping_el)*uint64_t(n_dims);
    group_sWs[grouping_el] = pmn_f_stat_sW_one<TFloat>(mat,n_dims,grouping,inv_group_sizes);
 } 
#else
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
#endif
}

