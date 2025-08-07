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

#if defined(SKBB_CUDA)

#include <cuda.h>
#include <cuda_runtime.h>
#include <stdexcept>

#elif defined(SKBB_HIP)

#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>
#include <stdexcept>

#elif defined(_OPENACC)

#include <openacc.h>

#elif !(defined(_OPENACC) || defined(OMPGPU))

#include <omp.h>

#define SKBB_CPU Y

#endif

static inline int pmn_get_max_parallelism_T() {
#if defined(SKBB_CPU)
  // No good reason to do more than max threads
  // (but use 2x to reduce thread spawning overhead)
  // but we do use 16x blocking, so account for that, too
  return 2*omp_get_max_threads()*16;

#elif defined(SKBB_CUDA)
  int deviceID;
  cudaDeviceProp props;

  if (cudaGetDevice(&deviceID)!=cudaSuccess) return 4000; // should never get in here, but just in case
  cudaGetDeviceProperties(&props, deviceID);

  // GPUs typically need at least 64 blocks per SM to be fully loaded
  // We want a few multiples of that to deal with unbalanced load
  // Most permanovas are multiple of 100, so double that makes a good constant
  return 200*props.multiProcessorCount;

#elif defined(SKBB_HIP)
  int deviceID = 0;
  hipDeviceProp_t props;

  if (hipGetDevice(&deviceID)!=hipSuccess) return 4000; // should never get in here, but just in case
  if (hipGetDeviceProperties(&props, deviceID)!=hipSuccess) throw std::runtime_error("hipGetDeviceProperties failed");

  // GPUs typically need at least 64 blocks per SM to be fully loaded
  // We want a few multiples of that to deal with unbalanced load
  // Most permanovas are multiple of 100, so double that makes a good constant
  return 200*props.multiProcessorCount;

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
template<class TFloat, int NBLOCK>
static inline void pmn_f_stat_sW_block(
		const uint32_t n_dims,
		const TFloat * mat,
		const uint32_t *grouping_arr[],
		const TFloat *inv_group_sizes,
		TFloat out_sW_arr[]) {
  // Use full precision for intermediate compute, to minimize accumulation errors
  double s_W[NBLOCK];
  for (int i=0; i<NBLOCK; i++) s_W[i] = 0.0;

  constexpr uint32_t TILE = 128;  // 128 is big enough to speed up the code, without cache trashing

  for (uint32_t trow=0; trow < (n_dims-1); trow+=TILE) {   // no columns in last row
    for (uint32_t tcol=trow+1; tcol < n_dims; tcol+=TILE) { // diagonal is always zero
      const uint32_t max_row = std::min(trow+TILE,n_dims-1);
      const uint32_t max_col = std::min(tcol+TILE,n_dims);

      for (uint32_t row=trow; row < max_row; row++) {
        const uint32_t min_col = std::max(tcol,row+1);
        uint32_t group_idx[NBLOCK];
	for (int i=0; i<NBLOCK; i++) group_idx[i] = grouping_arr[i][row];

        // Use full precision for intermediate compute, to minimize accumulation errors
        double local_s_W[NBLOCK];
        for (int i=0; i<NBLOCK; i++) local_s_W[i] = 0.0;
        const TFloat * mat_row = mat + uint64_t(row)*uint64_t(n_dims);
        for (uint32_t col=min_col; col < max_col; col++) {
	    // speculatively read, we will likely use it at least in one of the ifs
	    // small penalty if we mis-predicted, no effect on semantics
            TFloat val = mat_row[col];  // mat[row,col];
	    val = val*val;
	    for (int i=0; i<NBLOCK; i++) if (grouping_arr[i][col] == group_idx[i]) local_s_W[i] += val;
        }
	for (int i=0; i<NBLOCK; i++) s_W[i] += local_s_W[i]*inv_group_sizes[group_idx[i]];
      }
    }
  }

  for (int i=0; i<NBLOCK; i++) out_sW_arr[i] = s_W[i];
}
#endif

// Compute PERMANOVA pseudo-F partial statistic
// mat is symmetric matrix of size n_dims x n_dims
// groupings is a matrix of size n_dims x n_grouping_dims
// inv_group_sizes is an array of size maxel(groupings)
// Results in group_sWs, and array of size n_grouping_dims

#if defined(SKBB_CPU)

template<class TFloat>
static inline void pmn_f_stat_sW_cpu(
		const uint32_t n_dims,
		const TFloat * mat,
		const uint32_t n_grouping_dims,
		const uint32_t *groupings,
		const TFloat *inv_group_sizes,
		TFloat *group_sWs) {
 constexpr int NBLOCK = 16;
// CPU version, call function
#pragma omp parallel for
 for (uint32_t gblock=0; gblock < n_grouping_dims; gblock+=NBLOCK) {
    const uint32_t *grouping_arr[NBLOCK];
    for (int i=0; i<NBLOCK; i++) grouping_arr[i] = groupings + uint64_t(gblock+i)*uint64_t(n_dims);
    if ((gblock+(NBLOCK-1))<n_grouping_dims) {
      // can process multiple permutations per pass
      const uint32_t grouping_el = gblock;
      pmn_f_stat_sW_block<TFloat,NBLOCK>(n_dims,mat,
		    grouping_arr,
		    inv_group_sizes,
		    group_sWs+grouping_el);
    } else {
      uint32_t gblock2=gblock;
      // Note: if we ever change NBLOCK, we need to update this logic, too
      if ((gblock2+(8-1))<n_grouping_dims) {
        // can process multiple permutations per pass
        const uint32_t grouping_el = gblock2;
        pmn_f_stat_sW_block<TFloat,8>(n_dims,mat,
		    grouping_arr,
		    inv_group_sizes,
		    group_sWs+grouping_el);
        gblock2+=8;
      }
      if ((gblock2+(4-1))<n_grouping_dims) {
        // can process multiple permutations per pass
        const uint32_t grouping_el = gblock2;
        pmn_f_stat_sW_block<TFloat,4>(n_dims,mat,
		    grouping_arr,
		    inv_group_sizes,
		    group_sWs+grouping_el);
        gblock2+=4;
      }
      if ((gblock2+(2-1))<n_grouping_dims) {
        // can process multiple permutations per pass
        const uint32_t grouping_el = gblock2;
        pmn_f_stat_sW_block<TFloat,2>(n_dims,mat,
		    grouping_arr,
		    inv_group_sizes,
		    group_sWs+grouping_el);
        gblock2+=2;
      }
      if (gblock2<n_grouping_dims) {
        // can process multiple permutations per pass
        const uint32_t grouping_el = gblock2;
        pmn_f_stat_sW_block<TFloat,1>(n_dims,mat,
		    grouping_arr,
		    inv_group_sizes,
		    group_sWs+grouping_el);
      }
    } // if can use big block
 } 
}

#elif (defined(SKBB_CUDA) || defined(SKBB_HIP))

template<class TFloat>
__global__ void pmn_f_stat_sW_cuda_one(
		const uint32_t n_dims,
		const TFloat * mat,
		const uint32_t n_grouping_dims,
		const uint32_t *groupings,
		const TFloat *inv_group_sizes,
		TFloat *group_sWs) {
    const uint32_t grouping_el = blockIdx.x;
    const uint32_t icol = threadIdx.x;
    const uint32_t TILE = blockDim.x;
    __shared__ double s_W_arr[128]; // we will not have more than 128 threads
    __shared__ uint32_t row_grouping[128]; // we will not have more than 128 threads

    const uint32_t *grouping = groupings + uint64_t(grouping_el)*uint64_t(n_dims);
    // Use full precision for intermediate compute, to minimize accumulation errors
    double s_W = 0.0;

    for (uint32_t trow=0; trow < (n_dims-1); trow+=TILE) {   // no columns in last row
     for (uint32_t tcol=trow+1; tcol < n_dims; tcol+=TILE) { // diagonal is always zero
      const uint32_t max_row = min(trow+TILE,n_dims-1);
      const uint32_t max_col = min(tcol+TILE,n_dims);

      __syncthreads();
      // read grouping[row] in advance into shared memory, using the whole block
      // helps with read coalesence here, and shared memory is OK to read sparse
      {
	 const uint32_t irow = threadIdx.x;
         uint32_t row=trow+irow;
	 if (row < max_row) row_grouping[irow] = grouping[row];
      }
      __syncthreads();

      TFloat local_s_W = 0.0;
      for (uint32_t row=trow; row < max_row; row++) {
        const uint32_t min_col = max(tcol,row+1);
	const uint32_t irow = row-trow;
        uint32_t group_idx = row_grouping[irow];

        const TFloat * mat_row = mat + uint64_t(row)*uint64_t(n_dims);
	uint32_t col= min_col + icol;
        if ( col < max_col) { // since TILE==col_stride, we get at most one thread picking this up
            if (grouping[col] == group_idx) {
                TFloat val = mat_row[col];  // mat[row,col];
                local_s_W += val * val * inv_group_sizes[group_idx];
            }
        }
      }
      s_W += local_s_W;
     } // for tcol
    } // for trow

    // now we need to collect data from all the threads
    s_W_arr[icol] = s_W;
    __syncthreads();
    // first sum across workers (4*32)
    if (threadIdx.x<32) {
      s_W = 0.0;
      for (uint32_t t=threadIdx.x; t < blockDim.x; t+=32) {
        s_W += s_W_arr[t]; 
      }
      s_W_arr[threadIdx.x] = s_W;
    }
    __syncthreads();
    // then sum using partial warps (4*8)
    if (threadIdx.x<8) {
      s_W = 0.0;
      for (uint32_t t=threadIdx.x; t < 32; t+=8) {
        s_W += s_W_arr[t]; 
      }
      s_W_arr[threadIdx.x] = s_W;
    }
    __syncthreads();
    // finally sum the remainder and write to global memory
    if (threadIdx.x==0) {
      s_W = 0.0;
      for (uint32_t t=threadIdx.x; t < 8; t++) {
        s_W += s_W_arr[t]; 
      }
      group_sWs[grouping_el] = s_W;
    }
}

template<class TFloat>
static inline void pmn_f_stat_sW_gpu(
		const uint32_t n_dims,
		const TFloat * mat,
		const uint32_t n_grouping_dims,
		const uint32_t *groupings,
		const TFloat *inv_group_sizes,
		TFloat *group_sWs) {
  pmn_f_stat_sW_cuda_one<TFloat><<<n_grouping_dims,128>>>(n_dims,mat,n_grouping_dims,groupings,inv_group_sizes,group_sWs);
#if defined(SKBB_CUDA)
  if (cudaDeviceSynchronize()!=cudaSuccess) throw std::runtime_error("cudaDeviceSynchronize failed");
#else
  if (hipDeviceSynchronize()!=hipSuccess) throw std::runtime_error("hipDeviceSynchronize failed");
#endif
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
#if defined(SKBB_CPU)
  pmn_f_stat_sW_cpu(n_dims, mat, n_grouping_dims, groupings, inv_group_sizes, group_sWs);
#else
  pmn_f_stat_sW_gpu(n_dims, mat, n_grouping_dims, groupings, inv_group_sizes, group_sWs);
#endif
}
