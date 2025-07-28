/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2025, UniFrac development team.
 * Copyright (c) 2025--, scikit-bio development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

// Implementation specific acceleration primitives

/*
 *
 * This file is used to create the necessary interfaces
 * by means of
 *   generate_skbb_accapi.py
 *
 * Anything ending in _T will get a acc-specific function wrapper.
 *
 */

#include "util/skbb_accapi.hpp"
#include <cstdlib>

#if defined(CUDA)

#include <cuda_runtime_api.h>

#elif defined(OMPGPU)

#include <omp.h>

#elif defined(_OPENACC)

#include <openacc.h>

#endif

static inline bool acc_found_gpu_T() {
#if defined(CUDA)
  int deviceCount;
  cudaError_t error = cudaGetDeviceCount(&deviceCount);
  if (error != cudaSuccess) return false;
  return deviceCount != 0;
#elif defined(OMPGPU)
  return omp_get_num_devices() > 0;
#elif defined(_OPENACC)
  return acc_get_device_type() != acc_device_host;
#else
  return false;
#endif
}

// is the implementation async, and need the alt structures?
static inline bool acc_need_alt_T() {
#if defined(_OPENACC) || defined(OMPGPU) || defined(CUDA)
   return true;
#else
   return false;
#endif
}

static inline void acc_wait_T() {
#if defined(CUDA)
  cudaDeviceSynchronize();
#elif defined(OMPGPU)
    // TODO: Change if we ever implement async in OMPGPU
#elif defined(_OPENACC)
#pragma acc wait
#endif
}

template<class TNum>
static inline void acc_create_buf_T(
		TNum*  buf_host,
		TNum** buf_device,
		uint64_t size) {
#if defined(CUDA)
  cudaMalloc((void**)buf_device, sizeof(TNum) * size);
#elif defined(OMPGPU)
#pragma omp target enter data map(alloc:buf_host[0:size])
  *buf_device = buf_host;
#elif defined(_OPENACC)
#pragma acc enter data create(buf_host[0:size])
  *buf_device = buf_host;
#else
  *buf_device = buf_host;
#endif
}

template<class TNum>
static inline void acc_copyin_buf_T(
		TNum*  buf_host,
		TNum** buf_device,
		uint64_t size) {
#if defined(CUDA)
  cudaMalloc((void**)buf_device, sizeof(TNum) * size);
  cudaMemcpy(*buf_device, buf_host, sizeof(TNum) * size, cudaMemcpyHostToDevice);
#elif defined(OMPGPU)
#pragma omp target enter data map(to:buf_host[0:size])
  *buf_device = buf_host;
#elif defined(_OPENACC)
#pragma acc enter data copyin(buf_host[0:size])
  *buf_device = buf_host;
#else
  *buf_device = buf_host;
#endif    
}

template<class TNum>
static inline void acc_update_device_T(
		TNum *buf_host,
		TNum *buf_device,
		uint64_t start, uint64_t end) {
#if defined(CUDA)
  cudaMemcpy(buf_device+start, buf_host+start, sizeof(TNum) * (end-start), cudaMemcpyHostToDevice);
#elif defined(OMPGPU)
 // assert buf_host==buf_device
#pragma omp target update to(buf_host[start:end])
#elif defined(_OPENACC)
 // assert buf_host==buf_device
#pragma acc update device(buf_host[start:end])
#endif
}

template<class TNum>
static inline void acc_copyout_buf_T(
		TNum *buf_host,
		TNum *buf_device,
		uint64_t size) {
#if defined(CUDA)
  cudaMemcpy(buf_host, buf_device, sizeof(TNum) * size, cudaMemcpyDeviceToHost);
  cudaFree(buf_device);
#elif defined(OMPGPU)
 // assert buf_host==buf_device
#pragma omp target exit data map(from:buf_device[0:size])
#elif defined(_OPENACC)
 // assert buf_host==buf_device
#pragma acc exit data copyout(buf_device[0:size])
#endif
}

template<class TNum>
static inline void acc_destroy_buf_T(
		TNum *buf_device,
		uint64_t size) {
#if defined(CUDA)
  cudaFree(buf_device);
#elif defined(OMPGPU)
#pragma omp target exit data map(delete:buf_device[0:size])
#elif defined(_OPENACC)
#pragma acc exit data delete(buf_device[0:size])
#endif
}

