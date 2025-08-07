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

#if defined(SKBB_CUDA)

#include <cuda_runtime_api.h>
#include <stdexcept>

#elif defined(SKBB_HIP)

#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>
#include <stdexcept>

#elif defined(OMPGPU)

#include <omp.h>

#elif defined(_OPENACC)

#include <openacc.h>

#endif

static inline bool acc_found_gpu_T() {
#if defined(SKBB_CUDA)
  int deviceCount;
  cudaError_t error = cudaGetDeviceCount(&deviceCount);
  if (error != cudaSuccess) return false;
  return deviceCount != 0;
#elif defined(SKBB_HIP)
  int deviceCount;
  hipError_t error = hipGetDeviceCount(&deviceCount);
  if (error != hipSuccess) return false;
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
#if defined(_OPENACC) || defined(OMPGPU) || defined(SKBB_CUDA) || defined(SKBB_HIP)
   return true;
#else
   return false;
#endif
}

static inline void acc_wait_T() {
#if defined(SKBB_CUDA)
  if (cudaDeviceSynchronize()!=cudaSuccess) throw std::runtime_error("cudaDeviceSynchronize failed");
#elif defined(SKBB_HIP)
  if (hipDeviceSynchronize()!=hipSuccess) throw std::runtime_error("hipDeviceSynchronize failed");
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
#if defined(SKBB_CUDA)
  if (cudaMalloc((void**)buf_device, sizeof(TNum) * size)!=cudaSuccess) throw std::runtime_error("cudaMalloc failed");
#elif defined(SKBB_HIP)
  if (hipMalloc((void**)buf_device, sizeof(TNum) * size)!=hipSuccess) throw std::runtime_error("hipMalloc failed");
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
#if defined(SKBB_CUDA)
  if (cudaMalloc((void**)buf_device, sizeof(TNum) * size)!=cudaSuccess) throw std::runtime_error("cudaMalloc failed");
  if (cudaMemcpy(*buf_device, buf_host, sizeof(TNum) * size, cudaMemcpyHostToDevice)!=cudaSuccess) throw std::runtime_error("cudaMemcpy failed");
#elif defined(SKBB_HIP)
  if (hipMalloc((void**)buf_device, sizeof(TNum) * size)!=hipSuccess) throw std::runtime_error("hipMalloc failed");
  if (hipMemcpy(*buf_device, buf_host, sizeof(TNum) * size, hipMemcpyHostToDevice)!=hipSuccess) throw std::runtime_error("hipMemcpy failed");
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
#if defined(SKBB_CUDA)
  if (cudaMemcpy(buf_device+start, buf_host+start, sizeof(TNum) * (end-start), cudaMemcpyHostToDevice)!=cudaSuccess) throw std::runtime_error("cudaMemcpy failed");
#elif defined(SKBB_HIP)
  if (hipMemcpy(buf_device+start, buf_host+start, sizeof(TNum) * (end-start), hipMemcpyHostToDevice)!=hipSuccess) throw std::runtime_error("hipMemcpy failed");
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
#if defined(SKBB_CUDA)
  if (cudaMemcpy(buf_host, buf_device, sizeof(TNum) * size, cudaMemcpyDeviceToHost)!=cudaSuccess) throw std::runtime_error("cudaMemcpy failed");
  cudaFree(buf_device);
#elif defined(SKBB_HIP)
  if (hipMemcpy(buf_host, buf_device, sizeof(TNum) * size, hipMemcpyDeviceToHost)!=hipSuccess) throw std::runtime_error("hipMemcpy failed");
  if (hipFree(buf_device)!=hipSuccess) {} // ignore any errors, not critical
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
#if defined(SKBB_CUDA)
  cudaFree(buf_device);
#elif defined(SKBB_HIP)
  if (hipFree(buf_device)!=hipSuccess) {} // ignore any errors, not critical
#elif defined(OMPGPU)
#pragma omp target exit data map(delete:buf_device[0:size])
#elif defined(_OPENACC)
#pragma acc exit data delete(buf_device[0:size])
#endif
}

