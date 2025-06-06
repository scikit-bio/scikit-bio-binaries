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

#if defined(OMPGPU)

#include <omp.h>

#elif defined(_OPENACC)

#include <openacc.h>

#endif

static inline bool acc_found_gpu_T() {
#if defined(OMPGPU)
  return omp_get_num_devices() > 0;
#elif defined(_OPENACC)
  return acc_get_device_type() != acc_device_host;
#else
  return false;
#endif
}

// is the implementation async, and need the alt structures?
static inline bool acc_need_alt_T() {
#if defined(_OPENACC) || defined(OMPGPU)
   return true;
#else
   return false;
#endif
}

static inline void acc_wait_T() {
#if defined(OMPGPU)
    // TODO: Change if we ever implement async in OMPGPU
#elif defined(_OPENACC)
#pragma acc wait
#endif
}

template<class TNum>
static inline void acc_create_buf_T(
		TNum *buf,
		uint64_t start, uint64_t end) {
#if defined(OMPGPU)
#pragma omp target enter data map(alloc:buf[start:end])
#elif defined(_OPENACC)
#pragma acc enter data create(buf[start:end])
#endif
}

template<class TNum>
static inline void acc_copyin_buf_T(
		TNum *buf,
		uint64_t start, uint64_t end) {
#if defined(OMPGPU)
#pragma omp target enter data map(to:buf[start:end])
#elif defined(_OPENACC)
#pragma acc enter data copyin(buf[start:end])
#endif    
}

template<class TNum>
static inline void acc_update_device_T(
		TNum *buf,
		uint64_t start, uint64_t end) {
#if defined(OMPGPU)
#pragma omp target update to(buf[start:end])
#elif defined(_OPENACC)
#pragma acc update device(buf[start:end])
#endif
}

template<class TNum>
static inline void acc_copyout_buf_T(
		TNum *buf,
		uint64_t start, uint64_t end) {
#if defined(OMPGPU)
#pragma omp target exit data map(from:buf[start:end])
#elif defined(_OPENACC)
#pragma acc exit data copyout(buf[start:end])
#endif
}

template<class TNum>
static inline void acc_destroy_buf_T(
		TNum *buf,
		uint64_t start, uint64_t end) {
#if defined(OMPGPU)
#pragma omp target exit data map(delete:buf[start:end])
#elif defined(_OPENACC)
#pragma acc exit data delete(buf[start:end])
#endif
}

