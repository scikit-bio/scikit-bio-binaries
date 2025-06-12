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
 *  from  skbb_accapi_impl.hpp
 * by means of
 *   generate_skbb_accapi.py
 *
 */

#ifdef SKBB_ACC_NM
// do nothing if SKBB_ACC_NM is not defined

#include <stdint.h>

namespace SKBB_ACC_NM {

    // do we have access to a GPU?
    bool acc_found_gpu();

    // is the implementation async, and need the alt structures?
    bool acc_need_alt();

    // wait for the async compute to finish
    void acc_wait();

    // create the equivalent buffer in the device memory space, if partitioned
    // the content in undefined
    template<class T>
    void acc_create_buf(T *buf, uint64_t start, uint64_t end);

    // create the equivalent buffer in the device memory space, if partitioned
    // also copy the buffer over
    template<class T>
    void acc_copyin_buf(T *buf, uint64_t start, uint64_t end);

    // make a copy from host to device buffer, if partitioned
    template<class T>
    void acc_update_device(T *buf, uint64_t start, uint64_t end);

    // make a copy from device to host buffer, if partitioned
    // destroy the equivalent buffer in the device memory space, if partitioned
    template<class T>
    void acc_copyout_buf(T *buf, uint64_t start, uint64_t end);

    // destroy the equivalent buffer in the device memory space, if partitioned
    template<class T>
    void acc_destroy_buf(T *buf, uint64_t start, uint64_t end);

}

#endif
