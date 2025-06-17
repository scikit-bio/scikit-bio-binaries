/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2025--, scikit-bio development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 * These are the public functions accessible from the
 *  libskbb
 * shared library,
 * the source of which is maintained in the scikit-bio-binaries repo.
 *  https://github.com/scikit-bio/scikit-bio-binaries
 *
 */

#ifndef SKBB_EXTERN_UTIL_H
#define SKBB_EXTERN_UTIL_H

#ifdef __cplusplus
#define EXTERN extern "C"
#else
#include <stdbool.h>
#define EXTERN
#endif

/* ====================================================== */

/*
 * Version of the API this header refers to.
 *
 * Note to developers:
 * Increase the number every time there is a change in the API.
 *
 */
#define SKBB_API_CURRENT_VERSION 1

/*
 * All other functions have a min API constant defined
 * When supporting flexible dependencies,
 * use them to determine if the a fuction
 * has its implementation in the shared library in use
 */


/* What API version is the shared library providing */
EXTERN unsigned int skbb_get_api_version();

/* ====================================================== */

/*
 * Set random seed used by any and all the functions
 * in this library
 */

#define SKBB_RANDOM_API_MIN_VERSION 1

EXTERN void skbb_set_random_seed(unsigned int new_seed);

/* ====================================================== */


/*
 * If the library supports accelerated compute
 * these functions control its logic.
 *
 * By default, the accelerator implementation is auto-detected, 
 * but the the environment variable
 *   SKBB_USE_GPU = N
 * allows for disabling of accelerated compute altogether.
 *
 * To allow monitoring of the selection logic, the env. variable
 *   SKBB_GPU_INFO = Y
 * can be used to enable runtime reporting.
 *
 * The following functions can both
 * check what auto-detection discovered
 * and force a different outcome.
 */
#define SKBB_ACC_API_MIN_VERSION 1

/* no acceleration, always possible */
#define SKBB_ACC_CPU 0
/* true acceleration, may not be available */
#define SKBB_ACC_NV  1
#define SKBB_ACC_AMD 2

/* What variant of accelerated code is currently in use */
EXTERN unsigned int skbb_get_acc_mode();

/*
 *  Override the acceleration mode.
 *  Use with caution.
 *
 *  Returns:
 *    true if acc_type should be one of the above supported ones
 *    false otherwise, and the internal value is not changed
 */
EXTERN bool skbb_set_acc_mode(unsigned int acc_type);

/* Is reporting of accelerated compute use enabled? */
EXTERN bool skbb_is_acc_reporting();

/* Change accelerated compute use reporting flag */
EXTERN void skbb_set_acc_reporting_flag(bool new_value);

#endif
