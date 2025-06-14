/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2025--, scikit-bio development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
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
 * Set random seed used by any and all the functions
 * in this library
 */

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

extern unsigned int SKBB_ACC_CPU;  /* no acceleration, always possible */
extern unsigned int SKBB_ACC_NV;
extern unsigned int SKBB_ACC_AMD;

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
