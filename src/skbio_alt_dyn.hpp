/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2025, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 *
 * This file is used to create the necessary interfaces between
 *   unifrac_task.hpp and unifrac_task_impl.hpp
 * by means of
 *   generate_unifrac_task_noclass.py
 *
 */

#include <stdint.h>
#include "unifrac_accapi.hpp"

#ifndef SUCMP_NM
/* create a default */
#define SUCMP_NM su_cpu
#endif

namespace SUCMP_NM {

	// Return the recommended max parallelism to use in pmn_f_stat_sW
	int pmn_get_max_parallelism();

	// Compute PERMANOVA pseudo-F partial statistic
	// mat is symmetric matrix of size n_dims x n_dims
	// groupings is a matrix of size n_dims x n_grouping_dims
	// inv_group_sizes is an array of size maxel(groupings)
	// Results in group_sWs, and array of size n_grouping_dims
	// Note: Best results when n_grouping_dims fits in L1 cache
	template<class TFloat>
	void pmn_f_stat_sW(
			const TFloat * mat,
			const uint32_t n_dims,
		        const uint32_t *groupings,
			const uint32_t n_grouping_dims,
			const TFloat *inv_group_sizes,
			TFloat *group_sWs);
}

