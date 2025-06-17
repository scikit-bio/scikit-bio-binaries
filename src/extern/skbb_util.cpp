/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2025--, scikit-bio development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#include "extern/util.h"

#include "util/rand.hpp"
#include "util/skbb_detect_acc.hpp"

void skbb_set_random_seed(unsigned int new_seed) {
  skbb::set_random_seed(new_seed);
}

unsigned int skbb_get_acc_mode() {
  return skbb::check_use_acc();
}

bool skbb_set_acc_mode(unsigned int acc_type) {
  return skbb::set_use_acc(acc_type);
}

bool skbb_is_acc_reporting() {
  return skbb::check_report_acc();
}

void skbb_set_acc_reporting_flag(bool new_value) {
  skbb::set_report_acc(new_value);
}

