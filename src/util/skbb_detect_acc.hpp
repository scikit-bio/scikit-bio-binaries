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
 * Classes, methods and functions that accelerator selection logic
 */

#ifndef SKBB_DETECT_ACC_HPP
#define SKBB_DETECT_ACC_HPP

#include <stdint.h>

namespace skbb {
  static constexpr unsigned int ACC_CPU=0;
  static constexpr unsigned int ACC_NV=1;
  static constexpr unsigned int ACC_AMD=2;

  // is reporting of ACC use enabled?
  bool check_report_acc();

  // override the default value
  void set_report_acc(bool new_value);

  // =======

  // what variant of accelerated code should I use
  unsigned int check_use_acc();

  // override the default value
  // Returns:
  //   true if acc_type should be one of the above supported ones
  //   false otherwise, and the internal value is not changed
  bool set_use_acc(unsigned int acc_type);

#ifdef SKBB_ENABLE_CPU_X86_LEVELS
  // optional, internal, x86 levels interface
  static constexpr unsigned int CPU_X86_BASE=0;
  static constexpr unsigned int CPU_X86_V3=3;
  static constexpr unsigned int CPU_X86_V4=4;

  // is reporting of ACC use enabled?
  bool check_report_cpu_x86();

  // override the default value
  void set_report_cpu_x86(bool new_value);

  // =======

  // what variant of accelerated code should I use
  unsigned int check_use_cpu_x86();

  // override the default value
  // Returns:
  //   true if acc_type should be one of the above supported ones
  //   false otherwise, and the internal value is not changed
  bool set_use_cpu_x86(unsigned int cpu_x86_level);

#endif /* SKBB_ENABLE_CPU_X86_LEVELS */
}

#endif
