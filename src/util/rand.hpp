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
 * Classes, methods and functions that provide random number handling
 */

#ifndef SKBB_RAND_HPP
#define SKBB_RAND_HPP

#include <stdint.h>
#include <random>

namespace skbb {

  typedef std::mt19937 TRAND;

  // Set random seed used by any and all the functions
  // in this library
  void set_random_seed(uint32_t new_seed);

  // return the global random generator
  TRAND& get_random_generator();

  // List of random generators, e.g. for use in threads
  class RandomGeneratorArray {
  private:
    TRAND *randomGenerators;
  public:
    // Arguments:
    //   array_size - number of random generators
    //   seed       - Optional random seed, if non-negative. Use system random seed if <0
    RandomGeneratorArray(uint32_t array_size, int seed);
    RandomGeneratorArray(uint32_t array_size); // equivalent to passing seed<0
    ~RandomGeneratorArray() {
      delete[] randomGenerators;
    }

    // forbid copying the object
    RandomGeneratorArray(const RandomGeneratorArray &other) = delete;
    RandomGeneratorArray& operator=(const RandomGeneratorArray &other) = delete;

    // return the index-specific random generatora
    // Note: No boundary checking, caller must ensure i is in a valid range
    TRAND& get_random_generator(uint32_t i) {
      return randomGenerators[i];
    }
  };
}

#endif
