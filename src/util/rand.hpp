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
    RandomGeneratorArray(uint32_t array_size);
    ~RandomGeneratorArray() {
      delete[] randomGenerators;
    }

    // forbid copying the object
    RandomGeneratorArray(const RandomGeneratorArray &other) = delete;
    RandomGeneratorArray& operator=(const RandomGeneratorArray &other) = delete;

    // return the index-specific random generator
    TRAND& get_random_generator(uint32_t i) {
      return randomGenerators[i];
    }
  };
}

#endif
