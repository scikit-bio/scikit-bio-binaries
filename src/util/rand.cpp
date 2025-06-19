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
 * Classes, methods and functions that provide library-wide random number generator
 */

#include "rand.hpp"

// need a library-global random generator
static skbb::TRAND myRandomGenerator;

void skbb::set_random_seed(uint32_t new_seed) {
  myRandomGenerator.seed(new_seed);
}

skbb::TRAND& skbb::get_random_generator() {
  return myRandomGenerator;
}

static inline void rga_init(const uint32_t array_size, skbb::TRAND randomGenerators[], skbb::TRAND& rand) {
  // use the provided random generator to get a unique starting seed
  for (uint32_t el=0; el < array_size; el++) {
    auto new_seed = rand();
    randomGenerators[el].seed(new_seed);
  }
}


skbb::RandomGeneratorArray::RandomGeneratorArray(uint32_t array_size)
  : randomGenerators(new skbb::TRAND[array_size]) {
  rga_init(array_size,randomGenerators,myRandomGenerator);
}

skbb::RandomGeneratorArray::RandomGeneratorArray(uint32_t array_size, int seed)
  : randomGenerators(new skbb::TRAND[array_size]) {
  if (seed<0) {
    // no seed, just use global rand
    rga_init(array_size,randomGenerators,myRandomGenerator);
  } else {
    // use the global random generator to get a unique starting seed
    skbb::TRAND local_random;
    local_random.seed((unsigned int)(seed));
    // now use the loca random for the rest of the logic
    rga_init(array_size,randomGenerators,local_random);
  }
}

