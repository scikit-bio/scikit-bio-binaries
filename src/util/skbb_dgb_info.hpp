/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2025--, scikit-bio development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#ifndef SKBB_DBG_INFO_HPP
#define SKBB_DBG_INFO_HPP

#include <chrono>

// To be used once per function
#define SETUP_TDBG(method) const char *tdbg_method=method; \
                          bool print_tdbg = false;\
                          if (const char* env_p = std::getenv("SKBB_TIMING_INFO")) { \
                            print_tdbg = true; \
                            std::string env_s(env_p); \
                            if ((env_s=="NO") || (env_s=="N") || (env_s=="no") || (env_s=="n") || \
                               (env_s=="NEVER") || (env_s=="never")) print_tdbg = false; \
                          } \
                          std::chrono::time_point<std::chrono::high_resolution_clock> tgdb_t0; \
                          if(print_tdbg) {\
                            tgdb_t0 = std::chrono::high_resolution_clock::now(); \
                            printf("INFO (skbio_bins): Starting %s\n",tdbg_method); \
                          }

// Can be invoked many time in a function
// But only after SETUP_TDBG was called once
#define TDBG_STEP(sname) if(print_tdbg) {\
                           auto tgdb_t1 = std::chrono::high_resolution_clock::now(); \
                           auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(tgdb_t1 - tgdb_t0); \
                           printf("INFO (skbio_bins): dt %6.2f : Completed %s.%s\n",time_span.count(),tdbg_method,sname); \
                           tgdb_t0 = tgdb_t1; \
                         }

#endif

