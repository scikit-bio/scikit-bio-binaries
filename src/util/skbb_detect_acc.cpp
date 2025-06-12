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

#include "util/skbb_detect_acc.hpp"

#ifdef SKBB_ENABLE_ACC_NV
#define SKBB_ACC_NM  skbb_acc_nv
#include "util/skbb_accapi.hpp"
#undef SKBB_ACC_NM
#endif

#ifdef SKBB_ENABLE_ACC_AMD
#define SKBB_ACC_NM  skbb_acc_amd
#include "util/skbb_accapi.hpp"
#undef SKBB_ACC_NM
#endif

#include <stdlib.h> 
#include <string> 

static constexpr unsigned int ACC_REPORT_FALSE=0;
static constexpr unsigned int ACC_REPORT_TRUE=1;

// test only once, then use persistent value
static int skbb_report_acc = -1;
static int skbb_use_acc = -1;

// test only once, then use persistent value
bool skbb::check_report_acc() {
 if (skbb_report_acc>=0) return skbb_report_acc!=ACC_REPORT_FALSE; // keep the cached version

 bool print_info = false;

 if (const char* env_p = std::getenv("SKBB_GPU_INFO")) {
   print_info = true;
   std::string env_s(env_p);
   if ((env_s=="NO") || (env_s=="N") || (env_s=="no") || (env_s=="n") ||
       (env_s=="NEVER") || (env_s=="never")) {
     print_info = false;
   }
 }

 skbb_report_acc = print_info ? ACC_REPORT_TRUE : ACC_REPORT_FALSE;

 return print_info;
}

void skbb::set_report_acc(bool new_value) {
 skbb_report_acc = new_value ? ACC_REPORT_TRUE : ACC_REPORT_FALSE;
}


unsigned int skbb::check_use_acc() {
 if (skbb_use_acc>=0) return skbb_use_acc; // keep the cached version

 bool print_info = skbb::check_report_acc();

 int detected_acc = skbb::ACC_CPU;
#if defined(SKBB_ENABLE_ACC_NV)
 bool detected_nv_acc = skbb_acc_nv::acc_found_gpu();
 if (print_info) {
   if (detected_nv_acc) {
     printf("INFO (skbio_bins): NVIDIA GPU detected\n");
   } else {
     printf("INFO (skbio_bins): NVIDIA GPU not detected\n");
   }
 }
 if ((detected_acc==skbb::ACC_CPU) && detected_nv_acc) {
   detected_acc = skbb::ACC_NV;
   if (const char* env_p = std::getenv("SKBB_USE_NVIDIA_GPU")) {
     std::string env_s(env_p);
     if ((env_s=="NO") || (env_s=="N") || (env_s=="no") || (env_s=="n") ||
         (env_s=="NEVER") || (env_s=="never")) {
       if (print_info) printf("INFO (skbio_bins): NVIDIA GPU was detected but use explicitly disabled\n");
       detected_acc = skbb::ACC_CPU;
     }
   }
 }
#endif

#if defined(SKBB_ENABLE_ACC_AMD)
 bool detected_amd_acc = skbb_acc_amd::acc_found_gpu();
 if (print_info) {
   if (detected_amd_acc) {
     printf("INFO (skbio_bins): AMD GPU detected\n");
   } else {
     printf("INFO (skbio_bins): AMD GPU not detected\n");
   }
 }
 if ((detected_acc==skbb::ACC_CPU) && detected_amd_acc) {
   detected_acc = skbb::ACC_AMD;
   if (const char* env_p = std::getenv("SKBB_USE_AMD_GPU")) {
     std::string env_s(env_p);
     if ((env_s=="NO") || (env_s=="N") || (env_s=="no") || (env_s=="n") ||
         (env_s=="NEVER") || (env_s=="never")) {
       if (print_info) printf("INFO (skbio_bins): AMD GPU was detected but use explicitly disabled\n");
       detected_acc = skbb::ACC_CPU;
     }
   }
 }
#endif

 if (const char* env_p = std::getenv("SKBB_USE_GPU")) {
   std::string env_s(env_p);
   if ((env_s=="NO") || (env_s=="N") || (env_s=="no") || (env_s=="n") ||
       (env_s=="NEVER") || (env_s=="never")) {
     if (detected_acc!=skbb::ACC_CPU) {
        if (print_info) printf("INFO (skbio_bins): GPU was detected but use explicitly disabled\n");
         detected_acc = skbb::ACC_CPU;
     }
   }
 }

 if (print_info) {
   if (detected_acc == skbb::ACC_CPU) {
     printf("INFO (skbio_bins): Using CPU (not GPU)\n");
#if defined(SKBB_ENABLE_ACC_NV)
   } else if (detected_acc == skbb::ACC_NV) {
     printf("INFO (skbio_bins): Using NVIDIA GPU\n");
#endif
#if defined(SKBB_ENABLE_ACC_AMD)
   } else if (detected_acc == skbb::ACC_AMD) {
     printf("INFO (skbio_bins): Using AMD GPU\n");
#endif
   } else {
     printf("WARNING (skbio_bins): Logic error in picking GPU\n");
   }
 }
 // we can assume int is atomic
 skbb_use_acc = detected_acc;

 return skbb_use_acc;
}

// acc_type must be one of the supported ones
bool skbb::set_use_acc(unsigned int acc_type) {
 if (acc_type == skbb::ACC_CPU) {
   skbb_use_acc = acc_type; 
#if defined(SKBB_ENABLE_ACC_NV)
 } else if (acc_type == skbb::ACC_NV) {
   skbb_use_acc = acc_type; 
#endif
#if defined(SKBB_ENABLE_ACC_AMD)
 } else if (acc_type == skbb::ACC_AMD) {
   skbb_use_acc = acc_type; 
#endif
 } else {
   return false;
 }

 bool print_info = skbb::check_report_acc();
 if (print_info) {
   if (skbb_use_acc == skbb::ACC_CPU) {
     printf("INFO (skbio_bins): Using CPU (not GPU)\n");
#if defined(SKBB_ENABLE_ACC_NV)
   } else if (skbb_use_acc == skbb::ACC_NV) {
     printf("INFO (skbio_bins): Using NVIDIA GPU\n");
#endif
#if defined(SKBB_ENABLE_ACC_AMD)
   } else if (skbb_use_acc == skbb::ACC_AMD) {
     printf("INFO (skbio_bins): Using AMD GPU\n");
#endif
   } else {
     printf("WARNING (skbio_bins): Logic error in picking GPU\n");
   }
 }

 return true;
}

