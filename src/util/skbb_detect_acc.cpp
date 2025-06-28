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
#include <string.h>

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

#ifdef SKBB_ENABLE_CPU_X86_LEVELS

// test only once, then use persistent value
static int skbb_report_cpu_x86 = -1;
static int skbb_use_cpu_x86 = -1;

// test only once, then use persistent value
bool skbb::check_report_cpu_x86() {
 if (skbb_report_cpu_x86>=0) return skbb_report_cpu_x86!=ACC_REPORT_FALSE; // keep the cached version

 bool print_info = false;

 if (const char* env_p = std::getenv("SKBB_CPU_INFO")) {
   print_info = true;
   std::string env_s(env_p);
   if ((env_s=="NO") || (env_s=="N") || (env_s=="no") || (env_s=="n") ||
       (env_s=="NEVER") || (env_s=="never")) {
     print_info = false;
   }
 }

 skbb_report_cpu_x86 = print_info ? ACC_REPORT_TRUE : ACC_REPORT_FALSE;

 return print_info;
}

void skbb::set_report_cpu_x86(bool new_value) {
 skbb_report_cpu_x86 = new_value ? ACC_REPORT_TRUE : ACC_REPORT_FALSE;
}


unsigned int skbb::check_use_cpu_x86() {
 if (skbb_use_cpu_x86>=0) return skbb_use_cpu_x86; // keep the cached version

 bool print_info = skbb::check_report_cpu_x86();

 int detected_cpu_x86 = skbb::CPU_X86_BASE;

 __builtin_cpu_init ();

 const char* env_max_cpu = getenv("SKBB_MAX_CPU");

 bool allow_v3 = true;
 bool allow_v4 = true;
 if (env_max_cpu!=NULL) {
    if ((strcmp(env_max_cpu,"basic")==0) ||
        (strcmp(env_max_cpu,"base")==0) ||
        (strcmp(env_max_cpu,"x86-64-v2")==0) ||
        (strcmp(env_max_cpu,"sse")==0) ||
        (strcmp(env_max_cpu,"sse3")==0) ||
        (strcmp(env_max_cpu,"sse4")==0) ||
        (strcmp(env_max_cpu,"avx")==0) ) {
      allow_v3 = false;
      allow_v4 = false;
    } else if ((strcmp(env_max_cpu,"x86-64-v3")==0) || (strcmp(env_max_cpu,"avx2")==0)) {
      allow_v4 = false;
    }
    // ignore any unknown strings
 }

 // do the check in priority order, since they are supersets
#if defined(SKBB_ENABLE_CPU_X86V4)
 if (detected_cpu_x86 == skbb::CPU_X86_BASE) {
   bool has_v4  = __builtin_cpu_supports ("x86-64-v4");
   if (has_v4) {
     if (allow_v4) {
       if (print_info) printf("INFO (skbio_bins): x86-64-v4 CPU detected\n");
       detected_cpu_x86 = skbb::CPU_X86_V4;
     } else {
       if (print_info) printf("INFO (skbio_bins): x86-64-v4 CPU was detected but use explicitly disabled\n");
     }
   } else {
     if (print_info) printf("INFO (skbio_bins): x86-64-v4 CPU not detected\n");
   }
 }
#endif

#if defined(SKBB_ENABLE_CPU_X86V3)
 if (detected_cpu_x86 == skbb::CPU_X86_BASE) {
   bool has_v3  = __builtin_cpu_supports ("x86-64-v3");
   if (has_v3) {
     if (allow_v3) {
       if (print_info) printf("INFO (skbio_bins): x86-64-v3 CPU detected\n");
       detected_cpu_x86 = skbb::CPU_X86_V3;
     } else {
       if (print_info) printf("INFO (skbio_bins): x86-64-v3 CPU was detected but use explicitly disabled\n");
     }
   } else {
     if (print_info) printf("INFO (skbio_bins): x86-64-v3 CPU not detected\n");
   }
 }
#endif

 if (print_info) {
   if (detected_cpu_x86 == skbb::CPU_X86_BASE) {
     printf("INFO (skbio_bins): Using base x86 CPU\n");
#if defined(SKBB_ENABLE_CPU_X86V3)
   } else if (detected_cpu_x86 == skbb::CPU_X86_V3) {
     printf("INFO (skbio_bins): Using x86-64-v3 CPU\n");
#endif
#if defined(SKBB_ENABLE_CPU_X86V4)
   } else if (detected_cpu_x86 == skbb::CPU_X86_V4) {
     printf("INFO (skbio_bins): Using x86-64-v4 CPU\n");
#endif
   } else {
     printf("WARNING (skbio_bins): Logic error in picking x86 CPU level\n");
   }
 }
 // we can assume int is atomic
 skbb_use_cpu_x86 = detected_cpu_x86;

 return skbb_use_cpu_x86;
}

// cpu_x86_levele must be one of the supported ones
bool skbb::set_use_cpu_x86(unsigned int cpu_x86_level) {
 if (cpu_x86_level == skbb::CPU_X86_BASE) {
   skbb_use_cpu_x86 = cpu_x86_level;
#if defined(SKBB_ENABLE_CPU_X86V3)
 } else if (cpu_x86_level == skbb::CPU_X86_V3) {
   skbb_use_cpu_x86 = cpu_x86_level;
#endif
#if defined(SKBB_ENABLE_CPU_X86V4)
 } else if (cpu_x86_level == skbb::CPU_X86_V4) {
   skbb_use_cpu_x86 = cpu_x86_level;
#endif
 } else {
   return false;
 }

 bool print_info = skbb::check_report_cpu_x86();
 if (print_info) {
   if (skbb_use_cpu_x86 == skbb::CPU_X86_BASE) {
     printf("INFO (skbio_bins): Using base x86 CPU\n");
#if defined(SKBB_ENABLE_CPU_X86V3)
   } else if (skbb_use_cpu_x86 == skbb::CPU_X86_V3) {
     printf("INFO (skbio_bins): Using x86-64-v3 CPU\n");
#endif
#if defined(SKBB_ENABLE_CPU_X86V4)
   } else if (skbb_use_cpu_x86 == skbb::CPU_X86_V4) {
     printf("INFO (skbio_bins): Using x86-64-v4 CPU\n");
#endif
   } else {
     printf("WARNING (skbio_bins): Logic error in picking x86 CPU level\n");
   }
 }

 return true;
}

#endif /* SKBB_ENABLE_CPU_X86_LEVELS */

