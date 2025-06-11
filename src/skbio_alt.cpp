/*
 * Classes, methods and unction that provide skbio-like unctionality
 */

#include "skbio_alt.hpp"
#include "util/rand.hpp"
#include <stdlib.h> 
#include <omp.h>

#include <algorithm>

// We will always have the CPU version
#define SUCMP_NM  su_cpu
#include "skbio_alt_dyn.hpp"
#undef SUCMP_NM
static constexpr int ACC_CPU=0;

#ifdef UNIFRAC_ENABLE_ACC_NV
#define SUCMP_NM  su_acc_nv
#include "skbio_alt_dyn.hpp"
#undef SUCMP_NM
static constexpr int ACC_NV=1;
#endif

#ifdef UNIFRAC_ENABLE_ACC_AMD
#define SUCMP_NM  su_acc_amd
#include "skbio_alt_dyn.hpp"
#undef SUCMP_NM
static constexpr int ACC_AMD=2;
#endif

// test only once, then use persistent value
static int skbio_use_acc = -1;

inline void skbio_check_acc() {
 if (skbio_use_acc!=-1) return; // keep the cached version

 bool print_info = false;

 if (const char* env_p = std::getenv("UNIFRAC_GPU_INFO")) {
   print_info = true;
   std::string env_s(env_p);
   if ((env_s=="NO") || (env_s=="N") || (env_s=="no") || (env_s=="n") ||
       (env_s=="NEVER") || (env_s=="never")) {
     print_info = false;
   }
 }

 int detected_acc = ACC_CPU;
#if defined(UNIFRAC_ENABLE_ACC_NV)
 bool detected_nv_acc = su_acc_nv::acc_found_gpu();
 if (print_info) {
   if (detected_nv_acc) {
     printf("INFO (skbio_alt): NVIDIA GPU detected\n");
   } else {
     printf("INFO (skbio_alt): NVIDIA GPU not detected\n");
   }
 }
 if ((detected_acc==ACC_CPU) && detected_nv_acc) {
   detected_acc = ACC_NV;
   if (const char* env_p = std::getenv("UNIFRAC_SKBIO_USE_NVIDIA_GPU")) {
     std::string env_s(env_p);
     if ((env_s=="NO") || (env_s=="N") || (env_s=="no") || (env_s=="n") ||
         (env_s=="NEVER") || (env_s=="never")) {
       if (print_info) printf("INFO (skbio_alt): NVIDIA GPU was detected but use explicitly disabled\n");
       detected_acc = ACC_CPU;
     }
   } else if (const char* env_p = std::getenv("UNIFRAC_USE_NVIDIA_GPU")) {
     std::string env_s(env_p);
     if ((env_s=="NO") || (env_s=="N") || (env_s=="no") || (env_s=="n") ||
         (env_s=="NEVER") || (env_s=="never")) {
       if (print_info) printf("INFO (skbio_alt): NVIDIA GPU was detected but use explicitly disabled\n");
       detected_acc = ACC_CPU;
     }
   }
 }
#endif
#if defined(UNIFRAC_ENABLE_ACC_AMD)
 bool detected_amd_acc = su_acc_amd::acc_found_gpu();
 if (print_info) {
   if (detected_amd_acc) {
     printf("INFO (skbio_alt): AMD GPU detected\n");
   } else {
     printf("INFO (skbio_alt): AMD GPU not detected\n");
   }
 }
 if ((detected_acc==ACC_CPU) && detected_amd_acc) {
   detected_acc = ACC_AMD;
   if (const char* env_p = std::getenv("UNIFRAC_SKBIO_USE_AMD_GPU")) {
     std::string env_s(env_p);
     if ((env_s=="NO") || (env_s=="N") || (env_s=="no") || (env_s=="n") ||
         (env_s=="NEVER") || (env_s=="never")) {
       if (print_info) printf("INFO (skbio_alt): AMD GPU was detected but use explicitly disabled\n");
       detected_acc = ACC_CPU;
     }
   } else if (const char* env_p = std::getenv("UNIFRAC_USE_AMD_GPU")) {
     std::string env_s(env_p);
     if ((env_s=="NO") || (env_s=="N") || (env_s=="no") || (env_s=="n") ||
         (env_s=="NEVER") || (env_s=="never")) {
       if (print_info) printf("INFO (skbio_alt): AMD GPU was detected but use explicitly disabled\n");
       detected_acc = ACC_CPU;
     }
   }
 }
#endif

 if (const char* env_p = std::getenv("UNIFRAC_SKBIO_USE_GPU")) {
   std::string env_s(env_p);
   if ((env_s=="NO") || (env_s=="N") || (env_s=="no") || (env_s=="n") ||
       (env_s=="NEVER") || (env_s=="never")) {
     if (detected_acc!=ACC_CPU) {
        if (print_info) printf("INFO (skbio_alt): GPU was detected but use explicitly disabled\n");
         detected_acc = ACC_CPU;
     }
   }
 } else if (const char* env_p = std::getenv("UNIFRAC_USE_GPU")) {
   std::string env_s(env_p);
   if ((env_s=="NO") || (env_s=="N") || (env_s=="no") || (env_s=="n") ||
       (env_s=="NEVER") || (env_s=="never")) {
     if (detected_acc!=ACC_CPU) {
        if (print_info) printf("INFO (skbio_alt): GPU was detected but use explicitly disabled\n");
         detected_acc = ACC_CPU;
     }
   }
 }

 if (print_info) {
   if (detected_acc == ACC_CPU) {
     printf("INFO (skbio_alt): Using CPU (not GPU)\n");
#if defined(UNIFRAC_ENABLE_ACC_NV)
   } else if (detected_acc == ACC_NV) {
     printf("INFO (skbio_alt): Using NVIDIA GPU\n");
#endif
#if defined(UNIFRAC_ENABLE_ACC_AMD)
   } else if (detected_acc == ACC_AMD) {
     printf("INFO (skbio_alt): Using AMD GPU\n");
#endif
   }
 }
 // we can assume int is atomic
 skbio_use_acc = detected_acc;
}


//
// ======================= permanova ========================
//

// Compute PERMANOVA pseudo-F partial statistic using permutations
// mat is symmetric matrix of size n_dims x n_dims
// grouping is an array of size n_dims
// group_sizes is an array of size maxel(grouping)
//
// Results in permutted_sWs, and array of size (n_perm+1)
template<class TFloat>
inline void permanova_perm_fp_sW_T(const TFloat * mat, const uint32_t n_dims,
                                   const uint32_t *grouping, 
                                   const uint32_t *group_sizes, uint32_t n_groups,
                                   const uint32_t n_perm,
                                   TFloat *permutted_sWs) {
  const uint64_t mat_size = uint64_t(n_dims)*uint64_t(n_dims);

  // There is acc-specific logic here, initialize skbio_use_acc ASAP
  skbio_check_acc();

  uint32_t PERM_CHUNK = 1; // just a dummy default
  if (skbio_use_acc==ACC_CPU) {
    PERM_CHUNK = su_cpu::pmn_get_max_parallelism();
#if defined(UNIFRAC_ENABLE_ACC_NV)
  } else if (skbio_use_acc==ACC_NV) {
    PERM_CHUNK = su_acc_nv::pmn_get_max_parallelism();
#endif
#if defined(UNIFRAC_ENABLE_ACC_AMD)
  } else if (skbio_use_acc==ACC_AMD) {
    PERM_CHUNK = su_acc_amd::pmn_get_max_parallelism();
#endif
  }

  // need temp bufffer for bulk processing
  const uint32_t step_perms = std::min(n_perm+1,PERM_CHUNK);
  const uint64_t permutted_groupings_size = uint64_t(n_dims)*uint64_t(step_perms);
  uint32_t *permutted_groupings = new uint32_t[permutted_groupings_size];

  // first copy the original in all of the buffer rows
#pragma omp parallel for schedule(static,1) default(shared)
  for (uint32_t grouping_el=0; grouping_el < step_perms; grouping_el++) {
    uint32_t *my_grouping = permutted_groupings + uint64_t(grouping_el)*uint64_t(n_dims);
    for (uint32_t i=0; i<n_dims; i++) my_grouping[i] = grouping[i];
  }

  // will also use dedicated generators for deterministic behavior in threaded environment
  skbb::RandomGeneratorArray randomGenerators(step_perms);

  // We will use only 1/N, so pre-process
  TFloat *inv_group_sizes = new TFloat[n_groups];
  for (uint32_t i=0; i<n_groups; i++) {
    inv_group_sizes[i] = TFloat(1.0)/group_sizes[i];
  }

  // for acc implementations, make a copy of mat into the GPU memory
#if defined(UNIFRAC_ENABLE_ACC_NV)
  if (skbio_use_acc==ACC_NV) {
    // must remove const as it will indeed write to GPU memory
    su_acc_nv::acc_copyin_buf(const_cast<TFloat*>(mat),0,mat_size);
    su_acc_nv::acc_copyin_buf(inv_group_sizes,0,n_groups);
  }
#endif
#if defined(UNIFRAC_ENABLE_ACC_AMD)
  if (skbio_use_acc==ACC_AMD) {
    // must remove const as it will indeed write to GPU memory
    su_acc_amd::acc_copyin_buf(const_cast<TFloat*>(mat),0,mat_size);
    su_acc_amd::acc_copyin_buf(inv_group_sizes,0,n_groups);
  }
#endif

  // now permute and compute sWs
  for (uint32_t tp=0; tp < (n_perm+1); tp+=step_perms) {
      const uint32_t max_p = std::min(tp+step_perms,n_perm+1);
      // Split in chunks for data locality

      // first, permute the buffers
#pragma omp parallel for schedule(static,1) default(shared)
      for (uint32_t p=tp; p<max_p; p++) {
         if (p!=0) { // do not permute the first one
           const uint32_t grouping_el = p-tp;
           uint32_t *my_grouping = permutted_groupings + uint64_t(grouping_el)*uint64_t(n_dims);
           std::shuffle(my_grouping, my_grouping+n_dims, randomGenerators.get_random_generator(grouping_el));
         }
      }
      // now call the actual permanova
      if (skbio_use_acc==ACC_CPU) {
        su_cpu::pmn_f_stat_sW<TFloat>(mat, n_dims,
                                      permutted_groupings, max_p-tp,
                                      inv_group_sizes,
                                      permutted_sWs+tp);
#if defined(UNIFRAC_ENABLE_ACC_NV)
      } else if (skbio_use_acc==ACC_NV) {
        su_acc_nv::pmn_f_stat_sW<TFloat>(mat, n_dims,
                                         permutted_groupings, max_p-tp,
                                         inv_group_sizes,
                                         permutted_sWs+tp);
#endif
#if defined(UNIFRAC_ENABLE_ACC_AMD)
      } else if (skbio_use_acc==ACC_AMD) {
        su_acc_amd::pmn_f_stat_sW<TFloat>(mat, n_dims,
                                          permutted_groupings, max_p-tp,
                                          inv_group_sizes,
                                          permutted_sWs+tp);
#endif
      }
  }

#if defined(UNIFRAC_ENABLE_ACC_NV)
  if (skbio_use_acc==ACC_NV) {
    su_acc_nv::acc_destroy_buf(inv_group_sizes,0,n_groups);
    // must remove const, as it will indeed destroy the copy in GPU memory
    su_acc_nv::acc_destroy_buf(const_cast<TFloat*>(mat),0,mat_size);
  }
#endif
#if defined(UNIFRAC_ENABLE_ACC_AMD)
  if (skbio_use_acc==ACC_AMD) {
    su_acc_amd::acc_destroy_buf(inv_group_sizes,0,n_groups);
    // must remove const, as it will indeed destroy the copy in GPU memory
    su_acc_amd::acc_destroy_buf(const_cast<TFloat*>(mat),0,mat_size);
  }
#endif

  delete[] inv_group_sizes;
  delete[] permutted_groupings;
}

// Compute the square sum of the upper triangle
template<class TFloat>
inline TFloat sum_upper_square(const TFloat * mat, const uint32_t n_dims) {
  // Use full precision for intermediate compute, to minimize accumulation errors
  double sum = 0.0;
  // not optimal parallelism, but this is cheap anyway
#pragma omp parallel for shared(mat) reduction(+:sum)
  for (uint32_t row=0; row < (n_dims-1); row++) {   // no columns in last row
    const TFloat * mat_row = mat + uint64_t(row)*uint64_t(n_dims);
    for (uint32_t col=row+1; col < n_dims; col++) { // diagonal is always zero
          TFloat val = mat_row[col]; // mat[row,col];
          sum+=val*val;
    } // for tcol
  } // for trow
  return sum;
}

// Compute the permanova values for all of the permutations
// mat is symmetric matrix of size n_dims x n_dims
// grouping is an array of size n_dims
//
// Results in permutted_fstats, and array of size (n_perm+1)
template<class TFloat>
inline void permanova_all_T(const TFloat * mat, const uint32_t n_dims,
                            const uint32_t *grouping, 
                            const uint32_t n_perm,
                            TFloat *permutted_fstats) {
  // first count the elements in the grouping
  uint32_t n_groups = (*std::max_element(grouping,grouping+n_dims)) + 1;
  uint32_t *group_sizes = new uint32_t[n_groups];
  for (uint32_t i=0; i<n_groups; i++) group_sizes[i] = 0;
  for (uint32_t i=0; i<n_dims; i++) group_sizes[grouping[i]]++;

  // compute the pseudo-F partial statistics
  {
    // Use the same buffer as the output
    TFloat *permutted_sWs = permutted_fstats;
    permanova_perm_fp_sW_T<TFloat>(mat,n_dims,grouping,
                                   group_sizes,n_groups,
                                   n_perm,
                                   permutted_sWs);
  }

  // get the normalization factor
  TFloat s_T = sum_upper_square<TFloat>(mat,n_dims)/n_dims;
 
  TFloat inv_ngroups_1 = TFloat(1.0)/ (n_groups - 1);
  TFloat inv_dg =   TFloat(1.0)/  (n_dims - n_groups);
  for (uint32_t i=0; i<(n_perm+1); i++) {
     // reminder permutted_sWs == permutted_fstats
     TFloat s_W = permutted_fstats[i];
     TFloat s_A = s_T - s_W;
     permutted_fstats[i] = (s_A * inv_ngroups_1) / (s_W * inv_dg);
  }

  delete[] group_sizes;
}

// Compute the permanova values for all of the permutations
// mat is symmetric matrix of size n_dims x n_dims
// grouping is an array of size n_dims
//
// Results in permutted_fstats, and array of size (n_perm+1)
template<class TFloat>
inline void permanova_T(const TFloat * mat, const uint32_t n_dims,
                        const uint32_t *grouping, 
                        const uint32_t n_perm,
                        TFloat &fstat, TFloat &pvalue) {
  // First compute all the permutations
  TFloat *permutted_fstats = new TFloat[n_perm+1];
  permanova_all_T<TFloat>(mat,n_dims,grouping,n_perm,permutted_fstats);

  // keep the first one and compute p_value, too
  TFloat myfstat = permutted_fstats[0];
  fstat = myfstat;
  if (n_perm>0) {
    uint32_t count_larger = 0;
    for (uint32_t i=0; i<n_perm; i++) {
      if (permutted_fstats[i+1] >= myfstat) count_larger++;
    }
    pvalue = (TFloat(1.0)*(count_larger+1))/(n_perm+1);
  } else {
    pvalue = 0.0; // just to have a deterministic value
  }

  delete[] permutted_fstats;
}

void su::permanova(const double * mat, unsigned int n_dims,
                   const uint32_t *grouping,
                   unsigned int n_perm,
                   double &fstat_out, double &pvalue_out) {
  permanova_T<double>(mat, n_dims, grouping, n_perm,
                      fstat_out, pvalue_out);
}

void su::permanova(const float * mat, unsigned int n_dims,
                   const uint32_t *grouping,
                   unsigned int n_perm,
                   float &fstat_out, float &pvalue_out) {
  permanova_T<float>(mat, n_dims, grouping, n_perm,
                     fstat_out, pvalue_out);
}

