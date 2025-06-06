/*
 * Classes, methods and unction that provide skbio-like unctionality
 */

#include "skbio_alt.hpp"
#include <stdlib.h> 
#include <omp.h>

#include <random>
#include <algorithm>

#include <cblas.h>
#include <lapacke.h>

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

static std::mt19937 myRandomGenerator;


void su::set_random_seed(uint32_t new_seed) {
  myRandomGenerator.seed(new_seed);
}

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
// ======================= PCoA ========================
//

// Compute the E_matrix with means
// centered must be pre-allocated and same size as mat (n_samples*n_samples)...will work even if centered==mat
// row_means must be pre-allocated and n_samples in size
template<class TRealIn, class TReal>
inline void E_matrix_means(const TRealIn * mat, const uint32_t n_samples,               // IN
                           TReal * centered, TReal * row_means, TReal &global_mean) {   // OUT
  /*
    Compute E matrix from a distance matrix and store in temp centered matrix.

    Squares and divides by -2 the input elementwise. Eq. 9.20 in
    Legendre & Legendre 1998.

    Compute sum of the rows at the same time.
  */

  TReal global_sum = 0.0;

#pragma omp parallel for shared(mat,centered,row_means) reduction(+: global_sum)
  for (uint32_t row=0; row<n_samples; row++) {
    const TRealIn * mat_row = mat + uint64_t(n_samples)*row;
    TReal         * centered_row = centered + uint64_t(n_samples)*row;

    TReal row_sum = 0.0;

    const TReal mhalf = -0.5;
    uint32_t col=0;

#ifdef __AVX2__
    // group many together when HW supports vecotrization
    for (; (col+7)<n_samples; col+=8) {
       TReal el0 = mat_row[col  ];
       TReal el1 = mat_row[col+1];
       TReal el2 = mat_row[col+2];
       TReal el3 = mat_row[col+3];
       TReal el4 = mat_row[col+4];
       TReal el5 = mat_row[col+5];
       TReal el6 = mat_row[col+6];
       TReal el7 = mat_row[col+7];
       el0 =  mhalf*el0*el0;
       el1 =  mhalf*el1*el1;
       el2 =  mhalf*el2*el2;
       el3 =  mhalf*el3*el3;
       el4 =  mhalf*el4*el4;
       el5 =  mhalf*el5*el5;
       el6 =  mhalf*el6*el6;
       el7 =  mhalf*el7*el7;
       centered_row[col  ] = el0;
       centered_row[col+1] = el1;
       centered_row[col+2] = el2;
       centered_row[col+3] = el3;
       centered_row[col+4] = el4;
       centered_row[col+5] = el5;
       centered_row[col+6] = el6;
       centered_row[col+7] = el7;
      
       row_sum += el0 + el1 + el2 + el3 + el4 + el5 + el6 + el7; 
    }

// else if not __AVX2__
#else

#ifdef __AVX__
    for (; (col+3)<n_samples; col+=4) {
       TReal el0 = mat_row[col  ];
       TReal el1 = mat_row[col+1];
       TReal el2 = mat_row[col+2];
       TReal el3 = mat_row[col+3];
       el0 =  mhalf*el0*el0;
       el1 =  mhalf*el1*el1;
       el2 =  mhalf*el2*el2;
       el3 =  mhalf*el3*el3;
       centered_row[col  ] = el0;
       centered_row[col+1] = el1;
       centered_row[col+2] = el2;
       centered_row[col+3] = el3;
      
       row_sum += el0 + el1 + el2 + el3; 
    }

// endif __AVX__
#endif

// endif __AVX2__
#endif

    // in case there are any leftovers
    for (; col<n_samples; col++) {
       TReal el0 = mat_row[col  ];
       el0 =  mhalf*el0*el0;
       centered_row[col  ] = el0;
       row_sum += el0;
    }

    global_sum += row_sum;
    row_means[row] = row_sum/n_samples;
  }

  global_mean = (global_sum/n_samples)/n_samples;
}

// centered must be pre-allocated and same size as mat
template<class TReal>
inline void F_matrix_inplace(const TReal * __restrict__ row_means, const TReal global_mean, TReal * __restrict__ centered, const uint32_t n_samples) {
  /*
    Compute F matrix from E matrix.

    Centring step: for each element, the mean of the corresponding
    row and column are substracted, and the mean of the whole
    matrix is added. Eq. 9.21 in Legendre & Legendre 1998.
    Pseudo-code:
    row_means = E_matrix.mean(axis=1, keepdims=True)
    col_means = Transpose(row_means)
    matrix_mean = E_matrix.mean()
    return E_matrix - row_means - col_means + matrix_mean
  */

  // use a tiled pattern to maximize locality of row_means
#pragma omp parallel for shared(centered,row_means)
  for (uint32_t trow=0; trow<n_samples; trow+=512) {
    uint32_t trow_max = std::min(trow+512, n_samples);

    for (uint32_t tcol=0; tcol<n_samples; tcol+=512) {
      uint32_t tcol_max = std::min(tcol+512, n_samples);

      for (uint32_t row=trow; row<trow_max; row++) {
        TReal *  __restrict__ centered_row = centered + uint64_t(n_samples)*row;
        const TReal gr_mean = global_mean - row_means[row];

        for (uint32_t col=tcol; col<tcol_max; col++) {
          centered_row[col] += gr_mean - row_means[col];
        }
      }
    }

  }
}

// Center the matrix
// mat and center must be nxn and symmetric
// centered must be pre-allocated and same size as mat...will work even if centered==mat
template<class TRealIn, class TReal>
inline void mat_to_centered_T(const TRealIn * mat, const uint32_t n_samples, TReal * centered) {

   TReal global_mean;
   TReal *row_means = (TReal *) malloc(uint64_t(n_samples)*sizeof(TReal));
   E_matrix_means(mat, n_samples, centered, row_means, global_mean);
   F_matrix_inplace(row_means, global_mean, centered, n_samples);
   free(row_means);
}

void su::mat_to_centered(const double * mat, const uint32_t n_samples, double * centered) {
  mat_to_centered_T(mat,n_samples,centered);
}

void su::mat_to_centered(const float  * mat, const uint32_t n_samples, float  * centered) {
  mat_to_centered_T(mat,n_samples,centered);
}

void su::mat_to_centered(const double * mat, const uint32_t n_samples, float  * centered) {
  mat_to_centered_T(mat,n_samples,centered);
}

// Matrix dot multiplication
// Expects FORTRAN-style ColOrder
// mat must be   cols x rows
// other must be cols x rows (ColOrder... rows elements together)
template<class TReal>
inline void mat_dot_T(const TReal *mat, const TReal *other, const uint32_t rows, const uint32_t cols, TReal *out);

template<>
inline void mat_dot_T<double>(const double *mat, const double *other, const uint32_t rows, const uint32_t cols, double *out)
{
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans, rows , cols, rows, 1.0, mat, rows, other, rows, 0.0, out, rows);
}

template<>
inline void mat_dot_T<float>(const float *mat, const float *other, const uint32_t rows, const uint32_t cols, float *out)
{
  cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans, rows , cols, rows, 1.0, mat, rows, other, rows, 0.0, out, rows);
}

// Expects FORTRAN-style ColOrder
// Based on N. Halko, P.G. Martinsson, Y. Shkolnisky, and M. Tygert.
//     Original Paper: https://arxiv.org/abs/1007.5510
// Step 1
// centered == n x n
// randomized = k*2 x n (ColOrder... n elements together)
template<class TReal>
inline void centered_randomize_T(const TReal * centered, const uint32_t n_samples, const uint32_t k, TReal * randomized) {
  uint64_t matrix_els = uint64_t(n_samples)*uint64_t(k);
  TReal * tmp = (TReal *) malloc(matrix_els*sizeof(TReal));

  // Form a real n x k matrix whose entries are independent, identically
  // distributed Gaussian random variables of zero mean and unit variance
  TReal *G = tmp;
  {
    std::normal_distribution<TReal> distribution;
    for (uint64_t i=0; i<matrix_els; i++) G[i] = distribution(myRandomGenerator);
  }

  // Note: Using the transposed version for efficiency (COL_ORDER)
  // Since centered is symmetric, it works just fine

  //First compute the top part of H
  mat_dot_T<TReal>(centered,G,n_samples,k,randomized);

  // power method... single iteration.. store in 2nd part of output
  // Reusing tmp buffer for intermediate storage
  mat_dot_T<TReal>(centered,randomized,n_samples,k,tmp);
  mat_dot_T<TReal>(centered,tmp,n_samples,k,randomized+matrix_els);

  free(tmp);
}

// templated LAPACKE wrapper

// Compute QR
// H is in,overwritten by Q on out
// H is (r x c), Q is (r x qc), with rc<=c
template<class TReal>
inline int qr_i_T(const uint32_t rows, const uint32_t cols, TReal *H, uint32_t &qcols);

template<>
inline int qr_i_T<double>(const uint32_t rows, const uint32_t cols, double *H, uint32_t &qcols) {
  qcols= std::min(rows,cols);
  double *tau= new double[qcols];
  int rc = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, rows, cols, H, rows, tau);
  if (rc==0) {
    qcols= std::min(rows,cols);
    rc = LAPACKE_dorgqr(LAPACK_COL_MAJOR, rows, qcols, qcols, H, rows, tau);
  }
  delete[] tau;
  return rc;
}

template<>
inline int qr_i_T<float>(const uint32_t rows, const uint32_t cols, float *H, uint32_t &qcols) {
  qcols= std::min(rows,cols);
  float *tau= new float[qcols];
  int rc = LAPACKE_sgeqrf(LAPACK_COL_MAJOR, rows, cols, H, rows, tau);
  if (rc==0) {
    qcols= std::min(rows,cols);
    rc = LAPACKE_sorgqr(LAPACK_COL_MAJOR, rows, qcols, qcols, H, rows, tau);
  }
  delete[] tau;
  return rc;
}

namespace su {

// helper class, since QR ops are multi function
template<class TReal>
class QR {
  public:
    uint32_t rows;
    uint32_t cols;

    TReal *Q;

    // will take ownership of _H
    QR(const uint32_t _rows, const uint32_t _cols, TReal *_H) 
    : rows(_rows)
    , Q(_H)
    {
      int rc = qr_i_T<TReal>(_rows, _cols, Q, cols);
      if (rc!=0) {
        fprintf(stderr, "qr_i_T(_rows,_cols, H, cols) failed with %i\n", rc);
        exit(1); // should never fail
      }
    }

    ~QR() {
      free(Q);
    }

    // res = mat * Q
    // mat must be  rows x rows
    // res will be rows * cols
    void qdot_r_sq(const TReal *mat, TReal *res);

    // res = Q * mat
    // mat must be cols * cols
    // res will be rows * cols
    void qdot_l_sq(const TReal *mat, TReal *res);

};

}

template<>
inline void su::QR<double>::qdot_r_sq(const double *mat, double *res) {
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans, rows , cols, rows, 1.0, mat, rows, Q, rows, 0.0, res, rows);
}

template<>
inline void su::QR<float>::qdot_r_sq(const float *mat, float *res) {
  cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans, rows , cols, rows, 1.0, mat, rows, Q, rows, 0.0, res, rows);
}

template<>
inline void su::QR<double>::qdot_l_sq(const double *mat, double *res) {
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans, rows , cols, cols, 1.0, Q, rows, mat, cols, 0.0, res, rows);
}

template<>
inline void su::QR<float>::qdot_l_sq(const float *mat, float *res) {
  cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans, rows , cols, cols, 1.0, Q, rows, mat, cols, 0.0, res, rows);
}

// compute svd, and return S and V
// T = input
// S output
// T is Vt on output
template<class TReal>
inline int svd_it_T(const uint32_t rows, const uint32_t cols, TReal *T, TReal *S);

template<>
inline int svd_it_T<double>(const uint32_t rows, const uint32_t cols, double *T, double *S) {
  double *superb = (double *) malloc(sizeof(double)*rows);
  int res =LAPACKE_dgesvd(LAPACK_COL_MAJOR, 'N', 'O', rows, cols, T, rows, S, NULL, rows, NULL, cols, superb);
  free(superb);

  return res;
}

template<>
inline int svd_it_T<float>(const uint32_t rows, const uint32_t cols, float *T, float *S) {
  float *superb = (float *) malloc(sizeof(float)*rows);
  int res =LAPACKE_sgesvd(LAPACK_COL_MAJOR, 'N', 'O', rows, cols, T, rows, S, NULL, rows, NULL, cols, superb);
  free(superb);

  return res;
}

// square matrix transpose, with org not alingned
template<class TReal>
inline void transpose_sq_st_T(const uint64_t n, const uint64_t stride, const TReal *in, TReal *out) {
  // n expected to be small, so simple single-thread perfect
  // org_n>=n guaranteed
  for (uint64_t i=0; i<n; i++)
    for (uint64_t j=0; j<n; j++)
       out[i*n+j] = in[i + j*stride];
}

// arbitrary matrix transpose, with copy
// in  is cols x rows
// out is rows x cols
template<class TReal>
inline void transpose_T(const uint64_t rows, const uint64_t cols, const TReal *in, TReal *out) {
  // To be optimizedc
  for (uint64_t i=0; i<rows; i++)
    for (uint64_t j=0; j<cols; j++)
       out[i*cols+j] = in[i + j*rows];
}


// Based on N. Halko, P.G. Martinsson, Y. Shkolnisky, and M. Tygert.
//     Original Paper: https://arxiv.org/abs/1007.5510
// centered == n x n, must be symmetric, Note: will be used in-place as temp buffer
template<class TReal>
inline void find_eigens_fast_T(const uint32_t n_samples, const uint32_t n_dims, TReal * centered, TReal * &eigenvalues, TReal * &eigenvectors) {
  const uint32_t k = n_dims+2;

  int rc;

  TReal *S = (TReal *) malloc(uint64_t(n_samples)*sizeof(TReal));  // take worst case size as a start
  TReal *Ut = NULL;

  {
    TReal *H = (TReal *) malloc(sizeof(TReal)*uint64_t(n_samples)*uint64_t(k)*2);

    // step 1
    centered_randomize_T<TReal>(centered, n_samples, k, H);

    // step 2
    // QR decomposition of H 

    su::QR<TReal> qr_obj(n_samples, k*2, H); // H is now owned by qr_obj, as Q

    // step 3
    // T = centered * Q (since centered^T == centered, due to being symmetric)
    // centered = n x n
    // T = n x ref
    
    TReal *T = (TReal *) malloc(sizeof(TReal)*uint64_t(qr_obj.rows)*uint64_t(qr_obj.cols));
    qr_obj.qdot_r_sq(centered,T);

    // step 4
    // compute svd
    // update T in-place, Wt on output (Vt according to the LAPACK nomenclature)
    rc=svd_it_T<TReal>(qr_obj.rows,qr_obj.cols, T, S);
    if (rc!=0) {
      fprintf(stderr, "svd_it_T<TReal>(n_samples, T, S) failed with %i\n",rc);
      exit(1); // should never fail
    }

    // step 5
    // Compute U = Q*Wt^t
    {
      // transpose Wt -> W, Wt uses n_samples strides
      TReal * W = (TReal *) malloc(sizeof(TReal)*uint64_t(qr_obj.cols)*uint64_t(qr_obj.cols));
      transpose_sq_st_T<TReal>(qr_obj.cols, qr_obj.rows, T, W);  // Wt == T on input

      Ut = T; // Ut takes ownership of the T buffer
      qr_obj.qdot_l_sq(W, Ut);

      free(W);
    }

  } // we don't need qr_obj anymore, release memory

  // step 6
  // get the interesting subset, and return
  
  // simply truncate the values, since it is a vector
  eigenvalues  = (TReal *) realloc(S, sizeof(TReal)*n_dims);

  // *eigenvectors = U = Vt
  // use only the truncated part of W, then transpose
  TReal *U = (TReal *) malloc(uint64_t(n_samples)*uint64_t(n_dims)*sizeof(TReal));

  transpose_T<TReal>(n_samples, n_dims, Ut, U);
  eigenvectors = U;

  free(Ut);
}

void su::find_eigens_fast(const uint32_t n_samples, const uint32_t n_dims, double * centered, double * &eigenvalues, double * &eigenvectors) {
  find_eigens_fast_T<double>(n_samples, n_dims, centered, eigenvalues, eigenvectors);
}

void su::find_eigens_fast(const uint32_t n_samples, const uint32_t n_dims, float * centered, float * &eigenvalues, float * &eigenvectors) {
  find_eigens_fast_T<float>(n_samples, n_dims, centered, eigenvalues, eigenvectors);
}

// helper class

namespace su {

template<class TReal>
class NewCentered {
private:
   const uint32_t n_samples;
   const uint32_t n_dims;
   TReal * centered_buf;
public:
   NewCentered(const uint32_t _n_samples, const uint32_t _n_dims) 
   : n_samples(_n_samples)
   , n_dims(_n_dims)
   , centered_buf(NULL)
   {}

   TReal * get_buf() {
     if (centered_buf==NULL) centered_buf = (TReal *) malloc(sizeof(TReal)*uint64_t(n_samples)*uint64_t(n_samples));
     return centered_buf;
   }

   void release_buf() {
     if (centered_buf!=NULL) free(centered_buf);
     centered_buf=NULL;
   }

   ~NewCentered() {
     if (centered_buf!=NULL) release_buf();
   }

private:
   NewCentered(const NewCentered<TReal> &other) = delete;
   NewCentered<TReal>& operator=(const NewCentered<TReal> &other) = delete;
}; 

template<class TReal>
class InPlaceCentered {
private:
   TReal * mat;
public:
   InPlaceCentered(TReal * _mat)
   : mat(_mat)
   {}

   TReal * get_buf() { return mat; }

   void release_buf() {}

   ~InPlaceCentered() {}
};

}

/*
    Perform Principal Coordinate Analysis.

    Principal Coordinate Analysis (PCoA) is a method similar
    to Principal Components Analysis (PCA) with the difference that PCoA
    operates on distance matrices, typically with non-euclidian and thus
    ecologically meaningful distances like UniFrac in microbiome research.

    In ecology, the euclidean distance preserved by Principal
    Component Analysis (PCA) is often not a good choice because it
    deals poorly with double zeros (Species have unimodal
    distributions along environmental gradients, so if a species is
    absent from two sites at the same site, it can't be known if an
    environmental variable is too high in one of them and too low in
    the other, or too low in both, etc. On the other hand, if an
    species is present in two sites, that means that the sites are
    similar.).

    Note that the returned eigenvectors are not normalized to unit length.
*/

// mat       - in, result of unifrac compute
// inplace   - in, if true, use mat as a work buffer 
// n_samples - in, size of the matrix (n x n)
// n_dims    - in, Dimensions to reduce the distance matrix to. This number determines how many eigenvectors and eigenvalues will be returned.
// eigenvalues - out, alocated buffer of size n_dims
// samples     - out, alocated buffer of size n_dims x n_samples
// proportion_explained - out, allocated buffer of size n_dims

template<class TRealIn, class TReal, class TCenter>
inline void pcoa_T(TRealIn * mat, TCenter &center_obj, const uint32_t n_samples, const uint32_t n_dims, TReal * &eigenvalues, TReal * &samples,TReal * &proportion_explained) {
  proportion_explained = (TReal *) malloc(sizeof(TReal)*n_dims);

  TReal diag_sum = 0.0;
  TReal *eigenvectors = NULL;

  {
    TReal *centered = center_obj.get_buf();

    // First must center the matrix
    mat_to_centered_T<TRealIn,TReal>(mat,n_samples,centered);

    // get the sum of the diagonal, needed later
    // and centered will be updated in-place in find_eigen
    for (uint32_t i=0; i<n_samples; i++) diag_sum += centered[i*uint64_t(n_samples)+i];

    // Find eigenvalues and eigenvectors
    // Use the Fast method... will return the allocated buffers
    eigenvalues = NULL;
    eigenvectors = NULL;
    find_eigens_fast_T<TReal>(n_samples,n_dims,centered,eigenvalues,eigenvectors);

    center_obj.release_buf();
  }

  // expects eigenvalues to be ordered and non-negative
  // The above unction guarantees that


  // Scale eigenvalues to have length = sqrt(eigenvalue). This
  // works because np.linalg.eigh returns normalized
  // eigenvectors. Each row contains the coordinates of the
  // objects in the space of principal coordinates. Note that at
  // least one eigenvalue is zero because only n-1 axes are
  // needed to represent n points in a euclidean space.
  // samples = eigvecs * np.sqrt(eigvals) 
  // we will  just update in place and pass out
  samples = eigenvectors;

  // use proportion_explained as tmp buffer here
  {
    TReal *sqvals = proportion_explained;
    for (uint32_t i=0; i<n_dims; i++) sqvals[i]= sqrt(eigenvalues[i]);

    // we will  just update in place and pass out
    samples = eigenvectors;

#pragma omp parallel for default(shared)
    for (uint32_t row=0; row<n_samples; row++) {
      TReal *prow = samples+(row*uint64_t(n_dims));
      for (uint32_t i=0; i<n_dims; i++) prow[i] *= sqvals[i];
    }
  }

  // now compute the real proportion_explained
  for (uint32_t i=0; i<n_dims; i++) proportion_explained[i] = eigenvalues[i]/diag_sum;

}

void su::pcoa(const double * mat, const uint32_t n_samples, const uint32_t n_dims, double * &eigenvalues, double * &samples, double * &proportion_explained) {
  su::NewCentered<double> cobj(n_samples, n_dims);
  pcoa_T(mat, cobj , n_samples, n_dims, eigenvalues, samples, proportion_explained);
}

void su::pcoa(const float  * mat, const uint32_t n_samples, const uint32_t n_dims, float  * &eigenvalues, float  * &samples, float  * &proportion_explained) {
  su::NewCentered<float> cobj(n_samples, n_dims);
  pcoa_T(mat, cobj, n_samples, n_dims, eigenvalues, samples, proportion_explained);
}

void su::pcoa(const double * mat, const uint32_t n_samples, const uint32_t n_dims, float  * &eigenvalues, float  * &samples, float  * &proportion_explained) {
  su::NewCentered<float> cobj(n_samples, n_dims);
  pcoa_T(mat, cobj, n_samples, n_dims, eigenvalues, samples, proportion_explained);
}

void su::pcoa_inplace(double * mat, const uint32_t n_samples, const uint32_t n_dims, double * &eigenvalues, double * &samples, double * &proportion_explained) {
  su::InPlaceCentered<double> cobj(mat);
  pcoa_T(mat, cobj, n_samples, n_dims, eigenvalues, samples, proportion_explained);
}

void su::pcoa_inplace(float  * mat, const uint32_t n_samples, const uint32_t n_dims, float  * &eigenvalues, float  * &samples, float  * &proportion_explained) {
  su::InPlaceCentered<float> cobj(mat);
  pcoa_T(mat, cobj, n_samples, n_dims, eigenvalues, samples, proportion_explained);
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
  std::mt19937 *randomGenerators = new std::mt19937[step_perms];

  for (uint32_t grouping_el=0; grouping_el < step_perms; grouping_el++) {
    auto new_seed = myRandomGenerator();
    randomGenerators[grouping_el].seed(new_seed);
  }

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
           std::shuffle(my_grouping, my_grouping+n_dims, randomGenerators[grouping_el]);
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
  delete[] randomGenerators;
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

// ======================= skbio_biom_subsampled  ================================


su::skbio_biom_subsampled::skbio_biom_subsampled(const biom_inmem &parent, const bool w_replacement, const uint32_t n)
 : su::biom_subsampled(parent, w_replacement, n, uint32_t(myRandomGenerator()))
{}

