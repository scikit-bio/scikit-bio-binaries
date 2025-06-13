/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2025, UniFrac development team.
 * Copyright (c) 2025--, scikit-bio development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#include "ordination/principal_coordinate_analysis.hpp"
#include <unistd.h>
#include "tests/test_helper.hpp"

void test_center_mat() {
    SUITE_START("test center mat");

    // unweighted unifrac of test.biom
    double matrix[] = { 
      0.0000000000, 0.2000000000, 0.5714285714, 0.6000000000, 0.5000000000, 0.2000000000,
      0.2000000000, 0.0000000000, 0.4285714286 ,0.6666666667, 0.6000000000, 0.3333333333,
      0.5714285714, 0.4285714286, 0.0000000000, 0.7142857143, 0.8571428571, 0.4285714286,
      0.6000000000, 0.6666666667, 0.7142857143, 0.0000000000, 0.3333333333, 0.4000000000,
      0.5000000000, 0.6000000000, 0.8571428571, 0.3333333333, 0.0000000000, 0.6000000000,
      0.2000000000, 0.3333333333, 0.4285714286, 0.4000000000, 0.6000000000, 0.0000000000};

    const uint32_t n_samples = 6;

    double exp[] = { 0.05343726,  0.04366213, -0.0329743 , -0.07912698, -0.00495654, 0.01995843,
                     0.04366213,  0.073887  ,  0.04867914, -0.11112434, -0.04973167,-0.00537226,
                    -0.0329743 ,  0.04867914,  0.20714475, -0.07737528, -0.17044974, 0.02497543,
                    -0.07912698, -0.11112434, -0.07737528,  0.14830877,  0.11192366, 0.00739418,
                    -0.00495654, -0.04973167, -0.17044974,  0.11192366,  0.18664966,-0.07343537,
                     0.01995843, -0.00537226,  0.02497543,  0.00739418, -0.07343537, 0.02647959 };

    // fp64
    {
      double *centered = (double *) malloc(6*6*sizeof(double)); 

      skbb::mat_to_centered(n_samples, matrix, centered);

      for(int i = 0; i < (6*6); i++) {
        //printf("%i %f %f\n",i,float(centered[i]),float(exp[i]));
        ASSERT(fabs(centered[i] - exp[i]) < 0.000001);
      }

      free(centered);
    }

    // mixed
    {
      float *centered_fp32 = (float *) malloc(6*6*sizeof(float));

      skbb::mat_to_centered(n_samples, matrix, centered_fp32);

      for(int i = 0; i < (6*6); i++) {
        //printf("%i %f %f\n",i,float(centered_fp32[i]),float(exp[i]));
        ASSERT(fabs(centered_fp32[i] - exp[i]) < 0.000001);
      }

      free(centered_fp32);
    }

    float *matrix_fp32 = (float *) malloc(6*6*sizeof(float));
    for(int i = 0; i < (6*6); i++) matrix_fp32[i] = matrix[i];


    // fp32
    {
      float *centered_fp32 = (float *) malloc(6*6*sizeof(float));

      skbb::mat_to_centered(n_samples, matrix_fp32, centered_fp32);

      for(int i = 0; i < (6*6); i++) {
        //printf("%i %f %f\n",i,float(centered_fp32[i]),float(exp[i]));
        ASSERT(fabs(centered_fp32[i] - exp[i]) < 0.000001);
      }

      free(centered_fp32);
    }

    // fp32 in-place
    {
      skbb::mat_to_centered(n_samples, matrix_fp32, matrix_fp32);

      for(int i = 0; i < (6*6); i++) {
        //printf("%i %f %f\n",i,float(matrix_fp32[i]),float(exp[i]));
        ASSERT(fabs(matrix_fp32[i] - exp[i]) < 0.000001);
      }
    }

    free(matrix_fp32);


    // fp64 in-place
    {
      skbb::mat_to_centered(n_samples, matrix, matrix);

      for(int i = 0; i < (6*6); i++) {
        //printf("%i %f %f\n",i,float(matrix_fp32[i]),float(exp[i]));
        ASSERT(fabs(matrix[i] - exp[i]) < 0.000001);
      }
    }
    SUITE_END();
}

void test_pcoa_fsvd() {
    SUITE_START("test pcoa<fsvd>");

    // unweighted unifrac of crawford.biom
    double matrix[] = { 
      0.         , 0.71836067 , 0.71317361 , 0.69746044 , 0.62587207 , 0.72826674
    , 0.72065895 , 0.72640581 , 0.73606053, 
      0.71836067 , 0.         , 0.70302967 , 0.73407301 , 0.6548042  , 0.71547381
    , 0.78397813 , 0.72318399 , 0.76138933,
      0.71317361 , 0.70302967 , 0.         , 0.61041275 , 0.62331299 , 0.71848305
    , 0.70416337 , 0.75258475 , 0.79249029,
      0.69746044 , 0.73407301 , 0.61041275 , 0.         , 0.6439278  , 0.70052733
    , 0.69832716 , 0.77818938 , 0.72959894,
      0.62587207 , 0.6548042  , 0.62331299 , 0.6439278  , 0.         , 0.75782689
    , 0.71005144 , 0.75065046 , 0.78944369,
      0.72826674 , 0.71547381 , 0.71848305 , 0.70052733 , 0.75782689 , 0.
    , 0.63593642 , 0.71283615 , 0.58314638,
      0.72065895 , 0.78397813 , 0.70416337 , 0.69832716 , 0.71005144 , 0.63593642
    , 0.         , 0.69200762 , 0.68972056,
      0.72640581 , 0.72318399 , 0.75258475 , 0.77818938 , 0.75065046 , 0.71283615
    , 0.69200762 , 0.         , 0.71514083,
      0.73606053 , 0.76138933 , 0.79249029 , 0.72959894 , 0.78944369 , 0.58314638
    , 0.68972056 , 0.71514083 , 0. };


    const uint32_t n_samples = 9;

    // Test centering

    double exp2[] = { 
        0.22225336 , -0.025481   , -0.03491711 , -0.02614685  , 0.01899659 , -0.05103589
     , -0.03971812 , -0.02699392 , -0.03695706,
       -0.025481   ,  0.24282669 , -0.01744751 , -0.04206624  , 0.0107569  , -0.03151438
     , -0.07706765 , -0.0143721  , -0.0456347,
       -0.03491711 , -0.01744751 ,  0.21652901 ,  0.02791465  , 0.0177328  , -0.04682078
     , -0.03082866 , -0.04921529 , -0.08294711,
       -0.02614685 , -0.04206624 ,  0.02791465 ,  0.21190401  , 0.00235833 , -0.03639361
     , -0.02904855 , -0.07112525 , -0.03739649,
        0.01899659 ,  0.0107569  ,  0.0177328  ,  0.00235833  , 0.20745566 , -0.08039931
     , -0.03952883 , -0.05229812 , -0.08507403,
       -0.05103589 , -0.03151438 , -0.04682078 , -0.03639361  , -0.08039931 , 0.20604732
     ,  0.00964595 , -0.02533192 ,  0.05580262,
       -0.03971812 , -0.07706765 , -0.03082866 , -0.02904855  , -0.03952883 , 0.00964595
     ,  0.21765972 , -0.00489531 , -0.00621856,
       -0.02699392 , -0.0143721  , -0.04921529 , -0.07112525  , -0.05229812 , -0.02533192
     , -0.00489531 ,  0.2514242  , -0.00719229,
       -0.03695706 , -0.0456347  , -0.08294711 , -0.03739649  , -0.08507403 , 0.05580262
     , -0.00621856 , -0.00719229 , 0.24561762 };

    double *centered = (double *) malloc(9*9*sizeof(double));
    float *centered_fp32 = (float *) malloc(9*9*sizeof(float));

    skbb::mat_to_centered(n_samples, matrix, centered);

    for(int i = 0; i < (9*9); i++) {
      //printf("%i %f %f\n",i,float(centered[i]),float(exp[i]));
      ASSERT(fabs(centered[i] - exp2[i]) < 0.00001);
    }

    skbb::mat_to_centered(n_samples, matrix, centered_fp32);

    for(int i = 0; i < (9*9); i++) {
      //printf("%i %f %f\n",i,float(centered_fp32[i]),float(exp[i]));
      ASSERT(fabs(centered_fp32[i] - exp2[i]) < 0.00001);
    }

    // Test eigens

    double exp3a[] = {0.45752162, 0.3260088 , 0.2791141 , 0.26296948, 0.20924533};
    double exp3b[] = {
       -0.17316152,  0.17579996,  0.23301609, -0.74519625, -0.05194624,
       -0.19959264,  0.53235665, -0.53370018,  0.2173474 ,  0.26736004,
       -0.35794942, -0.27956624,  0.01114096,  0.40488848, -0.13121464,
       -0.2296467 , -0.47494333, -0.12571292,  0.02313551, -0.46916459,
       -0.44501584,  0.05597451,  0.07717711, -0.15881922,  0.24594442,
        0.40335552, -0.16290597, -0.30327343,  0.03778646,  0.21664806,
        0.23769142, -0.29034629,  0.46813757,  0.13858945,  0.58624834,
        0.2407584 ,  0.51300752,  0.48211607,  0.34422672, -0.42416046,
        0.52356078, -0.0693768 , -0.30890127, -0.26195855, -0.23971493};

    {
      double eigenvalues[5];
      double eigenvectors[5*9];
      skbb::find_eigens_fast(n_samples, centered, 5, eigenvalues, eigenvectors);

      for(int i = 0; i < 5; i++) {
        //printf("%i %f %f\n",i,float(eigenvalues[i]),float(exp3a[i]));
        ASSERT(fabs(eigenvalues[i] - exp3a[i]) < 0.00001);
      }

      // signs may flip, that's normal
      for(int i = 0; i < (5*9); i++) {
        //printf("%i %f %f %f\n",i,float(eigenvectors[i]),float(exp3b[i]),float(fabs(eigenvectors[i]) - fabs(exp3b[i])));
        ASSERT( fabs(fabs(eigenvectors[i]) - fabs(exp3b[i])) < 0.00001);
      }
    
    }
    {
      float eigenvalues[5];
      float eigenvectors[5*9];
      skbb::find_eigens_fast(n_samples, centered_fp32, 5, eigenvalues, eigenvectors);

      for(int i = 0; i < 5; i++) {
        //printf("%i %f %f\n",i,float(eigenvalues[i]),float(exp3a[i]));
        ASSERT(fabs(eigenvalues[i] - exp3a[i]) < 0.00001);
      }

      // signs may flip, that's normal
      for(int i = 0; i < (5*9); i++) {
        //printf("%i %f %f %f\n",i,float(eigenvectors[i]),float(exp3b[i]),float(fabs(eigenvectors[i]) - fabs(exp3b[i])));
        ASSERT( fabs(fabs(eigenvectors[i]) - fabs(exp3b[i])) < 0.00001);
      }
    
    }

    free(centered_fp32);
    free(centered);


    // Test PCoA (incudes the above calls

    double *exp4a = exp3a; // same eigenvals;
    double exp4b[] = {
       -0.11712705,  0.10037682, -0.12310531, -0.38214073, -0.02376195,
       -0.13500515,  0.30396064,  0.28196047,  0.11145694,  0.12229942,
       -0.24211822, -0.15962444, -0.00588591,  0.20762904, -0.06002196,
       -0.15533382, -0.27117925,  0.06641571,  0.01186402, -0.21461156,
       -0.30101024,  0.03195987, -0.04077363, -0.08144337,  0.1125032 ,
        0.27283106, -0.09301471,  0.16022314,  0.0193771 ,  0.09910206,
        0.16077529, -0.16577955, -0.24732293,  0.07106943,  0.26816958,
        0.16284981,  0.29291283, -0.25470794,  0.17652136, -0.19402517,
        0.35413832, -0.0396122 ,  0.1631964 , -0.13433378, -0.10965362};
    double exp4c[] = {0.22630343, 0.16125338, 0.13805791, 0.13007231, 0.10349879};

    {
      double *eigenvalues;
      double *samples;
      double *proportion_explained;

      skbb::pcoa_fsvd(n_samples, matrix, 5, eigenvalues, samples, proportion_explained); 

      for(int i = 0; i < 5; i++) {
        //printf("%i %f %f\n",i,float(eigenvalues[i]),float(exp4a[i]));
        ASSERT(fabs(eigenvalues[i] - exp4a[i]) < 0.00001);
      }

      // signs may flip, that's normal
      for(int i = 0; i < (5*9); i++) {
        //printf("%i %f %f %f\n",i,float(samples[i]),float(exp4b[i]),float(fabs(samples[i]) - fabs(exp4b[i])));
        ASSERT( fabs(fabs(samples[i]) - fabs(exp4b[i])) < 0.00001);
      }

      for(int i = 0; i < 5; i++) {
        //printf("%i %f %f\n",i,float(proportion_explained[i]),float(exp4c[i]));
        ASSERT(fabs(proportion_explained[i] - exp4c[i]) < 0.00001);
      } 

      free(eigenvalues);
      free(samples);
      free(proportion_explained);
    }

    // Test PCoA mixed mode
    {
      float *eigenvalues_fp32;
      float *samples_fp32;
      float *proportion_explained_fp32;

      skbb::pcoa_fsvd(n_samples, matrix, 5, eigenvalues_fp32, samples_fp32, proportion_explained_fp32);

      for(int i = 0; i < 5; i++) {
        //printf("%i %f %f\n",i,float(eigenvalues_fp32[i]),float(exp4a[i]));
       ASSERT(fabs(eigenvalues_fp32[i] - exp4a[i]) < 0.00001);
      }

      // signs may flip, that's normal
      for(int i = 0; i < (5*9); i++) {
        //printf("%i %f %f %f\n",i,float(samples_fp32[i]),float(exp4b[i]),float(fabs(samples_fp32[i]) - fabs(exp4b[i])));
        ASSERT( fabs(fabs(samples_fp32[i]) - fabs(exp4b[i])) < 0.00001);
      }

      for(int i = 0; i < 5; i++) {
        //printf("%i %f %f\n",i,float(proportion_explained_fp32[i]),float(exp4c[i]));
        ASSERT(fabs(proportion_explained_fp32[i] - exp4c[i]) < 0.00001);
      }

      free(eigenvalues_fp32);
      free(samples_fp32);
      free(proportion_explained_fp32);
    }

    // test in-place
    {
      double *eigenvalues;
      double *samples;
      double *proportion_explained;

      skbb::pcoa_fsvd_inplace(n_samples, matrix, 5, eigenvalues, samples, proportion_explained);
      // Note: matrix content has been destroyed

      for(int i = 0; i < 5; i++) {
        //printf("%i %f %f\n",i,float(eigenvalues[i]),float(exp4a[i]));
        ASSERT(fabs(eigenvalues[i] - exp4a[i]) < 0.00001);
      }

      // signs may flip, that's normal
      for(int i = 0; i < (5*9); i++) {
        //printf("%i %f %f %f\n",i,float(samples[i]),float(exp4b[i]),float(fabs(samples[i]) - fabs(exp4b[i])));
        ASSERT( fabs(fabs(samples[i]) - fabs(exp4b[i])) < 0.00001);
      }

      for(int i = 0; i < 5; i++) {
        //printf("%i %f %f\n",i,float(proportion_explained[i]),float(exp4c[i]));
        ASSERT(fabs(proportion_explained[i] - exp4c[i]) < 0.00001);
      }

      free(eigenvalues);
      free(samples);
      free(proportion_explained);
    }

    SUITE_END();
}


int main(int argc, char** argv) {
    test_center_mat();
    test_pcoa_fsvd();

    printf("\n");
    printf(" %i / %i suites failed\n", suites_failed, suites_run);
    printf(" %i / %i suites empty\n", suites_empty, suites_run);
    printf(" %i / %i tests failed\n", tests_failed, tests_run);
  
    printf("\n THE END.\n");
    
    return tests_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
