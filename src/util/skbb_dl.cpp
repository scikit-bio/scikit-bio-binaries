/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2023-2025, UniFrac development team.
 * Copyright (c) 2025--, scikit-bio development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

// This file is intended to be included in other .cpp files

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <dlfcn.h>
#include <pthread.h>

/* Handle pointing to the approriate shared library implementing the functionality
 * Initialized on first use. */
static void *dl_handle = NULL;

static void dl_load(const char *fncname,
                     void **dl_ptr) {
   char *error;

   if (dl_handle==NULL) {
       const char* lib_name = dl_get_lib_name();
       const char* env_cpu_info = getenv("SKBB_GPU_INFO");
       if ((env_cpu_info!=NULL) && (env_cpu_info[0]=='Y')) {
           printf("INFO (skbio_bins): Using shared library %s\n",lib_name);
       }
       dl_handle = dlopen(lib_name, RTLD_LAZY);
       if (!dl_handle) {
          fputs(dlerror(), stderr);
          exit(1);
       }
   }

   *dl_ptr = dlsym(dl_handle, fncname);
   if ((error = dlerror()) != NULL)  {
       fputs(error, stderr);
       exit(1);
   }
}

static pthread_mutex_t dl_mutex = PTHREAD_MUTEX_INITIALIZER;

static void cond_dl_load(const char *fncname,
                     void **dl_ptr) {

   pthread_mutex_lock(&dl_mutex);
   if ((*dl_ptr)==NULL) dl_load(fncname,dl_ptr);
   pthread_mutex_unlock(&dl_mutex);
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
/* This function is sometimes used when included, so prevent the compiler warning */
static bool dl_load_check() {

   pthread_mutex_lock(&dl_mutex);
   if (dl_handle==NULL) {
       const char* lib_name = dl_get_lib_name();
       dl_handle = dlopen(lib_name, RTLD_LAZY);
       if (!dl_handle) {
          // no such shared library
          pthread_mutex_unlock(&dl_mutex);
	  return false;
       }
       // only print out if the library exists
       const char* env_cpu_info = getenv("SKBB_GPU_INFO");
       if ((env_cpu_info!=NULL) && (env_cpu_info[0]=='Y')) {
           printf("INFO (skbio_bins): Using shared library %s\n",lib_name);
       }
   }
   pthread_mutex_unlock(&dl_mutex);
   return true;
}

#pragma GCC diagnostic pop

