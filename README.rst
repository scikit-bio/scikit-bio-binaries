|License| |Version| |Downloads| |Platforms|

.. image:: https://scikit.bio/_images/logo.svg
   :width: 300 px
   :target: https://scikit.bio
   :alt: scikit-bio logo

*scikit-bio-binaries is an open-source, BSD-licensed package providing optimized algorithms for bioinformatics in support of scikit-bio.*


Installation
------------

You can install the latest release of scikit-bio-binaries using ``conda`` (`<https://www.anaconda.com/docs/getting-started/miniconda/main>`_)::

    conda install -c conda-forge scikit-bio-binaries

Alternatively, in the near future you will be able to install the latest release of scikit-bio-binaries using ``pip``::

    pip install scikit-bio-binaries

Local build instructions
------------------------

The package can be build from source on your local machine.
We recommend using the conda-provided compilers and libraries, but system-installed ones should work as well.

If you decide to create a dedicated build environment in ``conda`` (`<https://www.anaconda.com/docs/getting-started/miniconda/main>`_)::

    cplatform=`conda info |awk '/platform/{print $3}'`
    if [[ "$(uname -s)" == "Linux" ]];
    then
      conda create -n skbb-build -c conda-forge gxx_${cplatform}
    else
      conda create -n skbb-build -c conda-forge clangxx_${cplatform}
    fi 
    conda activate skbb-build
    conda install -c conda-forge libcblas liblapacke blas-devel make
    make clean && make clean_install && make all

To test that the build succeeded, run::

    make test

GPU support
~~~~~~~~~~~

To build NVIDIA-GPU-enabled code, you will also need the CUDA compiler. On Linux, you can install it from conda with::

    conda install -c conda-forge cuda-compiler

The CUDA build is not enabled by default. If you are interested in building the GPU-accelerated libraries, set::

    export NV_CUDA=Y

There is also support for OpenMP Offload and OpenACC builds, by setting either ``AMD_CXX`` or ``NV_CXX``, but it is still experimental.

Provided files and their usage
------------------------------

The functions provided by scikit-bio-binaries are packaged as a shared library and will be installed in::

    $CONDA_PREFIX/lib/libskbb.so

If GPU code was build, you will also have::

    $CONDA_PREFIX/lib/libskbb_acc_nv.so

The C header files that describe how to access them are avaialable in `<src/extern/>`_, but will also be installed in::

    $CONDA_PREFIX/include/scikit-bio-binaries/

Compiled languages can use the provided C headers during compilation, and link against the provided shared library.
See `<api_tests/>`_ for an example::

    $(CC) $(CFLAGS) my_code.c $(LDFLAGS) -lskbb -o my_exe

Python users can instead use the `ctypes <https://docs.python.org/3/library/ctypes.html>`_ module
to import the shared library at runtime. No header files are needed (or provided),
but you can use the C headers to guide your implementation.
For example::

    import ctypes
    dll = ctypes.CDLL("libskbb.so")
    skbb_version = dll.skbb_get_api_version()

Environment considerations
--------------------------

Multi-core support
~~~~~~~~~~~~~~~~~~

This package uses OpenMP to make use of multiple CPU cores.
By default, scikit-bio-binaries will use all the cores that are available on the system.
To restrict the number of cores used, set::

    export OMP_NUM_THREADS=nthreads

Older CPU support
~~~~~~~~~~~~~~~~~~

On x86_64 based CPU platforms, scikit-bio-binaries will auto-detect the CPU generation,
i.e. if it supports avx2 or avx512 vector instructions.
To force the most compatible binary variant, one can set::

    export SKBB_MAX_CPU=basic

To check which binary is used (scikit-bio-binaries will print it to standard output at runtime), set::

    export SKBB_CPU_INFO=Y

GPU support
~~~~~~~~~~~

If the code has been compiled for GPUs, scikit-bio-binaries will auto-detect the presence
of either NVIDIA or AMD GPUs, and use such a GPU for the GPU-enabled algorithms.
To force CPU-only compute, one can set::

    export SKBB_USE_GPU=N

To check if a GPU is used (scikit-bio-binaries will print it to standard output at runtime), set::

    export SKBB_GPU_INFO=Y

Additional timing information
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When evaluating the performance of scikit-bio-binaries it is sometimes necessary to distinguish
the time spent interacting with the data from the compute proper.
Additional informational messages can be enabled by setting::

    export SKBB_TIMING_INFO=Y

Adoption
--------

In the near future, ``scikit-bio-binaries`` will be used by ``scikit-bio`` and ``unifrac-binaries``.

License
-------

scikit-bio-binaries is available under the new BSD license. See `LICENSE.txt <LICENSE.txt>`_ for scikit-bio's license.


Team
----

The library is currently mainatined by **Igor Sfiligoi** at the University of California San Diego (UCSD) (@sfiligoi).
Guidance and support is also provided by 
**Dr. Qiyun Zhu** at Arizona State University (ASU) (@qiyunzhu),
**Dr. Daniel McDonald** at the University of California San Diego (UCSD) (@wasade), and
**Dr. Rob Knight** at the University of California San Diego (UCSD) (@rob-knight).


Credits
-------

The algorithms in this package are based on code developped as part of the **scikit-bio** (`<https://scikit.bio>`_) package.
See the main ``scikit-bio`` page for credits about the original algorithm contributers.


Funding
-------

The development of scikit-bio is currently supported by the U.S. Department of Energy, Office of Science under award number `DE-SC0024320 <https://genomicscience.energy.gov/compbioawards2023/#Expanding>`_, awarded to Dr. Qiyun Zhu at ASU (lead PI), Dr. James Morton at Gutz Analytics, and Dr. Rob Knight at UCSD.


Citation
--------

If you use scikit-bio derived code, including scikit-bio-binaries, for any published research, please see our `Zenodo page <https://zenodo.org/record/8209901>`_ for how to cite.


Branding
--------

The logo of scikit-bio was created by `Alina Prassas <https://cargocollective.com/alinaprassas>`_. Vector and bitmap image files are available at the `logos <logos>`_ directory.


.. |License| image:: https://anaconda.org/conda-forge/scikit-bio-binaries/badges/license.svg
   :target: https://anaconda.org/conda-forge/scikit-bio-binaries
.. |Version| image:: https://anaconda.org/conda-forge/scikit-bio-binaries/badges/version.svg
   :target: https://anaconda.org/conda-forge/scikit-bio-binaries
.. |Downloads| image:: https://anaconda.org/conda-forge/scikit-bio-binaries/badges/downloads.svg
   :target: https://anaconda.org/conda-forge/scikit-bio-binaries
.. |Platforms| image:: https://anaconda.org/conda-forge/scikit-bio-binaries/badges/platforms.svg
   :target: https://anaconda.org/conda-forge/scikit-bio-binaries
