|License| |Downloads| |Install|

.. image:: https://scikit.bio/_images/logo.svg
   :width: 300 px
   :target: https://scikit.bio
   :alt: scikit-bio logo

*scikit-bio-binaries is an open-source, BSD-licensed package providing optimized algorithms for bioinformatics in support of scikit-bio (`<https://scikit.bio>`_).*


Releases
--------

There are currently no releases available.

Installation
------------

In the near future you will be able to install the latest release of scikit-bio-binaries using ``conda`` (`<https://www.anaconda.com/docs/getting-started/miniconda/main>`_)::

    conda install -c conda-forge scikit-bio-binaries

Or using ``pip``::

    pip install scikit-bio-binaries

Local build instructions
------------------------

The package can be build from source on your local machine.
We recommend using the conda-provided compilers and libraries, but system-installed ones should work as well.

If you decide to create a dedicated build environmwnt in ``conda`` (`<https://www.anaconda.com/docs/getting-started/miniconda/main>`_)::

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
.. |Downloads| image:: https://anaconda.org/conda-forge/scikit-bio-binaries/badges/downloads.svg
   :target: https://anaconda.org/conda-forge/scikit-bio-binaries
.. |Install| image:: https://anaconda.org/conda-forge/scikit-bio-binaries/badges/installer/conda.svg
   :target: https://conda.anaconda.org/conda-forge-binaries
