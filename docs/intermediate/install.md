# Software and package installation on Bianca


!!! info "Objectives" 

    - We'll go through the methods to install packages and tools
    - We'll briefly get an overview of the hardware on Bianca

## The module system

???- tip "Forgot how to use the module system?"

    See [the basic Bianca course page 'Using the module system'](../modules1.md).

- There is a **lot of programs and tools installed as modules** on Bianca.
- These have typically been **installed on Rackham** and is **synced over to Bianca a couple of times per day**.
- **You can request installations** but that may take **several days or weeks** to be handled by the application experts at UPPMAX.
- But you may be able to do **installations yourself**. Here the use of Rackham comes handy because of the:
    - internet connection
    - the computer architecture is somewhat similar such that precompiled binaries or compiled programs (x86_64) on Rackham will most often work also on Bianca.
    - you can use the **``wharf`` to transfer source files and binaries to Bianca from Rackham**

## Install software yourself

- If not available on Bianca already (like Conda repositories) you may have to use the ``wharf`` to install your tools
    - Alternatively let an Application Expert install the tool as a module.

!!! note "Typical workflow for installation"

    - **Download the source code or binary** (Linux on x86 and 64-bit) to Rackham first
    - **Transfer** to the ``wharf``
    - Then, either 
        - You can install in your home directory.
            - This is handy for personal needs and low numbers of files — i.e. not Conda.
        - Usually better to install in project directory.
            - This way the project contains both data and software — good for reproducibility, collaboration, and everyone's general sanity.
    - Binaries for Linux on x86 and 64-bit should be able to be run directly as it is, see the software specific installation documentation.
    - or build from source, see next session.
     

### Build from source
- To build from source use a **compiler module**
- We have several compiler versions from GNU and INTEL
     - Check with: ``$ ml avail gcc`` and ``$ ml avail intel``
- ``make`` is installed on the system
    - :warning: It could happen that the "Makefile" contains web fetching, which will not work from Bianca.
    - Usually it is not a problem to build on Rackham and move to Bianca.
- ``cmake`` is available as module
     - Check with: ``$ ml avail cmake``
- [Guide for compiling **serial** programs](https://www.uppmax.uu.se/support/user-guides/compiling-source-code/){:target="_blank"}
- [Guide for compiling **parallel** programs](https://www.uppmax.uu.se/support/user-guides/mpi-and-openmp-user-guide/){:target="_blank"}
    - [Available **combinations** of compilers and parallel libraries](https://www.uppmax.uu.se/support/user-guides/mpi-and-openmp-user-guide/#tocjump_48302061903476823_2){:target="_blank"}


???- info "About CPU hardware on Bianca"

    - Architecture:          **x86_64**
        - Intel Xeon E5-2630 v3 Huawei XH620 V3 nodes
        - Advanced Vector Extensions 2 (**AVX2**)
    - CPU op-mode(s):        32-bit, 64-bit
    - Byte Order:            Little Endian
    - CPU(s):                16
    - Thread(s) per core:    1
    - Core(s) per socket:    8
    - Socket(s):             2
    - NUMA node(s):          2
    - Model name:            Intel Core Processor (Haswell, no TSX, IBRS)
    - CPU MHz:               2394.446
    - For more info, type: ``lscpu`` in the terminal 


!!! info "Own development and Git"

    - [Own development and git](https://uppmax.github.io/bianca_workshop/extra/devel/)



## Packages and libraries to scripting programs

- Python, R and Julia all have some **centrally installed packages** that are available from the modules. 
- R has a special module called ``R_packages``, and some Machine Learning python packages are included in the ``python_ml_packages`` module.
- If not found there you can try to install those by yourself.

## Install packages or not? Check it!

### Python 

- Check python versions: ``ml avail python``
- load a python version: ``ml python/3.10.8``
- from the Python shell with the ``import`` command
- from BASH shell with the 
	
- ``pip list`` command 
- ``ml help python/3.9.5`` at UPPMAX

**Is it not there? Then proceed!**

!!! info "Tip Python packages"

    - Try Conda first directly on Bianca. We have mirrored all _major_ Conda repositories directly on UPPMAX, on both Rackham and Bianca. These are updated every third day.
    - If you want to keep number of files down, use PyPI (pip), but then you need to use Rackham and the ``wharf``.

### R

- On UPPMAX the module ``R_packages`` is an omnibus package library containing almost all packages in the CRAN and BioConductor repositories. 
- As of 2023-05-31, there were a total of 23100 R packages installed in ``R_packages/4.2.1``.
    -  A total of 23109 packages were available in CRAN and BioConductor, and 23000 of these were installed in ``R_packages/4.2.1``
    -  The additional 100 R packages available in this module were installed from the CRAN/BioConductor archives, or were hosted on github, gitlab or elsewhere.

Chances are good the R packages you need are already available once you load this module.  You can quickly check by loading it:

``$ ml R_packages/4.2.1``

Then within R, try loading the package you want, like ``glmnet``:

``library(glmnet)``

**Is it not there? Then proceed!**

## Install packages! 

### Conda

- We have mirrored all major Conda repositories directly on UPPMAX, on both Rackham and Bianca. These are updated every third day.

!!! info "Available Conda channels"
      
    - bioconda
    - biocore
    - conda-forge
    - dranew
    - free
    - main
    - pro
    - qiime2
    - r
    - r2018.11
    - scilifelab-lts

!!! info "More info"

    - [Extra material: Installing Conda packages](https://uppmax.github.io/bianca_workshop/extra/conda/)
    - [Conda user guide](https://www.uppmax.uu.se/support/user-guides/conda-user-guide/)
    - [UPPMAX Python user guide: Conda](https://www.uppmax.uu.se/support/user-guides/python-user-guide/#tocjump_9332829429720808_6)



### Python packages with pip

!!! info "Installation principle"

    - install on Rackham
        - ``pip install --user <package>``
        - ``python setup.py install --user or --prefix=<path>``
    - sync to ``wharf``
    - move the files on Bianca to correct place
    - you may have to update ``$PYTHONPATH``

!!! info "More info"

    - [Extra material: Installing pip packages](https://uppmax.github.io/bianca_workshop/extra/pip/){:target="_blank"}
    - [UPPMAX Python user guide: Pip](https://www.uppmax.uu.se/support/user-guides/python-user-guide/#tocjump_9332829429720808_5){:target="_blank"}
    - [From Python course: packages](https://uppmax.github.io/R-python-julia-HPC/python/packages.html){:target="_blank"}
    - [From Python course: isolated environments](https://uppmax.github.io/R-python-julia-HPC/python/isolated.html){:target="_blank"}


### R packages

!!! info "Installation principle"

    - install on Rackham
    - sync to ``wharf``
    - move the files on Bianca

!!! info "More info"

    - [Extra material: Installing R packages](https://uppmax.github.io/bianca_workshop/extra/rpackages/)
    - [From R course: packages](https://uppmax.github.io/R-python-julia-HPC/R/packagesR.html){:target="_blank"}
    - [From R course: isolated environments](https://uppmax.github.io/R-python-julia-HPC/R/isolatedR.html){:target="_blank"}

### Julia packages

- At UPPMAX there is a central library with installed packages.
- This is good, especially when working on Bianca, since you then do not need to install via the ``wharf``.
- A selection of the Julia packages and libraries installed on UPPMAX are:

        CSV
        CUDA
        MPI
        Distributed
        IJulia
        Plots
        PyPlot
        DataFrames

!!! info "Installation principle"

    - install on Rackham
    - sync to ``wharf``
    - move the files on Bianca

!!! info "More info"

    - [Extra material: Installing Julia packages](https://uppmax.github.io/bianca_workshop/extra/julia/){:target="_blank"}
    - [Julia course: isolated environments](https://uppmax.github.io/R-python-julia-HPC/julia/isolatedJulia.html){:target="_blank"}

## "Containers"

!!! info
   
    - Containers let you install programs without needing to think about the computer environment, like    
        - operative system
        - dependencies (libraries and other programs) with correct versions
    - 2(3) types
        - Singularity/Apptainer
        - Docker that does not work on HPC-systems
            - But docker images can be used byt Singularity and Apptainer
    - Everything is included
    - Draw-backs
        - you install also things that may be already installed
        - therefore, probably more disk space is needed

!!! info "More info"

    - [Extra material: Containers](https://uppmax.github.io/bianca_workshop/extra/containers/)

## Demo session

!!! example "What do you want to type-along"

    - [Extra material: Installing Conda packages](https://uppmax.github.io/bianca_workshop/extra/conda/)
    - [Extra material: Installing pip packages](https://uppmax.github.io/bianca_workshop/extra/pip/){:target="_blank"}
    - [Extra material: Installing R packages](https://uppmax.github.io/bianca_workshop/extra/rpackages/)
    - [Extra material: Installing Julia packages (a little imature)](https://uppmax.github.io/bianca_workshop/extra/julia/){:target="_blank"}
    - [Extra material: Containers](https://uppmax.github.io/bianca_workshop/extra/containers/)

!!! abstract "Keypoints"
    - You have got an overview of the procedures to install packages/libraries and tools on Bianca through the ``wharf``
    - If you feel uncomfortable or think that many users would benefit from the software, ask the support to install it.
