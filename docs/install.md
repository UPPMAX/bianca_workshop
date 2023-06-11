# Software and package installation on Bianca


!!! info "Objectives" 

    - We'll go through the methods to install packages and tools
    - We'll briefly get an overview of the hardware on Bianca

## The module system

- As we have seen this morning, there is a lot of programs and tools installed as modules on Bianca.
- These have typically been installed on Rackham and is synced over to Bianca a couple of times per day.
- You can request installations but that may take several days or even weeks to be handled by the application experts at UPPMAX.
- But you may be able to do installations yourself. Here the use of Rackham comes handy because of the:
    - internet connection
    - the computer architecture is somewhat similar such that precompiled binaries or compiled programs (x86_64) on Rackham will most often work also on Bianca.
    - you can use the wharf to transfer source files and binaries to Bianca from Rackham

## Install software yourself
- You can install in your home directory.
    - This is handy for personal needs, low numbers of files (i.e. not Conda).
- Usually better to install in project directory.
    - This way the project contains both data and software — good for reproducibility, collaboration, and everyone's general sanity.
- If not available on Bianca already (like Conda repositories) you may have to use the Wharf to install your tools
    - Alternatively let an Application Expert install the tool as a module.



## Packages and libraries to scripting programs

- Python, R and Julia all have some centrally installed packages that are available from the modules. 
- R has a special module called ``R_packages``, and some Machine Learning python packages are included in the ``python_ml_packages`` module.
- If not found there you can try to install those by yourself.


!!! info "Tip Python packages"

    - Try Conda first directly on Bianca. We have mirrored all major conda repositories directly on UPPMAX, on both Rackham and Bianca. These are updated every third day.
    - If you want to keep number of files down, use PyPI (pip), but then you need to use Rackham and the wharf.

### Conda

- We have mirrored all major conda repositories directly on UPPMAX, on both Rackham and Bianca. These are updated every third day.

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

    - [Extra material: Installing Conda packages](https://uppmax.github.io/bianca_workshop/conda/)
    - [Conda user guide](https://www.uppmax.uu.se/support/user-guides/conda-user-guide/)
    - [UPPMAX Python user guide: Conda](https://www.uppmax.uu.se/support/user-guides/python-user-guide/#tocjump_9332829429720808_6)



### Python packages with pip

!!! info "Principle"

    - install on Rackham
        - pip install --user <package>
        - python setup.py install --user or --prefix=<path>
    - sync to wharf
    - move the files on Bianca
    - you may have to update $PYTHONPATH

!!! info "More info"

    - [Extra material: Installing pip packages](https://uppmax.github.io/bianca_workshop/pip/)
    - [UPPMAX Python user guide: Pip](https://www.uppmax.uu.se/support/user-guides/python-user-guide/#tocjump_9332829429720808_5)
    - [From Python course: packages](https://uppmax.github.io/R-python-julia-HPC/python/packages.html)
    - [From Python course: isolated environments](https://uppmax.github.io/R-python-julia-HPC/python/isolated.html)


### R packages

- On UPPMAX the module ``R_packages`` is an omnibus package library containing almost all packages in the CRAN and BioConductor repositories. 
- As of 2023-05-31, there were a total of 23100 R packages installed in ``R_packages/4.2.1``.
    -  A total of 23109 packages were available in CRAN and BioConductor, and 23000 of these were installed in ``R_packages/4.2.1``
    -  The additional 100 R packages available in this module were installed from the CRAN/BioConductor archives, or were hosted on github, gitlab or elsewhere.

Chances are good the R packages you need are already available once you load this module.  You can quickly check by loading it:

``$ ml R_packages/4.2.1``

Then within R, try loading the package you want:

``library(glmnet)``

Or a bit longer way, you can ``grep`` for the package after this module is loaded using the environment variable ``$R_LIBS_SITE``, which contains the locations of all R packages installed within the module.

```bash
$ ls -l $R_LIBS_SITE | grep glmnet
drwxrwsr-x  9 douglas sw  4096 May 28 16:59 EBglmnet
drwxrwsr-x 11 douglas sw  4096 May 25 01:22 glmnet
drwxrwsr-x  6 douglas sw  4096 May 25 04:03 glmnetSE
drwxrwsr-x  7 douglas sw  4096 May 25 04:04 glmnetUtils
drwxrwsr-x  8 douglas sw  4096 May 25 04:04 glmnetcr
drwxrwsr-x  7 douglas sw  4096 May 25 10:46 glmnetr
```

!!! info "More info"

    - [Extra material: Installing R packages](https://uppmax.github.io/bianca_workshop/rpackages/)
    - [From R course: packages](https://uppmax.github.io/R-python-julia-HPC/R/packagesR.html)
    - [From R course: isolated environments](https://uppmax.github.io/R-python-julia-HPC/R/isolatedR.html)

### Julia packages

- At UPPMAX there is a central library with installed packages.
- This is good, especially when working on Bianca, since you then do not need to install via the Wharf.
- A selection of the Julia packages and libraries installed on UPPMAX are:

        BenchmarkTools
        CSV
        CUDA
        MPI
        Distributed
        IJulia
        Plots
        PyPlot
        Gadfly
        DataFrames
        DistributedArrays
        PlotlyJS

!!! info "More info"

    - [Extra material: Installing Julia packages](https://uppmax.github.io/bianca_workshop/julia/)
    - [Julia course: isolated environments](https://uppmax.github.io/R-python-julia-HPC/julia/isolatedJulia.html)

## "Containers"
### Singularity
- [Singularity user guide](https://www.uppmax.uu.se/support/user-guides/singularity-user-guide/)

### Docker
- Docker will unfortunately not work on the clusters, since it requires root permission.
- However, Singularity may use Docker images.

## Build from source
- We have several compiler versions from GNU and INTEL
- check with: ``$ ml avail gcc`` and ``$ ml avail intel``
- The safest way is to transfer the source code to Bianca via the wharf.
- [Guide for compiling serial and parallel programs](https://www.uppmax.uu.se/support/user-guides/mpi-and-openmp-user-guide/)
- Available combinations of compilers and parallel libraries (openmpi): <https://hackmd.io/_IqCbOiyS8SZ0Uqpa3UpHg?view>



!!! abstract "Keypoints"
    - You have got an overview of the procedures to install packages/libraries and tools on Bianca through the wharf
    - If you feel uncomfortable or think that many users would benefit from the software, ask the support to install it.
