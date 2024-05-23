# Software and package installation on Bianca


!!! info "Objectives" 

    - Understand _principles_ how to install software and packages yourself
    - Understand what containers are
    	- understand difference of Docker and Apptainer/Singularity
    - Know by doing how to install software, packages and libraries from at least one of the following
        - ... using Conda
        - ... using Python packages with pip
        - ... using R packages
        - ... using Apptainer/Singularity
    - Understand how to build from source


???- info "Notes for teachers"

    Teaching goals:

    - The learners have explored the UPPMAX documentation
    - The learners have installed a package (R, python or julia)
    - The learners understand how to install own software

    Lesson plan:

    ```mermaid
    gantt
      title IDEs
      dateFormat X
      axisFormat %m
      Introduction: intro, 0, 10m
      Vote on wih IDE: vote, after intro, 5m
      Exercise with personal favorite package: crit, exercise, after vote, 20m
      Feedback: feedback, after exercise, 10m
      Installing software and developing: monologue, after feedback, 5m
      Break: milestone, after monologue
    ```

## The module system

???- tip "Forgot how to use the module system?"

    See [the basic Bianca course page 'Using the module system'](../modules1.md).

- There is a **lot of programs and tools installed as modules** on Bianca.
- These have typically been **installed on Rackham** and is **synced over to Bianca a couple of times per day**.
- **You can request installations** but that may take **several days or weeks** to be handled by the application experts at UPPMAX.
- But you may be able to do **installations yourself**. Here the use of Rackham comes handy because of the:
    - internet connection
    - the computer architecture is somewhat similar such that precompiled binaries or compiled programs (x86_64) on Rackham will most often work also on Bianca.
    - you can use the **``wharf`` to transfer installation files and packages to Bianca from Rackham**

## Content

- Principles of software installation on Bianca
    - From source
    - From binary
    - Containers
- Principles of packages on Bianca
- Exercise: test yourself in EITHER
    - R
    - Conda
    - Python/pip
    - Julia
- Feedback
- Development and Git on Bianca

## Install software yourself

- If not available on Bianca already (like Conda repositories) you may have to use the ``wharf`` to install your tools
- Alternatively let an Application Expert install the tool as a module.

!!! note "Typical workflow for installation"

    - **Download** the 
        - _source_ code or 
	- _binary_ (Linux on x86 and 64-bit) to Rackham first
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
- [Guide for compiling **serial** programs](http://docs.uppmax.uu.se/cluster_guides/compiling_serial/){:target="_blank"}
- [Guide for compiling **parallel** programs](http://docs.uppmax.uu.se/cluster_guides/compiling_parallel/){:target="_blank"}
    - [Available **combinations** of compilers and parallel libraries](http://docs.uppmax.uu.se/cluster_guides/compiling_parallel/#mpi-using-the-openmpi-library){:target="_blank"}

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


### Containers

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
    - [Singularity course](https://www.uu.se/centrum/uppmax/utbildning/kurser-och-workshops/basic-singularity)
    
## Packages and libraries to scripting programs

- Python, R and Julia all have some **centrally installed packages** that are available from the modules. 
- R has a special module called ``R_packages``, and some Machine Learning python packages are included in the ``python_ml_packages`` module.
- If not found there you can try to install those by yourself.

!!! info 

    - "Install packages or not? Check it!"

## Check and install packages! 

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

!!! info "Installation principle"

    - install on Rackham
    - sync to ``wharf``
    - move the files on Bianca

!!! info "More info"

    - [Extra material: Installing R packages](https://uppmax.github.io/bianca_workshop/extra/rpackages/)
    - [From R course: packages](https://uppmax.github.io/R-python-julia-HPC/r/packagesR.html){:target="_blank"}
    - [From R course: isolated environments](https://uppmax.github.io/R-python-julia-HPC/R/isolatedR.html){:target="_blank"}


### Python and other

- **Check** python versions: ``ml avail python``
- Check the **help** output from: ``ml help python/3.9.5`` at UPPMAX
- **Load** a python version, like: ``ml python/3.10.8``
- from the **Python shell** with the ``import`` command
- from **BASH shell** with the ``pip list`` command 

**Is it not there? Then proceed!**

!!! info "Tip Python packages"

    - Try Conda first directly on Bianca. We have mirrored all _major_ Conda repositories directly on UPPMAX, on both Rackham and Bianca. These are updated every third day.
    - If you want to keep number of files down, use PyPI (pip), but then you need to use Rackham and the ``wharf``.


### Conda

- We have mirrored all major Conda repositories directly on UPPMAX, on both Rackham and Bianca. These are updated every third day.

!!! info "More info"

    - [Extra material: Installing Conda packages](https://uppmax.github.io/bianca_workshop/extra/conda/)
    - [Conda user guide](http://docs.uppmax.uu.se/cluster_guides/conda/)

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
    - [UPPMAX Python user guide: Pip](http://docs.uppmax.uu.se/software/python_install_packages/#pip){:target="_blank"}
    - [From Python course: packages](https://uppmax.github.io/R-python-julia-HPC/python/packages.html){:target="_blank"}
    - [From Python course: isolated environments](https://uppmax.github.io/R-python-julia-HPC/python/isolated.html){:target="_blank"}


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
    - [UPPMAX julia user guide: Pip](http://docs.uppmax.uu.se/software/julia/){:target="_blank"}
    - [Julia course: isolated environments](https://uppmax.github.io/R-python-julia-HPC/julia/isolatedJulia.html){:target="_blank"}
    - :warning: contact support@uppmax.uu.se for individual help!



## Exercise 20 min

???+ question "Pick one of the following topics!"

    - Read the introduction with a demo and use it to solve the exercise in the end
    - "Containers" contains less material but may take time to install.

    - [Extra material: Installing Conda packages](https://uppmax.github.io/bianca_workshop/extra/conda/){:target="_blank"}
    - [Extra material: Installing pip packages](https://uppmax.github.io/bianca_workshop/extra/pip/){:target="_blank"}
    - [Extra material: Installing R packages](https://uppmax.github.io/bianca_workshop/extra/rpackages/){:target="_blank"}
    - [Extra material: Containers](https://uppmax.github.io/bianca_workshop/extra/containers/){:target="_blank"}

!!! discussion

    - Did it work out well?
    - Any questions?
    - Any input?


## Own development and Git

- [Own development and git](https://uppmax.github.io/bianca_workshop/extra/devel/)
    

!!! abstract "Keypoints"
    - You have got an overview of the procedures to install packages/libraries and tools on Bianca through the ``wharf``
    - If you feel uncomfortable or think that many users would benefit from the software, ask the support to install it.
