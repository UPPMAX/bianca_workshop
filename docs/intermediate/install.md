# Software and package installation on Bianca


!!! info "Learning Objectives"

    Learners
    
    - understand _principles_ how to install software and packages yourself
    - can install Python packages using conda
    - can install Python packages using pip
    - can install R packages
    - understand what containers are
    - can install software using a container
    - can build software from source

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

    See [the basic Bianca course page 'Using the module system'](../modules.md).

- **Lots of programs and tools installed as modules** on Bianca.
    - You can **request installations** but that **may take several days or weeks** to be handled by the application experts at UPPMAX.
    - Workflow: Application expert installs on Rackham and it is synced over to Bianca within a day.
- **Installations yourself**.
    - Workflow: use the **``wharf`` to transfer installation files and packages to Bianca from Rackham or other place**. Here the use of Rackham comes handy because
        - of the internet connection.
        - the computer architecture is somewhat similar such that precompiled binaries or compiled programs (x86_64) on Rackham will most often work also on Bianca.

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

- If not available on Bianca already (like Conda repositories) --> use the ``wharf`` to install your tools

!!! note "Typical workflow for installation"

    - **Download** the
        - _source_ code or
        - _binary_ (Linux on x86 and 64-bit)
    - **Transfer** to the ``wharf``
    - Move file(s) to either
        - ``$HOME`` directory.
            - Handy for personal needs and low numbers of files — i.e. not Conda.
            - Example python/R/julia packages.
        - Usually better to install in project directory.
            - This way the project contains both data and software
            - Good for reproducibility, collaboration, and everyone's general sanity.
    - Then, either:
        - Binaries for Linux on x86 and 64-bit should be able to be run directly as they are.
        - Install program following instructions from documentation of the software.

### Build from source (C/C++ and Fortran)

- To build from source use a **compiler module**
- We have several compiler versions from GNU and INTEL
- ``make`` is installed on the system
    - :warning: It could happen that the "Makefile" contains web fetching, which will not work from Bianca.
    - Usually it is not a problem to build on Rackham and move to Bianca.
- ``cmake`` is available as module

 
!!! info "More info"

    - [Extra material: Build from source](../extra/source_install.md)
    - [Singularity course](https://www.uu.se/centrum/uppmax/utbildning/kurser-och-workshops/basic-singularity){:target="_blank"}
    

### Containers

!!! info

    - Containers let you install programs without needing to think about the computer environment, like
        - operative system
        - dependencies (libraries and other programs) with correct versions
    - 2(3) types
        - Singularity/Apptainer
        - Docker that does not work on HPC-systems
            - But docker images can be used by Singularity and Apptainer
    - Everything is included
    - Draw-backs
        - you install also things that may be already installed
        - therefore, probably more disk space is needed

!!! info "More info"

    - [Extra material: Containers](https://uppmax.github.io/bianca_workshop/extra/containers/){:target="_blank"}
    - [Singularity course](https://www.uu.se/centrum/uppmax/utbildning/kurser-och-workshops/basic-singularity){:target="_blank"}
    
## Packages and libraries to scripting programs

- Python, R and Julia all have some **centrally installed packages** that are available from the modules.
- R has a special module called ``R_packages``, and some Machine Learning python packages are included in the ``python_ml_packages`` module.
- If not found there you can try to install those by yourself.

!!! info

    - "Install packages or not? Check it!"

## Check and install packages

### R

- On UPPMAX the module ``R_packages`` is an omnibus package library containing almost all packages in the CRAN and BioConductor Chances are good the R packages you need are already available once you load this module.  You can quickly check by loading it:

!!! info "Check if available"

    ``$ ml R_packages/4.3.1``

    Then within R, try loading the package you want, like ``glmnet``:

    ``library(glmnet)``

**Is it not there? Then proceed!**

!!! info "Installation principle"

    - install on Rackham or other place
    - sync to ``wharf``
    - move the files on Bianca to local Bianca R package path

!!! info "More info"

    - [Extra material: Installing R packages](https://uppmax.github.io/bianca_workshop/extra/rpackages/)
    - [From R course: packages](https://uppmax.github.io/R-python-julia-matlab-HPC/r/packagesR.html){:target="_blank"}
    - [From R course: isolated environments](https://uppmax.github.io/R-python-julia-matlab-HPC/r/isolatedR.html){:target="_blank"}

### Python and other

!!! info "Check if available"

    - Check **python versions**: ``ml avail python``
    - Check **python packages/modules**

        1. **help** output from: ``ml help python/3.9.5`` at UPPMAX
        2. In a loaded python

            - **Load** a python version, like: ``ml python/3.11.8``
            - from **Python shell** with the ``import`` command
            - from **BASH shell** with the ``pip list`` command

**Is it not there? Then proceed!**

!!! info "Tip Python packages"

    - Try Conda first directly on Bianca. We have mirrored all _major_ Conda repositories directly on UPPMAX, on both Rackham and Bianca. These are updated every third day.
    - If you want to keep number of files down, use PyPI (pip), but then you need to use Rackham and the ``wharf``.

### Conda

- We have mirrored all major Conda repositories directly on UPPMAX, on both Rackham and Bianca. These are updated every third day.

!!! info "More info"

    - [Extra material: Installing Conda packages](https://uppmax.github.io/bianca_workshop/extra/conda/)
    - [Conda user guide](http://docs.uppmax.uu.se/software/conda/)

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
    - [From Python course: packages](https://uppmax.github.io/R-python-julia-matlab-HPC/python/packages.html){:target="_blank"}
    - [From Python course: isolated environments](https://uppmax.github.io/R-python-julia-matlab-HPC/python/isolated.html){:target="_blank"}

### Julia packages

- At UPPMAX there is a central library with installed packages.
- This is good, especially when working on Bianca, since you then do not need to install via the ``wharf``.

!!! info "Installation principle"

    - install on Rackham
    - sync to ``wharf``
    - move the files on Bianca

!!! info "More info"

    - [Extra material: Installing Julia packages](https://uppmax.github.io/bianca_workshop/extra/julia/){:target="_blank"}
    - [UPPMAX julia user guide: Pip](http://docs.uppmax.uu.se/software/julia/){:target="_blank"}
    - [Julia course: isolated environments](https://uppmax.github.io/R-python-julia-matlab-HPC/julia/isolatedJulia.html){:target="_blank"}
    - :warning: contact [NAISS support](javascript:void(window.open('https://supr.naiss.se/support/?centre_resource=c4%27,%27_blank%27,%27toolbar=1,location=1,status=1,menubar=1,scrollbars=1,resizable=1%27));) for individual help!

## Exercise 20 min

???+ question "Pick one of the following topics!"

    - Read the introduction with a demo and use it to solve the exercise in the end.
    - "Containers" contains less material but may take time to install.

    - [Extra material: Installing Conda packages](https://uppmax.github.io/bianca_workshop/extra/conda/){:target="_blank"}
    - [Extra material: Installing pip packages](https://uppmax.github.io/bianca_workshop/extra/pip/){:target="_blank"}
    - [Extra material: Installing R packages](https://uppmax.github.io/bianca_workshop/extra/rpackages/){:target="_blank"}
    - [Extra material: Containers](https://uppmax.github.io/bianca_workshop/extra/containers/){:target="_blank"}

    - One breakout room per topic: Help each-other!

!!! discussion

    - Did it work out well?
    - Any questions?
    - Any input?

## Own development and Git

- [Own development and git](https://uppmax.github.io/bianca_workshop/extra/devel/)

!!! abstract "Keypoints"

    - You have got an overview of the procedures to install packages/libraries and tools on Bianca through the ``wharf``
    - If you feel uncomfortable or think that many users would benefit from the software, ask the support to install it.
