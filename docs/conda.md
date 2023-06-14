# Conda on Bianca

!!! info "Objectives"
    - This is a brief description of the necessary steps to use the local Conda repository at UPPMAX, and install things for yourself or your project using Conda. 

!!! info "Summary of ``conda`` commands"
    - ```module load conda```
    - ```export CONDA_ENVS_PATH=/a/path/to/a/place/in/your/project/```
    - ```conda create``` ... etc
    - Remember to run ```conda clean -a``` once in a while. When you load the module, there is also a reminder displayed, so you get this info there also.


# Conda basics

-   What does Conda do?
-   How to create a Conda environment

-   Learn pros and cons with Conda
-   Learn how to install packages and work with the Conda (isolated)
    environment

!!! info "Hint"

    -   On Bianca Conda is the first choice when installing packages, because there is a local mirror of most of the Conda repositories.

## Using Conda

!!! info "Conda cheat sheet"

    -   List packages in present environment: `conda list`

    -   List all environments: `conda info -e` or `conda env list`

    -   Install a package: `conda install somepackage`

    -   Install from certain channel (conda-forge):
       
        - `conda install -c conda-forge somepackage`

    -   Install a specific version: `conda install somepackage=1.2.3`

    -   Create a new environment: `conda create --name myenvironment`

    -   Create a new environment from requirements.txt:
   
        - `conda create --name myenvironment --file requirements.txt`

    -   On e.g. HPC systems where you don’t have write access to central installation directory: ```conda create --prefix /some/path/to/env```

    -   Activate a specific environment: `conda activate myenvironment`

    -   Deactivate current environment: `conda deactivate`



## Installing using Conda


We have mirrored all major Conda repositories directly on UPPMAX, on
both Rackham and Bianca. These are updated every third day. We have the
following channels available:

!!! info "Conda Channels"

    -   bioconda
    -   biocore
    -   conda-forge
    -   dranew
    -   free
    -   main
    -   pro
    -   qiime2
    -   r
    -   r2018.11
    -   scilifelab-lts
    -   nvidia
    -   pytorch

You reach them all by loading the ``conda`` module. You don't have to state the specific channel when using UPPMAX. Also, you are offline on Bianca which means that the default is `--offline`, which you can specify if you want to simulate the experience on Rackham.

If you need a channel that isn't in our repository, we can easily add it. Just send us a message and we will do it.


## First steps

!!! Tip

    There will be an exercise in the end!


1.  First load our conda module (there is no need to install you own
    miniconda, for instance)

    ``` bash 
        module load conda
    ```
    
    - This grants you access to the latest version of Conda and all major repositories on all UPPMAX systems.
    - Check the text output as ``conda`` is loaded, especially the first time, see below

    !!! info "Conda load output"

        -   The variable CONDA_ENVS_PATH contains the location of your environments. Set it to your project's environments folder if you have one.
        -   Otherwise, the default is ``~/.conda/envs``.
        -   You may run `source conda_init.sh` to initialise your shell to be able to run `conda activate` and `conda deactivate` etc.
        -   Just remember that this command adds stuff to your shell outside the scope of the module system.
        -   REMEMBER TO `conda clean -a` once in a while to remove unused and unnecessary files


2.  First time

    -   The variable CONDA_ENVS_PATH contains the location of your environments. Set it to your project's environments folder if you have one.
    -   Otherwise, the default is ``~/.conda/envs``.
    -   Example:

    ``` bash 
        export CONDA_ENVS_PATH=/proj/\<your-project-id\>/nobackup/\<username\>
    ```
    
    ??? info "By choice"

        Run `source conda_init.sh` to initialise your shell (bash) to be able to run `conda activate` and `conda deactivate` etcetera instead of `source activate`. It will modify (append) your `.bashrc` file.

    -   When conda is loaded you will by default be in the base environment, which works in the same way as other Conda environments. include a Python installation and some core system libraries and dependencies of Conda. It is a “best practice” to avoid installing additional packages into your base software environment.

3.  Create the Conda environment

    -   Example:

    ```bash
      conda create --name python36-env python=3.6 numpy=1.13.1
      matplotlib=2.2.2
    ```

    !!! info "The `mamba` alternative"
        -   `mamba` is a fast drop-in alternative to conda, using
            "libsolv" for dependency resolution. It is available from the
            `conda` module.

        -   Example:

            ```bash
            mamba create --name python37-env python=3.7 numpy=1.13.1
            matplotlib=2.2.2
            ```
        
4.  Activate the conda environment by:

    ```bash 
    source activate python36-env
    ```
    -   You will see that your prompt is changing to start with `(python-36-env)` to show that you are within an environment.

5.  Now do your work!

6.  Deactivate

    ```bash
    (python-36-env) $ conda deactivate
    ```

!!! warning

    -   Conda is known to create **many** *small* files. Your diskspace is not only limited in GB:s, but also in number of files (typically `300000` in ``$home``).
    -   Check your disk usage and quota limit with `uquota`
    -   Do a `conda clean -a` once in a while to remove unused and unnecessary files


-   [More info about Conda on UPPMAX](https://uppmax.uu.se/support/user-guides/conda-user-guide/)

## Working with Conda environments defined by files

-   Create an environment based on dependencies given in an environment
    file:

        $ conda env create --file environment.yml

-   Create file from present conda environment:

        $ conda env export > environment.yml

`environments.yml` (for conda) is a yaml-file which looks like this:

``` yaml
name: my-environment
channels:
- defaults
dependencies:
- numpy
- matplotlib
- pandas
- scipy
```

`environments.yml` with versions:

``` yaml
name: my-environment
channels:
- defaults
dependencies:
- python=3.7
- numpy=1.18.1
- matplotlib=3.1.3
- pandas=1.1.2
- scipy=1.6.2
```

!!! admonition "More on dependencies"

    - Dependency management from course [Python for Scientific computing](https://aaltoscicomp.github.io/python-for-scicomp/dependencies/)

# Exercises


???+ Create a conda environment and install some packages

    -   First check the current installed packages while having `python/3.9.5` loaded

    -   Open a new terminal and have the old one available for later comparison

    -   Unload ``python`` module

    -   Use the ``conda`` module on Rackham and create an environment with name `HPC-python23` with `python 3.7` and `numpy 1.15`

    -   Use your a path for `CONDA_ENVS_PATH` of your own choice or `/proj/py-r-jl/<user>/python` 
        
        -   (It may take a minute or so)

    -   Activate!

    -   Check with `pip list` what is there. Compare with the environment given from the python module in the first terminal window.

    -   Which version of Python did you get?

    -   Don't forget to deactivate the Conda environment before doing other exercises!


??? Solution for UPPMAX

    Write this in the terminal

    ``` sh
    $ module load conda
    $ export CONDA_ENVS_PATH=/proj/py-r-jl/<user>/python
    $ conda create --name HPC-python23 python=3.7 numpy=1.15
    $ source activate HPC-python23
    $ pip list
    $ python -V
    $ source deactivate
    ```





!!! abstract "keypoints"
    
    -   Conda is an installer of packages but also bigger toolkits

    -   Conda creates isolated environments (see next section) not clashing with other installations of python and other versions of packages

    -   Conda environment requires that you install all packages needed by yourself.
    
    -   That is, you cannot load the python module and use the  packages therein inside you Conda environment.

