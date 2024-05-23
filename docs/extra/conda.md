# Conda on Bianca

!!! info "Objectives"
    - This is a brief description of the necessary steps to use the local Conda repository at UPPMAX, and install things for yourself or your project using Conda. 

!!! info "Hint"

    - On Bianca, Conda is the first choice when installing packages, because there is a local mirror of most of the Conda repositories.

## Using Conda

!!! info "Conda cheat sheet"

    -   List all environments: `conda info -e` or `conda env list`

    -   Create a conda environment (it is good to directly define the packages included AND channels do not need to be explicitly mentioned)
    
        ```conda create --prefix /some/path/to/env <package1> [<package2> ... ] ```
       
        - On our systems the above should replace `conda create --name myenvironment ...`
       
    -   Create a new environment from requirements.txt:
   
        - `conda create --prefix /some/path/to/env --file requirements.txt`

    -   Activate a specific environment: `conda activate myenvironment`

    -   List packages in present environment: `conda list`

        - Also pip list will work

    -   Install additional package from an active environment: 
    
        - `conda install somepackage`

    -   Install from certain channel (conda-forge):
       
        - `conda install -c conda-forge somepackage`

    -   Install a specific version: `conda install somepackage=1.2.3`

        -   Install a specific version: `conda install somepackage=1.2.3`

    -   Deactivate current environment: `conda deactivate`


## Your conda settings on Bianca

- ```export CONDA_ENVS_PATH=/a/path/to/a/place/in/your/project/```

!!! tip

    - You may wnat to have the same path for all conda environments in the present project
    - ``echo "export CONDA_ENVS_PATH=/a/path/to/a/place/in/your/project/" >> ~/.bashrc`` 
        - Example: ``echo "export CONDA_ENVS_PATH=/proj/sens2023598/bjornc/conda" >> ~/.bashrc``

!!! warning

    - It seems you are required to use this path, ending with the name of your environment, together with ``--prefix`` when you install new envronments AND packages also after activating the conda environment!
    Like: ``conda install --prefix $CONDA_ENVS_PATH/<your-environment> ...``

!!! tip

    -   REMEMBER TO `conda clean -a` once in a while to remove unused and unnecessary files

??? info "By choice"

    - Run `source conda_init.sh` to initialise your shell (bash) to be able to run `conda activate` and `conda deactivate` etcetera instead of `source activate`. It will modify (append) your `.bashrc` file.

    -   When conda is loaded you will by default be in the ``base`` environment, which works in the same way as other Conda environments. It is a “best practice” to avoid installing additional packages into your base software environment unless it is very general packages

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

- You reach them all by loading the ``conda`` module.
- You don't have to state the specific channel when using UPPMAX.
- Also, you are offline on Bianca which means that the default is `--offline`, which you can specify if you want to simulate the experience on Rackham.

!!! tip

    - If you need a channel that isn't in our repository, we can easily add it. Just send us a message and we will do it.


## Make a new conda environment

!!! tip

    - Since python or other packages are dependent on each-other expect solving the versions takes some time.
    - use an interactive session!


1.  Load Conda

    ```bash 
    module load conda
    ```
    
    - This grants you access to the latest version of Conda and all major repositories on all UPPMAX systems.
    - Check the text output as ``conda`` is loaded, especially the first time, see below


3.  Create the Conda environment

    - Example:

        ```bash
        conda create --prefix  $CONDA_ENVS_PATH/python36-env python=3.6 numpy=1.13.1 matplotlib=2.2.2
        ```

    !!! info "The `mamba` alternative is not needed in newer versions of Conda!

    - It all worked if you get something like this:
  
        ```bash
        # To activate this environment, use
        #
        #     $ conda activate python36-env 
        #
        # To deactivate an active environment, use
        #
        #     $ conda deactivate
        ```

    
4.  Activate the conda environment by `source activate` if you have not enabled ``conda activate``, see above:

    ```bash 
    source activate python36-env
    ```

    - You will see that your prompt is changing to start with `(python-36-env)` to show that you are within an environment.

    - You can also see the installed packages by:
  
    ```bash
    conda list
    pip list
    ```
    - you can also add more packages within the environment by exact version (use `=`) or latest (?) compatible version:
    ```bash
    conda install --prefix   $CONDA_ENVS_PATH/python36-env pandas
    ```
    - that may have given you ``pandas=1.1.5`` which would be the newest version compatible with ``python3.6`` and ``numpy=1.13.1``

5.  Now do your work!

6.  Deactivate with ``conda deactivate`` (this will work in any case!)

    ```bash
    (python-36-env) $ conda deactivate
    ```

!!! warning

    -   Conda is known to create **many** *small* files. Your diskspace is not only limited in GB:s, but also in number of files (typically `300000` in ``$home``).
    -   Check your disk usage and quota limit with `uquota`
    -   Do a `conda clean -a` once in a while to remove unused and unnecessary files


-   [More info about Conda on UPPMAX](http://docs.uppmax.uu.se/cluster_guides/conda/)

## Working with Conda environments defined by files

-   Create an environment based on dependencies given in an environment
    file:

        $ conda env create --file environment.yml

-   Create file from present conda environment:

        $ conda env export > environment.yml

`environments.yml` (for conda) is a yaml-file which looks like this:

```yaml
name: my-environment
channels:        # not needed on bianca
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
channels:            #not needed on bianca
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


???+ question "Create a conda environment and install some packages"

    -   First check the current installed packages while having `python/3.9.5` loaded

    -   Open a new terminal and have the old one available for later comparison

    -   Unload ``python`` module

    -   Use the ``conda`` module on Rackham and create an environment with name `bianca-course` with `python 3.7` and `numpy 1.15`

    -   Use your a path for `CONDA_ENVS_PATH` of your own choice (not doing this is perfectly OK and isnatalls in your ``$HOME`` folder) or `/proj/sens2023531/<user>/` 
        
        -   (It may take a minute or so)

    -   Activate!

    -   Check with `pip list` what is there. Compare with the environment given from the python module in the first terminal window.

    -   Which version of Python did you get?

    -   Don't forget to deactivate the Conda environment before doing other exercises!


??? Solution for UPPMAX

    Write this in the terminal

    ``` sh
    $ module load conda
    ($ export CONDA_ENVS_PATH=/proj/sens2023598/$USER)
    $ conda create --name HPC-python23 python=3.7 numpy=1.15
    $ source activate HPC-python23
    $ pip list
    $ python -V
    $ conda deactivate
    ```





!!! abstract "keypoints"
    
    -   Conda is an installer of packages but also bigger toolkits

    - Conda on Bianca is easy since the repos in the most used channels are local.

    -   Conda creates isolated environments not clashing with other installations of python and other versions of packages

    -   Conda environment requires that you install all packages needed by yourself, although automatically.
    
    -   That is, you cannot load the python module and use the packages therein inside your Conda environment.

