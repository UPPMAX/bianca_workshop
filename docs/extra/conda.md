# Conda on Bianca

???+ queastion "Read [Conda user guide](http://docs.uppmax.uu.se/cluster_guides/conda/)"

    - Read [Conda user guide](http://docs.uppmax.uu.se/cluster_guides/conda/)
        - skip [Working with Conda environments defined by files](http://docs.uppmax.uu.se/cluster_guides/conda/#working-with-conda-environments-defined-by-files)
    
    - ONLY for the interested: [Working with Conda environments defined by files](http://docs.uppmax.uu.se/cluster_guides/conda/#working-with-conda-environments-defined-by-files)
        - On bianca you have to get the `environments.yml` to wharf first!


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




# Exercises

!!! tip

    - You may want to have the same path for all conda environments in the present project
    - ``echo "export CONDA_ENVS_PATH=/a/path/to/a/place/in/your/project/" >> ~/.bashrc`` 
        - Example: ``echo "export CONDA_ENVS_PATH=/proj/sens2023598/bjornc/conda" >> ~/.bashrc``

!!! warning

    - It seems you are required to use this path, ending with the name of your environment, together with ``--prefix`` when you install new envronments AND packages also after activating the conda environment!
    Like: ``conda install --prefix $CONDA_ENVS_PATH/<your-environment> ...``


???+ question "Create a conda environment and install some packages"

    -   First check the current installed packages while having `python/3.9.5` loaded

    -   Open a new terminal and have the old one available for later comparison

    -   Make sure ``python`` module is not active in the new terminal

    -   Use the ``conda`` module on Bianca and the ``conda-forge`` channel to create an environment with name `bianca-course` with `python 3.7` and `numpy 1.15`

    -   Use your a path for `CONDA_ENVS_PATH` of your own choice (not doing this is perfectly OK and installs in your ``$HOME`` folder), that is `/proj/sens2023531/<user>/` 
        
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

    - It should show numpy=1.15 among others and the python version 3.7.X





!!! abstract "keypoints"
    
    - Conda is an installer of packages but also bigger toolkits

    - Conda on Bianca is easy since the repos in the most used channels are local.

    -  Conda creates isolated environments not clashing with other installations of python and other versions of packages

    -   Conda environment requires that you install all packages needed by yourself, although automatically.
    
    -   That is, you _cannot_ load the python module and use the packages therein inside your Conda environment.

