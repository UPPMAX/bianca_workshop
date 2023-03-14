# Conda

!!! info "Objectives"
    - We'll ...

!!! abstract "keypoints"
    - bullet 1
    - bullet 2

This is a brief description of the necessary steps to use the local Conda repository at UPPMAX, and install things for yourself or your project using Conda. 

To use conda on UPPMAX Bianca you could use the following steps:

1. Run: module load conda
This loads the module conda and grants you access to the latest version of conda and all major repositories on all UPPMAX systems.

2. Run: export CONDA_ENVS_PATH=/a/path/to/a/place/in/your/project/
The variable CONDA_ENVS_PATH contains the location of your environments. As long as all members of a project set the variable CONDA_ENVS_PATH to the same place, you can all use the same environments. Set it to your project's environments folder. Otherwise, the default is ~/.conda/envs (warning: installing conda environments in your home directory like this can quickly consume all your space and/or files).

3: Run: conda create ... etc
This where you use conda as you normally would. On Bianca this is equal to running conda with --offline on Rackham since there is no internet. 

Customising
    - You may run 'source conda_init.sh' to initialise your shell to be able to run 'conda activate' and 'conda deactivate' etc. Just remember that this command ADDS STUFF TO YOUR SHELL outside the scope of the module system.

Remember to run 'conda clean -a' once in a while. When you load the module, there is also a reminder displayed, so you get this info there also.




# Install with pip

You use `pip` this way, in a Linux shell OR a python shell:

``` sh 
$ pip install --user <package>
```

Use `pip3` if you loaded python3.

Then the package ends up in
`~/.local/lib/python<version>/site-packages/` .

Note that python\<version\> is omitting the last number (bug fix), like
3.8 for python-3.8.7. We HIGHLY recommend using a virtual environment
during installation, since this makes it easier to install for different
versions of Python. More information will follow later in this course
(<https://uppmax.github.io/HPC-python/isolated.html>).

<div class="note" markdown="1">

<div class="title" markdown="1">

Note

</div>

You will test this in the separated sessions about isolated environments
in a while.

</div>

<div class="keypoints" markdown="1">

-   You can check for packages

> -   from the Python shell with the `import` command
>
> -   from BASH shell with the
>
>     > -   `pip list` command at both centers
>     > -   `ml help python/3.9.5` at UPPMAX

-   Installation of Python packages can be done either with **PYPI** or
    **Conda**
-   You install own packages with the `pip install` command (This is the
    recommended way on HPC2N)
-   At UPPMAX Conda is also available (See Conda section)

</div>

# Conda

<div class="questions" markdown="1">

-   What does Conda do?
-   How to create a Conda environment

</div>

<div class="objectives" markdown="1">

-   Learn pros and cons with Conda
-   Learn how to install packages and work with the Conda (isolated)
    environment

</div>

<div class="hint" markdown="1">

<div class="title" markdown="1">

Hint

</div>

-   On Bianca (with no internet), Conda is the first choice when
    installing packages, because there is a local mirror of most of the
    Conda repositories.

</div>

## Using Conda

<div class="admonition" markdown="1">

Conda cheat sheet

-   List packages in present environment: `conda list`
-   List all environments: `conda info -e` or `conda env list`
-   Install a package: `conda install somepackage`
-   Install from certain channel (conda-forge):
    `conda install -c conda-forge somepackage`
-   Install a specific version: `conda install somepackage=1.2.3`
-   Create a new environment: `conda create --name myenvironment`
-   Create a new environment from requirements.txt:
    `conda create --name myenvironment --file requirements.txt`
-   On e.g. HPC systems where you don’t have write access to central
    installation directory: conda create --prefix /some/path/to/env\`\`
-   Activate a specific environment: `conda activate myenvironment`
-   Deactivate current environment: `conda deactivate`

</div>

## Install with conda (UPPMAX)

<div class="note" markdown="1">

<div class="title" markdown="1">

Note

</div>

We have mirrored all major conda repositories directly on UPPMAX, on
both Rackham and Bianca. These are updated every third day. We have the
following channels available:

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

You reach them all by loading the conda module. You don't have to state
the specific channel when using UPPMAX. Otherwise you do with
`conda -c <channel> ...`

</div>

## First steps

<div class="tip" markdown="1">

<div class="title" markdown="1">

Tip

</div>

There will be an exercise in the end!

</div>

1.  First load our conda module (there is no need to install you own
    miniconda, for instance)

> <div class="prompt" markdown="1">
>
> bash \$
>
> module load conda
>
> </div>
>
> -   This grants you access to the latest version of Conda and all
>     major repositories on all UPPMAX systems.
> -   Check the text output as conda is loaded, especially the first
>     time, see below
>
> > <div class="admonition dropdown" markdown="1">
> >
> > Conda load output
> >
> > -   The variable CONDA_ENVS_PATH contains the location of your
> >     environments. Set it to your project's environments folder if
> >     you have one.
> > -   Otherwise, the default is \~/.conda/envs.
> > -   You may run `source conda_init.sh` to initialise your shell to
> >     be able to run `conda activate` and `conda deactivate` etc.
> > -   Just remember that this command adds stuff to your shell outside
> >     the scope of the module system.
> > -   REMEMBER TO `conda clean -a` once in a while to remove unused
> >     and unnecessary files
> >
> > </div>

1.  First time

> -   The variable CONDA_ENVS_PATH contains the location of your
>     environments. Set it to your project's environments folder if you
>     have one.
>
> -   Otherwise, the default is \~/.conda/envs.
>
> -   Example:
>
>     > <div class="prompt" markdown="1">
>     >
>     > bash \$
>     >
>     > export
>     > CONDA_ENVS_PATH=/proj/\<your-project-id\>/nobackup/\<username\>
>     >
>     > </div>
>
> > <div class="admonition dropdown" markdown="1">
> >
> > By choice
> >
> > Run `source conda_init.sh` to initialise your shell (bash) to be
> > able to run `conda activate` and `conda deactivate` etcetera instead
> > of `source activate`. It will modify (append) your `.bashrc` file.
> >
> > </div>
> >
> > -   When conda is loaded you will by default be in the base
> >     environment, which works in the same way as other conda
> >     environments. include a Python installation and some core system
> >     libraries and dependencies of Conda. It is a “best practice” to
> >     avoid installing additional packages into your base software
> >     environment.

1.  Create the conda environment

> -   Example:
>
>     <div class="prompt" markdown="1">
>
>     bash \$
>
>     conda create --name python36-env python=3.6 numpy=1.13.1
>     matplotlib=2.2.2
>
>     </div>
>
>     <div class="admonition dropdown" markdown="1">
>
>     The `mamba` alternative
>
>     </div>
>
>     -   `mamba` is a fast drop-in alternative to conda, using
>         "libsolv" for dependency resolution. It is available from the
>         `conda` module.
>
>     -   Example:
>
>         > <div class="prompt" markdown="1">
>         >
>         > bash \$
>         >
>         > </div>
>         >
>         > mamba create --name python37-env python=3.7 numpy=1.13.1
>         > matplotlib=2.2.2

1.  Activate the conda environment by:

    > <div class="prompt" markdown="1">
    >
    > bash \$
    >
    > </div>
    >
    > source activate python36-env
    >
    > -   You will see that your prompt is changing to start with
    >     `(python-36-env)` to show that you are within an environment.

2.  Now do your work!

3.  Deactivate

> <div class="prompt" markdown="1" language="bash"
> prompts="(python-36-env) $">
>
> conda deactivate
>
> </div>

<div class="warning" markdown="1">

<div class="title" markdown="1">

Warning

</div>

-   Conda is known to create **many** *small* files. Your diskspace is
    not only limited in GB, but also in number of files (typically
    `300000` in \$home).
-   Check your disk usage and quota limit with `uquota`
-   Do a `conda clean -a` once in a while to remove unused and
    unnecessary files

</div>

-   [More info about Conda on
    UPPMAX](https://uppmax.uu.se/support/user-guides/conda-user-guide/)

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

<div class="admonition" markdown="1">

More on dependencies

-   Dependency management from course [Python for Scientific
    computing](https://aaltoscicomp.github.io/python-for-scicomp/dependencies/)

</div>

# Exercises

<div class="challenge" markdown="1">

UPPMAX: Create a conda environment and install some packages

-   First check the current installed packages while having
    `python/3.9.5` loaded
-   Open a new terminal and have the old one available for later
    comparison
-   Use the conda module on Rackham and create an environment with name
    `HPC-python23` with `python 3.7` and `numpy 1.15`

> -   
>
>     Use your a path for `CONDA_ENVS_PATH` of your own choice or `/proj/py-r-jl/<user>/python`
>
>     :   -   (It may take a minute or so)

-   Activate!
-   Check with `pip list` what is there. Compare with the environment
    given from the python module in the first terminal window.

> -   Which version of Python did you get?

-   Don't forget to deactivate the Conda environment before doing other
    exercises!

</div>

<div class="solution dropdown" markdown="1">

Solution for UPPMAX

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

</div>

<div class="keypoints" markdown="1">

-   Conda is an installer of packages but also bigger toolkits

-   Conda creates isolated environments (see next section) not clashing
    with other installations of python and other versions of packages

-   Conda environment requires that you install all packages needed by
    yourself.

    > -   That is, you cannot load the python module and use the
    >     packages therein inside you Conda environment.

</div>
