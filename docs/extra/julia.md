# Using Julia packages on Bianca

## UPPMAX Central library

!!! info

    The Julia application at UPPMAX comes with several preinstalled packages.
    A selection of the Julia packages and libraries installed on UPPMAX are:

    - BenchmarkTools
    - CSV
    - CUDA
    - MPI (needs gcc and openmpi loaded)
    - Distributed
    - IJulia (needs python loaded)
    - Plots
    - PyPlot (needs python loaded)
    - Gadfly
    - DataFrames
    - DistributedArrays
    - PlotlyJS

- In Julia you can check if the packages are their by: ``using <package-name>``
    - The first time Julia will precompile the package for you!
- You may control the present "central library" by typing ``ml help julia/<version>`` in the BASH shell.
    - There you will also find which python, gcc and openmpi version that are compatible.
    - Or see the [Julia user guide at UPPMAX](http://docs.uppmax.uu.se/software/julia/){:target="_blank"}
- A possibly more up-to-date status can be found from the Julia shell:

``` julia

    using Pkg
    Pkg.activate(DEPOT_PATH[2]*"/environments/v1.8");     #change version (1.8) accordingly if you have another main version of Julia
    Pkg.status()
    Pkg.activate(DEPOT_PATH[1]*"/environments/v1.8");     #to return to user library

```

### First time run on Bianca

You may have to build some package the first time with ``Pkg.build(“<package>”)``. Since “<package>" is pre-installed centrally on UPPMAX you must activate the central environment by following these steps below. This should only be needed the first time like this

``` julia
> using Pkg
> Pkg.activate(DEPOT_PATH[2]*"/environments/v1.8");
> Pkg.build("<package>")
```

## Install yourself

If you have started Julia once you will get the folders like this in the ``~/.julia`` folder.

```bash
   $ tree .julia/ -d -L 1
   .
   ├── artifacts
   ├── bin
   ├── compiled
   ├── conda
   ├── environments
   ├── logs
   ├── packages
   ├── prefs
   ├── registries
   └── scratchspaces
```

**The plan is that what you install on Rackham should be moved here in the same manner**

- Make an installation of the package on Rackham in the Julia package manager
- Use a transfer method to move the package files to the ``wharf``
    - To be certain to include all files, you may transfer the whole ``.julia`` dir. However, that can grow rather big with time.

**Transfer to the ``wharf``**

``` bash
sftp sftp bjornc-sens2023531@bianca-sftp
sftp> cd bjornc-sens2023531/
sftp> dir
sftp> 
```

If you have not uploaded anything to your ``wharf``, this will be empty. It might have a few things in it.

- **Alt1: If you would like all your locally installed packages:**

``` bash
sftp> put -r ~/.julia
```

- That may take a while since you are sending thousands of files. Consider archiving the .julia directory (a copy) and then tranfer it with the ``put`` command.
    
- **Alt 2: Just transfer the latest installed julia package(s)**

    :warning: The following part still needs some update    

- Check what was installed. It may have been several dependency packages as well. Look at the times!

``` bash
sftp>  lls -lrt ~/.julia/packages
```

``` bash
sftp> put -r ~/.julia/packages/<package name 1>
# and if several packages
sftp> put -r ~/.julia/packages/<package name 2>
# and so on...
```

**Move to site-packages folder**
On Bianca

``` bash
cd /proj/sens2023531/nobackup/wharf/bjornc/bjornc-sens2023531/
mv –a  <file(s)> ~/.julia/packages/
```

!!! info "If many files or packages"

    you may want to tar before copying to include all possible symbolic links:

    ``` bash
    $ tar cfz <package.tar.gz> <package>
    ```

    and in target directory (wharf_mnt) on Bianca:

    ``` bash
    $ tar xfz <tarfile.tar.gz> #if there is a tar file!
    $ mv –a  <file(s)> .julia/packages/
    ```

!!! info "More info"

    - Julia course: isolated environments  <!-- link to https://uppmax.github.io/R-python-julia-matlab-HPC/julia/isolatedJulia.html is broken and hence removed -->
