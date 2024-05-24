# Install with pip to Bianca

???+ question "Read through the content below"

???+ question "Try to follow the steps with the package ``numpy==1.22.3`` for ``python/3.11.8`` "

## Check for packages 
   
- from the Python shell with the ``import`` command
- from BASH shell with the 
	
- ``pip list`` command 
- ``ml help python/3.9.5`` at UPPMAX

**Is it not there? Then proceed!**

!!! info
    **Methods:**
    
    - You can either just download a python package, transfer to ``wharf`` and Bianca and install there.
    - Install it on Rackham. Perhaps you need it here as well! Then transfer to ``wharf`` and Bianca local python library.
    - Make a virtual environment with one or several packages on Rackham. Then transfer to ``wharf`` and Bianca (any place).



## Only download on Rackham and install on Bianca

**Rackham**
``` sh 
$ pip download <package-name>
``` 

**Transfer to the ``wharf``**

``` bash
sftp douglas-sens2017625@bianca-sftp
sftp> cd douglas-sens2017625/
sftp> dir
sftp>
```
If you have not uploaded anything to your ``wharf``, this will be empty. It might have a few things in it.

Now, upload to the ``wharf`` the package <package-name> and all the dependency packages ``pip download`` got you.

``` bash
sftp> put -r <package-name>
```

**Install on Bianca**

On Bianca
install it (Yes, you can do it from this place) by telling pip where to look for packages and dependencies

``` sh 
$ ml python
$ pip install --user --no-index --find-links <path-to-packages> <package-name>
```
<path-to-packages> is where your packages are, if in present working directory it is ``.``

**Then the package ends up in ``~/.local/lib/python<version>/site-packages/`` .**

## Install on Rackham and then transfer to Bianca

!!! info

    **The package ends up on Rackham in ``~/.local/lib/python<version>/site-packages/`` .**

    - Note that python<version> is omitting the last number (bug fix), like 3.8 for python-3.8.7.
    

**Install on Rackham**

``` sh 
$ ml python/<version>		# this is to make use the correct python version and possible dependencies already available
$ pip install --user <package-name>
```
- If there is a requirements.txt file with the content of packages to be installed:

```bash
pip install --user -r requirements.txt
```

**Then the package(s) ends up in ``~/.local/lib/python<version>/site-packages/`` .**

**Transfer to the ``wharf``**

``` bash
sftp douglas-sens2017625@bianca-sftp
sftp> cd douglas-sens2017625/
sftp> dir
sftp>
```
If you have not uploaded anything to your ``wharf``, this will be empty. It might have a few things in it.

- **Alt1: If you would like all your locally installed packages:**

``` bash
sftp> put -r ~/.local/lib/python<version>/site-packages/
```

- **Alt 2: Just transfer the latest installed python package(s)**

- Check what was installed. It may have been several dependency packages as well. Look at the times!

``` bash
sftp>  lls -lrt ~/.local/lib/python<version>/site-packages/
```

``` bash
sftp> put -r ~/.local/lib/python<version>/site-packages/<package name 1>
# and if several packages
sftp> put -r ~/.local/lib/python<version>/site-packages/<package name 2>
# and so on...
```

**Move to site-packages folder**
On Bianca

``` bash
cd /proj/sens2023531/nobackup/wharf/bjornc/bjornc-sens2023531/
mv –a  <file(s)> ~/.local/lib/python<version>/site-packages/
```

!!! info "If many files or packages"

    you may want to tar before copying to include all possible symbolic links:

    ``` bash
    $ tar cfz <tarfile.tar.gz> <package>
    ```
    and in target directory (wharf_mnt) on Bianca:

    ``` bash 
    $ tar xfz <tarfile.tar.gz> #if there is a tar file!
    $ mv –a  <file(s)> ~/.local/lib/python<version>/site-packages/
    ```


## Isolated/virtual environments

- We HIGHLY recommend using a virtual environment during installation, since this makes it easier to install for different versions of Python.  

!!! note
   
    Isolated environments solve a couple of problems:
   
    - You can install specific package, also older, versions into them.
    - You can create one for each project and no problem if the two projects require different versions.
    - You can remove the environment and create a new one, if not needed or with errors.

- More information about [isolated environments](https://uppmax.github.io/HPC-python/isolated_deeper.html). 

**Example, where python packages from the loaded module are used (``--system-site-packages``)**

``` bash
$ module load python/3.6.8
$ python -m venv --system-site-packages <path>/projectB
```

“projectB” is the name of the virtual environment. The directory “projectB” is created in the present working directory. The ``-m`` flag makes sure that you use the libraries from the python version you are using.	

- Activate and install with pip (package one by one or from requirements.txt)

``` bash
$ source <path>/projectB/bin/activate
```
- Note that your prompt is changing to start with (analysis) to show that you are within an environment.
- Install the packages from the file::

```
$ pip install -r requirements.txt

$ pip list   # check
$ deactivate
```

- Virtual environments can be saved easily anywhere
	
**Transfer to the ``wharf``**

``` bash
sftp douglas-sens2017625@bianca-sftp
sftp> cd douglas-sens2017625/
sftp> dir
sftp>
```
If you have not uploaded anything to your ``wharf``, this will be empty. It might have a few things in it.

``` bash
sftp> put -r <path>/projectB
```

**Move to site-packages folder**
On Bianca

``` bash
cd /proj/sens2023531/nobackup/wharf/bjornc/bjornc-sens2023531/
mv –a  projectB <path to any place, like project folder>
```

<br><br><br>

!!! error

    If problems arise, send an email to support@uppmax.uu.se and we'll help you.
    

