# Running Jupyter on Bianca

!!! info

    - You can run Python in a **Jupyter-notebook**, i.e. in a web interface with possibility of inline figures and debugging.
    - **Jupyter-lab** is installed in the **python>=3.10.8 module**
    - You can install a personal version of juputer-lab with Conda for lower versions. An easy way to do this is to load the python module as well. In shell:

!!! warning

    Always start Jupyter in a **ThinLinc** session and preferably in an **interactive** session.


## Start

Start a notebook like this:

```bash
$ module load python/<version>
$ jupyter-notebook
```
or jupyter lab:

``` bash
$ jupyter-lab
```

A local Firefox session (not a internet web page!) should start with the Jupyter notebook/lab interface. If not,  copy-paste one of the addresses into the address files in an open firefox session (start with ``firefox &``).

## Jupyter in a virtual environment (venv)

You could also use jupyter- (lab or notebook) in a virtual environment.

If you decide to use the ``--system-site-packages`` configuration you will get jupyter from the python modules you created you virtual environment with.
However, you won't find your locally installed packages from that jupyter session. To solve this, reinstall jupyter within the virtual environment by force (option ``-I``):

```bash
$ pip install -I jupyter
```
and run it as above.

Be sure to start the kernel with the virtual environment name, like "project A", and not "Python 3 (ipykernel)".
      
