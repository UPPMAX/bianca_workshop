# Own development and Git on Bianca

## Git on Bianca

[More on Git on Bianca](https://www.uppmax.uu.se/support/faq/software-faq/git-on-bianca/)


## Install software yourself

- If not available on Bianca already (like Conda repositories) you may have to use the Wharf to install your tools
    - Alternatively let an Application Expert install the tool as a module.

!!! info Typical workflow for installation

    - Download the source code or binary (Linux on x86 and 64-bit) to Rackham first
    - Transfer to the wharf
    - Then, either 
        - You can install in your home directory.
            - This is handy for personal needs, low numbers of files (i.e. not Conda).
         - Usually better to install in project directory.
            - This way the project contains both data and software — good for reproducibility, collaboration, and everyone's general sanity.
    - Binaries for Linux on x86 and 64-bit should be able to be run directly as it is, see the software specific installation documentation.
    - or build from source, see next session.
     

### Build from source
- To build from source use a compiler module
- We have several compiler versions from GNU and INTEL
- check with: ``$ ml avail gcc`` and ``$ ml avail intel``
- [Guide for compiling serial and parallel programs](https://www.uppmax.uu.se/support/user-guides/mpi-and-openmp-user-guide/)
- Available combinations of compilers and parallel libraries (openmpi): <https://hackmd.io/_IqCbOiyS8SZ0Uqpa3UpHg?view>




