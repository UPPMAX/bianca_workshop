# Build from source

- To build from source use a **compiler module**
- We have several compiler versions from GNU and INTEL. Check what is available with:
    - ``$ ml avail gcc``
    - or for intel
        - ``$ ml avail intel`` or
        - ``module load intel-oneapi``
        - ``module avail compiler``
- ``make`` is installed on the system
    - :warning: It could happen that the "Makefile" contains web fetching, which will not work from Bianca.
    - Usually it is not a problem to build on Rackham and move to Bianca.
- ``cmake`` is available as module
    - Check with: ``$ ml avail cmake``
- [Guide for compiling **serial** programs](https://docs.uppmax.uu.se/software/compiling_serial/){:target="_blank"}
- [Guide for compiling **parallel** programs](https://docs.uppmax.uu.se/software/compiling_parallel/){:target="_blank"}

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
