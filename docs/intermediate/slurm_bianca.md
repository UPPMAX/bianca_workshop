# Slurm at Bianca

## More Slurm and other advanced UPPMAX techniques

- A closer look at Slurm
- Using the GPUs
- Job efficiency with the ``jobstats`` tool
- Advanced job submission

## The Slurm Workload Manager

- Free, popular, lightweight
- Open source: https://slurm.schedmd.com
- available at all SNIC centra
- [UPPMAX Slurm userguide](http://docs.uppmax.uu.se/cluster_guides/slurm/)

### More on sbatch
Recap:

sbatch | -A sens2023598  |   -t 10:00 | -p core | -n 10 | my_job.sh
-|-|-|-|-|-
slurm batch| project name | max runtime | partition ("job type") | #cores | job script


### More on time limits

- ``-t dd-hh:mm:ss``
- ``0-00:10:00 = 00:10:00 = 10:00 = 10``
- ``0-12:00:00 = 12:00:00``
- ``3-00:00:00 =                    3-0``
- ``3-12:10:15``

### Job walltime

???- question "When you have no idea how long a program will take to run, what should you book?"

    A: very long time, e.g. 10-00:00:00

???- question "When you have an idea of how long a program would take to run, what should you book?"

    A: overbook by 50%

### More on partitions

- ``-p core``
    - “core” is the default partition
    - ≤ 16 cores on Bianca
    - a script or program written without any thought on parallelism will use 1 core

- ``-p node`
    - if you wish to book full node(s)

### Quick testing

- The “devel” partition

  - max 2 nodes per job
  - up to 1 hour in length
  - only 1 at a time
  - ``-p devcore``, ``-p devel`
???- question "Any free nodes in the devel partition? Check status with"

    - ``sinfo -p devel``
    - ``jobinfo -p devel`
   
- more on these tools later
- High priority queue for short jobs

  - 4 nodes
  - up to 15 minutes
  - ``--qos=short``

### Debugging or complicated workflows
- Interactive jobs

  - handy for debugging a code or a script by executing it line by line or for using programs with a graphical user interface
  - ``salloc -n 80 -t 03:00:00 -A sens2023598``
  - ``interactive -n 80 -t 03:00:00 -A sens2023598`

  - up to 12 hours
  - useful together with the --begin=<time> flag
  - ``salloc -A snic2022-22-50 --begin=2022-02-17T08:00:00`

    - asks for an interactive job that will start earliest tomorrow at 08:00

### Parameters in the job script or the command line?

- Command line parameters override script parameters
- A typical script may be:

```bash
#!/bin/bash
#SBATCH -A sens2023598
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 24:00:00
```
Just a quick test:

```console
sbatch -p devcore -t 00:15:00 jobscript.sh
```

???+ question "Hands-on #1: sbatch/jobinfo"

    - login to Bianca
    - find out which projects you’re a member of using projinfo
    - submit a short (10 min) test job; note the job ID
    - find out if there are any free nodes in the devel partition
    - submit a new job to use the devel partition
    - write in the HackMD when you’re done

### Memory in core or devcore jobs

- ``-n X`
- Bianca: 8GB per core
- Slurm reports the available memory in the prompt at the start of an interactive job

### More flags
- ``-J <jobname>`
- email:

  - ``--mail-type=BEGIN,END,FAIL,TIME_LIMIT_80``
  - ``--mail-user``

    - Don’t use. Set your email correctly in SUPR instead.

- out/err redirection:

  - ``--output=slurm-%j.out`` and ``—-error=slurm-%j.err`

    -  by default, where %j will be replaced by the job ID

  - ``--output=my.output.file``
  - ``--error=my.error.file``


## Monitoring jobs
### Monitoring and modifying jobs
### When a job goes wrong
### Priority

???+ question "Hands-on #2: sbatch/squeue/scancel/scontrol/jobinfo"

    - submit a new job; note the job ID
    - check all your running jobs
    - what is the priority or your recently-submitted job?
    - submit a new job to run for 24h; note the job ID
    - modify the name of the job to “wrongjob” and the maximum runtime to 7days, for example
    - cancel your job with name “wrongjob”

## Determining job efficiency
### Job efficiency

???+ question "Hands-on #3: jobstats"

    - Generate jobstats plots for your jobs
        - Firstly, find some job IDs from this month
        -  finishedjobinfo -m username
        - Write down the IDs from some interesting jobs.
        - Generate the images:
        ```console
        $ jobstats -p ID1 ID2 ID3
        ```
    - Look at the images
       
    ```console
    $ eog *png &
    ```

    - Which of the plots
        - Show good CPU or memory usage?
        - Indicate that the job requires a fat node?

## Different flavours of Slurm: Job script examples and workflows

### Simple workflow

### Job dependencies

### I/O intensive jobs: $SNIC_TMP

### OpenMP or multi-threaded job

### GPU nodes on Bianca

### Running on several nodes: MPI jobs

### Job arrays

### Snakemake and Nextflow 

???+ question "Hands-on #4: make it your own"

    - use 2 or 3 of the sample job scripts as a starting point for your own job script
    - tweak them so that you run something closer to your research; or just feel free to experiment
    - paste at least one of the examples in the HackMD
    - great if you could add a comment what the job script is about

## Feedback on Slurm

## Where to go from here?

Code documentation
SNIC training newsletter - software-specific training events included
https://coderefinery.org/workshops/upcoming/
https://nbis.se/training/events.html (bio)
email support@uppmax.uu.se or https://supr.naiss.se/support/
