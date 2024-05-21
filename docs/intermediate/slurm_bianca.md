# Slurm at Bianca

## More Slurm and other advanced UPPMAX techniques

- A closer look at Slurm
- Using the GPUs
- Job efficiency with the ``jobstats`` tool
- Advanced job submission

## The Slurm Workload Manager
Free, popular, lightweight
Open source:
https://slurm.schedmd.com
available at al SNIC centra
UPPMAX Slurm userguide:
https://www.uppmax.uu.se/support/user-guides/slurm-user-guide/

### More on sbatch
Recap:

sbatch | -A sens2023598  |   -t 10:00 | -p core | -n 10 | my_job.sh
-|-|-|-|-|-
slurm batch| project name | max runtime | partition ("job type") | #cores | job script




### More on time limits
### Job walltime
### More on partitions
### Quick testing
### Debugging or complicated workflows
### Parameters in the job script or the command line?

???+ question "Hands-on #1: sbatch/jobinfo"

    - login to Bianca
    - find out which projects you’re a member of using projinfo
    - submit a short (10 min) test job; note the job ID
    - find out if there are any free nodes in the devel partition
    - submit a new job to use the devel partition
    - write in the HackMD when you’re done

### Memory in core or devcore jobs
### More flags

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
