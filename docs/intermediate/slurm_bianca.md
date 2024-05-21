# Slurm at Bianca

## More Slurm and other advanced UPPMAX techniques

- A closer look at Slurm
- Using the GPUs
- Job efficiency with the ``jobstats`` tool
- Advanced job submission

## The Slurm Workload Manager
## More on sbatch
## More on time limits
## Job walltime
## More on partitions
## Quick testing
## Debugging or complicated workflows
## Parameters in the job script or the command line?

???+ question "Hands-on #1: sbatch/jobinfo"

    - login to Bianca
    - find out which projects you’re a member of using projinfo
    - submit a short (10 min) test job; note the job ID
    - find out if there are any free nodes in the devel partition
    - submit a new job to use the devel partition
    - write in the HackMD when you’re done

## Memory in core or devcore jobs
## More flags
## Monitoring jobs
## Monitoring and modifying jobs
## When a job goes wrong
## Priority

???+ question "Hands-on #2: sbatch/squeue/scancel/scontrol/jobinfo"

    - submit a new job; note the job ID
    - check all your running jobs
    - what is the priority or your recently-submitted job?
    - submit a new job to run for 24h; note the job ID
    - modify the name of the job to “wrongjob” and the maximum runtime to 7days, for example
    - cancel your job with name “wrongjob”








