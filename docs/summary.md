# Summary

## Next steps

- Q/A session
- Pick topics from Extra material or ask more general questions, see below summary
- Exercises ready in the end of the week.

## Today's topics

- **Overview of UPPMAX systems**

    - Several clusters, like 

        - Bianca
        - Rackham

    - Login and compute nodes

- **Logging in**
 
    - ThinLinc, a must for graphics
    - ``ssh``, may be faster if you just use command-line
    - VPN
    - 2FA


- **Command-line intro**

    - navigation 
    - aliases
    - chmod


- **Module system and workflows**


    - modules adds paths to tools


        - module load ...
        - bioinfo-tools may be required to load you tool


    - workflows


- **Transferring files**


    - The ``wharf`` works like a dock at the harbour.
    - There are several ways to use the ``wharf`` to transfer files


      - copy
      - transit server
      - rsync, scp/sftp


- **NAISS-SENS and sensitive data**


    - Sensitive personal data
    - Pseudonymisation and anonymisation
    - Apply for project


- **Introduction to compute nodes**

    - Submitting jobs


        - Slurm, sbatch, the job queue
        - Jobs
        - Slurm parameters
        - Interactive jobs 

            - allows you to work interactively with command line and your tools but on a compute node.

    - How compute nodes are moved between project clusters

    - Other Slurm tools

- **Software and package installation**


    - Install software yourself
    - Packages and libraries to scripting programs

        - Conda
        - Python packages with pip
        - R packages
        - Julia packages

    - "Containers"

        - Singularity
        - Docker

    - Build from source
    - Summary about the Bianca Hardware
