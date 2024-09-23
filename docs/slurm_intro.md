# Introduction to compute nodes

!!! info "Learning objectives"
    - This is a short introduction in how to reach the calculation/compute/worker nodes
    - The learners should be able to:
        - Run simple jobs in the batch system
        - Run interactively on compute nodes
        - When interactive and when batch
    -     Check your jobs


???- info "Notes for teachers"

    Teaching goals:

    - The learners demonstrate to have run in interactive 
    - The learners demonstrate to have run batch job
    - The learners demonstrate to have understood when to use batch or interactive 
    - The learners demonstrate to have understood how to plan for jobs

    Schedule (45 minutes):

    - 5 minutes: lecturing
    - 15 minutes type-alongs x 2
    - 20 minutes: exercise + quiz
    - 5 minutes: discuss answers

## Exercises

???+ question "You are developing code on Bianca."

    - You write the code line-by-line and schedule a test run after each addition. 
    - However, after each new line, it takes a couple of minutes before you know your code worked yes/no. 
    - How could you develop your code quicker?"

    ??? tip "Answer"
    
        - This is the typical use-case to use an interactive node.
        - One could also consider to develop code on a local computer instead (which uses nonsensitive/simulated/fake testing data) and upload the final code instead.
    
??? question "Start an interactive session"

    The goal of this exercise is to make sure you know how to start an 
    interactive session. 

???- question "Why not always use an interactive session?"

     - Because it is an inefficient use of your core hours.

     - An interactive session means that you use a calculation node with low efficiency: only irregularly you will use such a node to its full
    capacity. 
     - However, the number of core hours are registered as if the node is used at full capacity, as it is *reserved* to be used at that capacity.

???+ question "Which approach is best in the following use cases? Batch jobs or interactive sessions?"

    1.  Long jobs
    1.  Short jobs with interactive "run-time"/interactive user input
    1.  Short jobs without interactive "run-time"/interactive user input
    1.  Test/debugging/developing code
    1.  Playing with and plotting large data

    ??? tip "Answer"

        1.  batch
        1.  interactice
        1.  batch
        1.  interactive
        1.  interactive



???+ question "Submit a Slurm job"

    - Make a batch job to run the [demo](https://uppmax.github.io/bianca_workshop/modules1/#bigger-exercises) "Hands on: Processing a BAM file to a VCF using GATK, and annotating the variants with snpEff". Ask for 2 cores for 1h.
        - You can copy the my_bio_workflow.sh file in ``/proj/sens2023598/workshop/slurm`` to your home folder and make the necessary changes.
    
    ??? tip "Answer"
        - edit a file using you preferred editor, named `my_bio_worksflow.sh`, for example, with the content
        - alternatively copy the ``/proj/sens2023598/workshop/slurm/my_bio_workflow.sh`` file and modify it
          ``cd ~`` 
          ``cp /proj/sens2023598/workshop/slurm/my_bio_workflow.sh .``
          - edit ``my_bio_workflow.sh`` and add the SBATCH commands
        
        ```bash
        #!/bin/bash
        #SBATCH -A sens2023598
        #SBATCH -J workflow
        #SBATCH -t 01:00:00
        #SBATCH -p core
        #SBATCH -n 2


        cd ~
        mkdir -p myworkflow
        cd myworkflow

        module load bioinfo-tools

        # load samtools
        module load samtools/1.17

        # copy and example BAM file
        cp -a /proj/sens2023598/workshop/data/ERR1252289.subset.bam .

        # index the BAM file
        samtools index ERR1252289.subset.bam

        # load the GATK module
        module load GATK/4.3.0.0

        # make symbolic links to the hg38 genomes
        ln -s /sw/data/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.* .

        # create a VCF containing inferred variants
        gatk HaplotypeCaller --reference genome.fa --input ERR1252289.subset.bam --intervals chr1:100300000-100800000 --output ERR1252289.subset.vcf

        # use snpEFF to annotate variants
        module load snpEff/5.1
        java -jar $SNPEFF_ROOT/snpEff.jar eff hg38 ERR1252289.subset.vcf > ERR1252289.subset.snpEff.vcf

        # compress the annotated VCF and index it
        bgzip ERR1252289.subset.snpEff.vcf
        tabix -p vcf ERR1252289.subset.snpEff.vcf.gz
        ```

        - make the job script executable
        ```bash
        $ chmod a+x my_bio_workflow.sh
        ```
        
        - submit the job
        ```bash
        $ sbatch my_bio_workflow.sh
        ```
## Links

- [Slurm documentation](https://slurm.schedmd.com/){:target="_blank"}
- [New Slurm user guide (needs updates)](https://uppmax.github.io/UPPMAX-documentation/cluster_guides/slurm/){:target="_blank"}
- [Discovering job resource usage with `jobstats`](http://docs.uppmax.uu.se/software/jobstats/){:target="_blank"} 
- [Plotting your core hour usage](http://docs.uppmax.uu.se/software/projplot/){:target="_blank"} 

!!! example "Discussion"

    - Any further thoughts?

!!! abstract "Keypoints"
    - Slurm is a job scheduler to handle the compute nodes
        - add flags to describe your job.
    - You are always in the login node unless you:
        - start an interactive session to do development or hands-on work
        - start a batch job to run jobs not needing any manual input
    - There is a job wall time limit of ten days (240 hours).
 
