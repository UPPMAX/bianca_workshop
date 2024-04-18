# Working with environment modules on Bianca

![Working with a computer cluster module system](./img/627409_working_with_a_computer_cluster_module_system_256_x_256.png)

!!! info "Objectives" 

    - Being able to search/load/unload modules
    - Create an executable Bash script that uses a module (without SLURM)

???- info "Notes for teachers"

    Teaching goals:

    - The learners demonstrate they can find a module
    - The learners demonstrate they can load a module of a specific version
    - The learners demonstrate they can unload a module
    - The learners demonstrate they can load a module in a script

    ```mermaid
    gantt
      title Lesson plan Command line
      dateFormat X
      axisFormat %s
      Prior knowledge: prior, 0, 5s
      Theory: theory, after prior, 5s
      Exercises: crit, exercise, after theory, 25s
      Feedback: feedback, after exercise, 10s
    ```

## Exercises

Use [the UPPMAX documentation on modules](http://docs.uppmax.uu.se/cluster_guides/modules/) 
to do these exercises.

???- tip "Video with solutions"

    There is a video that shows the solution of all these exercises: 
    [YouTube](https://youtu.be/lNlq2Eb-qgc), [Download (.ogv)](https://richelbilderbeek.nl/bianca_modules.ogv)

???- question "1a. Verify that the tool `cowsay` is not available by default"

    ```
    cowsay hello
    ```

    Gives the error message: `cowsay: command not found`.

???- question "1b. Search for the module providing `cowsay`"

    ```
    module spider cowsay
    ```

    You will find the `cowsay/3.03` module.

???- question "1c. Load a specific version of that module"

    ```
    module load cowsay/3.03
    ```

???- question "1d. Verify that the tool `cowsay` now works"

    ```
    cowsay hello
    ```

???- question "1e. Unload that module"


    ```
    module unload cowsay/3.03
    ```

???- question "1f. Verify that the tool `cowsay` is not available anymore"


    ```
    cowsay hello
    ```

    Gives the error message: `cowsay: command not found`.


???- question "2a. Create an executable script called `cow_says_hello.sh`. It should load a specific version of the `cowsay` module, after which it uses `cowsay` to do something"

    ```
    nano cow_says_hello.sh
    ```

    Copy-paste this example text:

    ```
    #!/bin/bash
    module load cowsay/3.03
    cowsay hello
    ```

    Make the script executable:

    ```
    chmod +x cow_says_hello.sh
    ```

    Run:

    ```
    ./cow_says_hello.sh
    ```

???- question "2b. Find out: if the `cowsay` module is not loaded, after running the script, is it loaded yes/no?"
    
    Running the script does not load the module beyond running the script.

    ```
    [richel@sens2023598-bianca ~]$ cowsay hello
    -bash: cowsay: command not found
    [richel@sens2023598-bianca ~]$ ./cow_says_hello.sh 
     ________ 
    < hello >
     -------- 
            \   ^__^
             \  (oo)\_______
                (__)\       )\/\
                    ||----w |
                    ||     ||
    [richel@sens2023598-bianca ~]$ cowsay hello
    -bash: cowsay: command not found
    ```

???- question "3. `module load samtools/1.17` gives the error `These module(s) or extension(s) exist but cannot be loaded as requested: "samtools/1.17`. How to fix this?"

    If you do `module load samtools/1.17` without 
    doing `module load bioinfo-tools` first, you'll get the error:

    ```
    $ module load samtools/1.17
    Lmod has detected the following error:  These module(s) or
    extension(s) exist but cannot be loaded as requested: "samtools/1.17"
       Try: "module spider samtools/1.17" to see how to load the module(s).
    ```

    The solution is to do `module load bioinfo-tools` first.

???- info "Want more complex/realistic exercises?" 

    The goal of this lesson is to work with the module system
    in a minimal/fast way. 
    These exercises do not achieve anything useful.
    See 'Bigger exercises' for more complex/realistic exercises


## Bigger exercises

!!! warning 
    - To access bioinformatics tools, load the **bioinfo-tools** module first.

???- question "Hands on: Processing a BAM file to a VCF using GATK, and annotating the variants with snpEff"

    This workflow uses a pre-made BAM file that contains a subset of reads from a sample from European Nucleotide Archive project [PRJEB6463](https://www.ebi.ac.uk/ena/browser/view/PRJEB6463) aligned to human genome build hg38. These reads are from the region `chr1:100300000-100800000`.

    1. Copy example BAM file to your working directory.
    ```
    $ cp -a /proj/sens2023598/workshop/data/ERR1252289.subset.bam .
    ```

    2. Take a quick look at the BAM file. First see if `samtools` is available.
    ```
    $ which samtools
    ```

    3. If `samtools` is not found, load `bioinfo-tools` then `samtools/1.17`
    ```
    $ ml bioinfo-tools samtools/1.17
    ```

    4. Now create an index for the BAM file, and examine the first 10 reads aligned within the BAM file.
    ```
    $ samtools index ERR1252289.subset.bam
    $ samtools view ERR1252289.subset.bam | head
    ```

    5. Looks good. Now load the `GATK/4.3.0.0` module.
    ```
    $ module load GATK/4.3.0.0
    ```

    6. Make symbolic links to hg38 genome resources already available on UPPMAX. This provides local symbolic links for the hg38 resources `genome.fa`, `genome.fa.fai` and `genome.dict`.
    ```
    $ ln -s /sw/data/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.* .
    ```

    7. Create a VCF containing inferred variants. Speed it up by confining the analysis to this region of chr1.
    ```
    $ gatk HaplotypeCaller --reference genome.fa --input ERR1252289.subset.bam --intervals chr1:100300000-100800000 --output ERR1252289.subset.vcf
    ```
    This produces as its output the files `ERR1252289.subset.vcf` and `ERR1252289.subset.vcf.idx`.

    8. Now use `snpEff/5.1` to annotate the variants. Loading `snpEff/5.1` results in a change of java prerequisite. Also take a quick look at the help for the module for help with running this tool at UPPMAX.
    ```
    $ ml snpEff/5.1

    The following have been reloaded with a version change:
      1) java/sun_jdk1.8.0_151 => java/OpenJDK_12+32

    $ ml help snpEff/5.1

    ------------------- Module Specific Help for "snpEff/5.1" --------------------
        snpEff - use snpEff 5.1
        Version 5.1

        Usage: java -jar $SNPEFF_ROOT/snpEff.jar ...

        Usage: java -jar $SNPEFF_ROOT/SnpSift.jar ...
        along with the desired command and possible java options for memory, etc

        Note that databases must be added by an admin -- request via support@uppmax.uu.se
        See http://snpeff.sourceforge.net/protocol.html for general help

    Every database that is provided by snpEff/5.1 as of this installation is installed.  This complete list
    can be generated with

        java -jar $SNPEFF_ROOT/snpEff.jar databases

    Three additional databases have been installed.

        Database name                  Description                                      Notes
        -------------                  -----------                                      -----
        c_elegans.PRJNA13758.WS283     Caenorhabditis elegans genome version WS283      MtDNA uses Invertebrate_Mitochondrial codon table
        canFam4.0                      Canis familiaris genome version 4.0
        fAlb15.e73                     Ficedula albicollis ENSEMBLE 73 release

    The complete list of locally installed databases is available at $SNPEFF_ROOT/data/databases_list.installed

    To add your own snpEff database, see the guide at http://pcingola.github.io/SnpEff/se_buildingdb/#option-1-building-a-database-from-gtf-files
    ```

    9. Annotate the variants.
    ```
    $ java -jar $SNPEFF_ROOT/snpEff.jar eff hg38 ERR1252289.subset.vcf > ERR1252289.subset.snpEff.vcf
    ```

    10. Take a quick look!
    ```
    $ less ERR1252289.subset.snpEff.vcf
    ```

    11. Compress the annotated VCF and index it, using `bgzip` and `tabix` provided by the `samtools/1.17` module, already loaded.
    ```
    $ bgzip ERR1252289.subset.snpEff.vcf
    $ tabix -p vcf ERR1252289.subset.snpEff.vcf.gz
    ```


???- question "Hands on: Running R within RStudio, use ggplot2 from R_packages/4.1.1"

    1. Load the `R_packages/4.1.1` module and the latest `RStudio` module, and start RStudio with `rstudio &`.
    ![ThinLinc load R_packages RStudio](./img/modules-1-ml-rstudio.png)

    2. Load the `ggplot2` R library, provided by `R_packages/4.1.1`, and produce an example plot.
    ![ThinLinc ggplot2](./img/modules-2-ggplot2.png)

    2. Save the plot using `ggsave`.
    ![ThinLinc ggsave](./img/modules-3-ggsave.png)


???- question "Hands on: Loading the conda/latest module"

    1. Load the `conda/latest` module.
    ```
    $ ml conda/latest
    The variable CONDA_ENVS_PATH contains the location of your environments. Set it to your project's environments folder if you have one.
    Otherwise, the default is ~/.conda/envs. Remember to export the variable with export CONDA_ENVS_PATH=/proj/...

    You may run "source conda_init.sh" to initialise your shell to be able
    to run "conda activate" and "conda deactivate" etc.
    Just remember that this command adds stuff to your shell outside the scope of the module system.

    REMEMBER TO USE 'conda clean -a' once in a while
    ```

    We want to set the `CONDA_ENVS_PATH` variable to a directory within our project, rather than use the default which is our home directory.
    If you do not set this variable, your home directory will easily exceed its quotas when creating even a single Conda environment.
    This will be covered in more detail in the afternoon.
