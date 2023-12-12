# RStudio

!!! info "Objectives" 

    - Observe different IDEs running on Bianca
    - Start your favorite IDE on Bianca

???- info "Notes for teachers"

    Teaching goals:

    - Show the three IDEs in action
    - The learners demonstrate to have started at least on IDE on Bianca

    Schedule (.. minutes):

    - ...

## Introduction

RStudio is an IDE specialized for the R programming language.

???- tip "What is an IDE?"

    See [the page on IDEs](ides.md).

In this session, we show how to use RStudio on Bianca,
using Bianca's remote desktop environment.

???- tip "Forgot how to login to a remote desktop environment?"

    See [the basic Bianca course page 'Logging in'](../login_bianca.md).

    Spoiler: go to [https://bianca.uppmax.uu.se/](https://bianca.uppmax.uu.se/)

As RStudio is a resource-heavy program,
it must be run on an interactive node.

???- tip "Forgot how to start an interactive node?"

    See [the basic Bianca course page 'Starting an interactive node'](../start_interactive_node.md).

## Starting RStudio

### ?Optional

Start an interactive node:

```
interactive --account sens2017625 --nodes 1 --ntasks 16 --time 8:00:00
interactive --account sens2017625 --nodes 1 --ntasks 16 --time 8:00:00
interactive -A sens2017625 -N 1 -n 16 -t 8:00:00
```

### Load the modules needed

```
R/4.3.1
module load R_packages/4.1.1 RStudio/2022.02.0-443
module load R_packages/4.3.1 RStudio/2023.06.2-561
```

### Start RStudio

```
rstudio &
```



We recommend using at least two cores for RStudio, and to get those resources, you must should start an interactive job.

!!! example "Type-along"
    Use **ThinLinc**

    - Start **interactive session** on compute node (2 cores)
    - If you already have an interactive session going on use that.
        - If you don't find it, do
        
            ``$ squeue``
            
        - find your session, ssh to it, like:
        
            ``$ ssh sens2023598-b9``

    - ``$ interactive -A sens2023598 -p devcore -n 2 -t 60:00`` 


    - Once the interactive job has begun you need to load needed modules, even if you had loaded them before in the login node
    - You can check which node you are on?

        `$ hostname`
    
    - Also try: 

        `$ srun hostname`

        - This will give several output lines resembling the number of cores you allocated.
        - How many in this case??
        
    - If the name before ``.bianca.uppmax.uu.se`` is ending with bXX you are on a compute node!
    - The login node has ``sens2023598-bianca``
    - You can also probably see this information in your prompt, like:
        ``[bjornc@sens2023598-b9 ~]$`` 
  
    - Load an RStudio module and an R_packages module (if not loading R you will have to stick with R/3.6.0) and run "rstudio" from there. 

        `$ ml R_packages/4.2.1`
  
        `$ ml RStudio/2022.07.1-554`


    - **Start rstudio**, keeping terminal active (`&`)

      `$ rstudio &`

    - Slow to start?
    - Depends on:
        - number of packages 
        - if you save a lot of data in your RStudio workspace, to be read during start up.

    - **Quit RStudio**!
    - **Log out** from interactive session with `<Ctrl>-D` or `logout` or `exit`
