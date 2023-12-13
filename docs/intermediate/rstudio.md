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

## Procedure to start RStudio

### 1. Get within SUNET

???- tip "Forgot how to get within SUNET?"

    See [the basic Bianca course page 'Login to Bianca'](../login_bianca.md).

### 2. Start the Bianca remote desktop environment

???- tip "Forgot how to start Bianca's remote desktop environment?"

    See [the basic Bianca course page 'Login to Bianca'](../login_bianca.md).

### 3. Start an interactive session

Within the Bianca remote desktop environment, start a terminal.
Within that terminal, start an interactive session with 2 cores.

!!!- info "Why two cores?"

    RStudio is a resource-heavy program.
    Due to this, we recommend using at least two cores 
    for a more pleasant user experience.

???- tip "Forgot how to start an interactive node?"

    See [the basic Bianca course page 'Starting an interactive node'](../start_interactive_node.md).

    Spoiler: use:

    ```
    interactive -A sens2023598 -n 2 -t 8:00:00
    ```

### 4. Load the modules needed

RStudio needs R and its R packages.
These should be loaded via the module system.

In the terminal of the interactive session, do:

```
module load R_packages/4.3.1 RStudio/2023.06.2-561
```

### 5. Start RStudio

With the modules loaded, start RStudio from the terminal (on the
interactive node):

```
rstudio
```

RStudio can be slow to startup, as R has thousands (!) of packages.
Additionally, at startup and if enabled, your saved RStudio workspace
(with potentially a lot of data!) is read.

## Questions

