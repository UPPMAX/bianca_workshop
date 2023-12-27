# IDE:s

![](./img/rstudio_in_action_480_x_270.png)

> RStudio is one of the IDEs that can be used on Bianca.

!!! info "Objectives" 

    - Observe different IDEs running on Bianca
    - Start your favorite IDE on Bianca

???- info "Notes for teachers"

    Teaching goals:

    - Show the three IDEs in action
    - The learners demonstrate to have started at least on IDE on Bianca

    Schedule (45 minutes):

    - 5 mins: Let the learners start an interactive node: this can take dozens of minutes!
    - 10 mins: discuss this page and its sub-pages
    - 15 mins: do the exercises
    - 15 mins: discuss the exercises

## Exercises

Read the UPPMAX documentation on IDEs on Bianca
[here](https://uppmax.github.io/UPPMAX-documentation/cluster_guides/ides_on_bianca/),
then do these exercises.

???- question "Exercise: Start your favorite IDE"

    The goal of this exercise is to make sure you can start
    at least 1 IDE.

???- question "Why use an IDE?"

    It makes new Bianca users feel comfortable,
    as an IDE is a recognizable environment.
    Also, the terminal can be daunting to some.

    Additionally, an IDE allows one to do runtime debugging,
    i.e. running through code line-by-line and/or up/down
    the so-called call stack.

???- question "Why not always use an IDE?"

    Using an IDE on Bianca is cumbersome and
    there are superior ways to develop code on Bianca,
    as -for example- taught in the 
    [UPPMAX 'Programming Formalisms' course](https://github.com/UPPMAX/programming_formalisms).

???- question "Why not always use an interactive session?"

    Because it is an inefficient use of your core hours.

    An interactive session means that you use a calculation node with low
    efficiency: only irregularly you will use such a node to its full
    capacity. 
    However, the number of core hours are registered as if the node is used
    at full capacity, as it is *reserved* to be used at that capacity.

???- question "How to find out if you are on a login or interactive node"

    In the terminal, type `hostname`

    - the login node has `[project]-bianca`, where `[project]` is the name of the project, e.g. `sens2023598`
    - the interactive node has `b[number]` in it, where `[number]` is the compute node number

