# IDEs

![RStudio running on Bianca](./img/rstudio_in_action_480_x_270.png)

> RStudio is one of the IDEs that can be used on Bianca.

!!! info "Objectives"

    - Find information in the UPPMAX documentation about IDEs on Bianca
    - Give a reasonable definition of what an IDE is
    - Remember that RStudio, Jupyter, VSCodium are IDEs that can be run on Bianca
    - Can run the voted-for IDE on Bianca
    - (optional) Can give a reason why to use an IDE
    - (optional) Can give a reason why not to use an IDE on Bianca
    - (optional) Can give a reason why not to run an IDE on a login node
    - (optional) Can give a reason when to use an interactive session
    - (optional) Can find out if an interactive session is active

???- info "Notes for teachers"

    Teaching goals:

    - The learners have explored the UPPMAX documentation
    - The learners have seen that there are different IDEs on Bianca
    - The learners have start at least one IDE on Bianca
    - The learners understand why to use an IDE
    - The learners understand why not to use an IDE on Bianca
    - The learners understand when to run an IDE on a login node
    - The learners understand when to use an interactive session
    - The learners understand how to find out an interactive session is active

    Lesson plan:

    ```mermaid
    gantt
      title IDEs
      dateFormat X
      axisFormat %s
      Introduction: intro, 0, 5s
      Vote on whcih IDE: vote, after intro, 5s
      Exercise with winning IDE: crit, exercise, after vote, 20s
      Feedback: feedback, after exercise, 10s
      Monologue other 2 IDEs: monologue, after feedback, 5s
      Break: milestone, after monologue
    ```

## Why?

You want to develop/modify code on Bianca in a program ...

- ... that you already use on your regular computer
- ... that is not the terminal
- ... that helps you do so by providing code completion,
      code hints, run-time debugging, etc.

Hence, you want to use an IDE.

IDE (pronounce `aj-dee-ee`) is short for 'Integrated Development Environment',
or 'a program in which you do programming'.
The goal of an IDE is to help develop code, with features
such as code completion, code hints and interactive debugging.

Using an IDE on Bianca is cumbersome and
there are superior ways to develop code on Bianca,
as -for example- taught in the
[UPPMAX 'Programming Formalisms' course](https://github.com/UPPMAX/programming_formalisms).

There are three [IDEs on Bianca](http://docs.uppmax.uu.se/software/ides_on_bianca/)
we can show: Jupyter, RStudio or VSCodium.
We will practice and discuss one, and briefly talk about the others.


## Procedure

We pick the winning IDE democratically:

- In the shared document, add a character before each IDE you'd be interested
  in, between the `[ ]` of each option. This will be messy.
  You can vote for 0, 1, 2 or 3 IDEs.

It may look like below, where there is a preference for Jupyter:

```text
- [,a/b*efgh!] Jupyter
- [a/d,h] RStudio
- [cd!h]VSCodium
```

## Exercises

### Exercise 1: Start the favorite IDE

???- info "Learning objectives"

    - Explore the UPPMAX documentation
    - Start the favorite IDE on Bianca

- Go to the UPPMAX documentation of the [IDEs on Bianca](http://docs.uppmax.uu.se/software/ides_on_bianca/)
- Get the IDE that got the most votes to run on Bianca, by following
  its documentation. If you really want to run another IDE,
  you may do so! When done, do exercise 2.

???- question "If you chose Jupyter"

    You can skip the _venv_ step!

    - If not done so earlier:

    ```console
    $ cd ~   # just 'cd' will work as well
    $ ln -s /proj/sens2023598/ proj
    ```
     
    - Start jupyter from your $HOME folder
    - Try to run the script `Test-01.ipynb` located in /proj/workshop/Jupyter-demo/Test-01.ipynb
    - You can browse in jupyter to a test notebook in /proj/workshop/Jupyter-demo/Test-01.ipynb

### Exercise 2: Understand IDEs on Bianca

???- info "Learning objectives"

    - Explore the UPPMAX documentation
    - Understand why to use an IDE
    - Understand why not to use an IDE on Bianca
    - Understand when to run an IDE on a login node
    - Understand when to use an interactive session
    - Understand how to find out an interactive session is active


- Try to answer the questions below.
  Be generous in accepting you answer.
  If you have no idea, use the UPPMAX documentation of the [IDEs on Bianca](http://docs.uppmax.uu.se/software/ides_on_bianca/).
  When done, run the other two IDEs on Bianca

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

???- question "Is it OK to run IDEs on a login node? Why yes/no?"

    No. IDEs are big programs, use an interactive session instead.

    You could argue if you are the only one on a Bianca project,
    you code use the login node. This only works if the IDE
    works fine on such limited computational resources.

???- question "Why not always use an interactive session?"

    Because it is an inefficient use of your core hours.

    An interactive session means that you use a calculation node with low
    efficiency: only irregularly you will use such a node to its full
    capacity. 
    However, the number of core hours are registered as if the node is used
    at full capacity, as it is *reserved* to be used at that capacity.

???- question "How to find out if you are on a login or interactive node?"

    In the terminal, type `hostname`

    - the login node has `[project]-bianca`, where `[project]` is the name of the project, e.g. `sens2023598`
    - the interactive node has `b[number]` in it, where `[number]` is the compute node number
