# Create an executable script

![Using the command line on a computer cluster](./img/610803_a_woman_using_the_command_line_on_a_computer_cluster_256_x_256.png)

!!! info "Objectives"

    - Create an executable bash script
    - Optional: Being able to navigate in/out folders
    - Optional: Being able to view/create/move/delete files

???- info "Notes for teachers"

    Teaching goals:

    - The learners demonstrate they can use a text editor
    - The learners demonstrate they can create an executable script
    - Optional: The learners demonstrate they can create, move and delete files
    - Optional: The learners demonstrate they can create and delete folders

    Schedule:

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

## Why?

You need the command-line to start calculations.

## Overview

Bianca is a cluster with the Linux operating system.
We must use a Linux terminal to work with Bianca,
therefore we must learn some Linux commands.

We will learn to:

- create an executable script

### Create an executable script

Creating an executable script has two steps:

- Create a script
- Allow the script to execute

As an example, we create a script, called `do_it.sh`:

```bash
nano do_it.sh
```

!!! info "Why use a `.sh` file extension?"
    Using `.sh` as a file extension a social convention
    for how a Bash script is called,
    as (1) `sh` is short for 'shell',
    (2) Bash is short for 'Bourne Again Shell'.
    A 'shell' in this context is a program
    that allows working with an operating system.

As an example, copy-paste this content into the script:

```bash
#!/bin/bash
echo "Hello!"
```

!!! info "What does this program do?"
    
    The first line is called the [shebang](https://en.wikipedia.org/wiki/Shebang_(Unix)),
    and indicates this is a Bash script.
    
    The second line displays the text between double quotes.

Save and close `nano`.

- Use `CTRL-X` to start to exit, then press `y` to start saving the file, then
  press enter to use the current filename

Use [chmod](https://en.wikipedia.org/wiki/Chmod) to make the file executable:

```bash
chmod +x do_it.sh
```

- `+x` can be read as: 'add the right to execute'

!!! info "Create read-only files"
    If you want to protect your data from being modified accidentally,
    `chmod` can create read-only files,
    by removing the writing rights using `chmod -w`.


### Exercise 1: create an executable script

- Create a file called `do_it.sh` using `nano`

???- question "Answer"

    Use `nano` to create it:

    ```bash
    nano do_it.sh
    ```

    Then do `CTRL + O` to save, `CTRL + X` to exit

- Edit the file `do_it.sh` to have the content below:

```bash
#!/bin/bash
echo "Hello!"
```

???- question "Answer"

    Use `nano` to edit it:

    ```bash
    nano do_it.sh
    ```

    Then do `CTRL + O` to save, `CTRL + X` to exit

- Write an executable script that displays a welcome message in text (e.g. `Hello!`)

???- question "Answer"

    Edit the script:

    ```bash
    nano do_it.sh
    ```

    Change the text to:

    ```bash
    #!/bin/bash
    echo "Hello!"
    ls | rev
    ```

    Make the script executable:

    ```bash
    chmod +x ./do_it.sh
    ```

- Run the script

???- question "Answer"

    Run the script:

    ```bash
    ./do_it.sh
    ```

