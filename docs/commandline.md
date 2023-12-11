# The command-line on Bianca

!!! info "Objectives"
    - Being able to navigate in/out folders
    - Being able to view/create/move/delete files
    - Create an executable bash script

???- info "Notes for teachers"

    Teaching goals:

    - The learners demonstrate they can use a text editor
    - The learners demonstrate they can create, move and delete files
    - The learners demonstrate they can create and delete folders
    - The learners demonstrate they can create an executable script

    Schedule (15 minutes):

    - x minutes: [...]

## Overview

![Using the command line on a computer cluster](./img/610803_a_woman_using_the_command_line_on_a_computer_cluster_256_x_256.png)

Bianca is a cluster with the Linux operating system.
We must use a Linux terminal to work with Bianca,
therefore we must learn some Linux commands.

We will learn to:

- read the manual
- navigate through the file system
- work with directories
- work with files
- create an executable script

!!!- tip "Video: using the command-line on the UPPMAX Bianca cluster"

    There exists a video called 'Using the command-line on the UPPMAX Bianca cluster': 
    [YouTube](https://youtu.be/kjqLAx2bgJI), [download (.ogv)](https://richelbilderbeek.nl/uppmax_bianca_command_line.ogv)

### Read the manual

Use `man` to see the help pages about a command:

```
man man
man cd
man ls
```

These command give the help pages about the programs `man`, `cd` and `ls` respectively.

Press `q` (short for 'quit') to exit `man`

### Navigate through the file system

Like any operating system, Linux has directories (also called 'folders').

Use `cd` to change directory:

Where to                           |Example command
-----------------------------------|---------------------
The root folder                    |`cd /`
The project folder                 |`cd /proj/sens2023598`
Your home folder, using full path  |`cd /home/richel`
Your home folder, using tilde      |`cd ~`
The wharf                          |`cd /proj/sens2023598/nobackup/wharf`
Up one folder                      |`cd ..`
Into a folder, using relative path |`cd myfolder`
The previous location              |`cd -`

!!! tip "See the content of a folder"

    Use `ls` to see the content of a folder

!!! tip "See the current location"

    Use `pwd` to see your current location

!!! info "The Silence Is Golden Rule"
    When your command 'just works' there is no output
    (try, for example `cd ~`). 
    This is due to [The Silence Is Golden Rule](https://www.linfo.org/rule_of_silence.html)

### Work with directories

Linux can create, move and delete folders.

Do what                            |Example command
-----------------------------------|---------------------
Create a folder                    |`mkdir myfolder`
Move a folder                      |`mv from_folder to_folder`
Delete an empty folder             |`rmdir myfolder`
Delete a folder                    |`rm -r myfolder`

!!! tip "See the content of a folder"

    Use `ls` to see the content of a folder

!!! tip "See the current location"

    Use `pwd` to see your current location

!!!- tip "See the current real location (advanced)"

    For sysadmins: use `pwd -P` to see your real current location on the hardware

### Work with files

Linux can create, view, rename, move and delete files.
Additionally, there are some text editors that
allow one to edit files.

Do what                            |Example command
-----------------------------------|---------------------
Create an empty file               |`touch myfile.txt`
View a file using `cat`            |`cat myfile.txt`
Edit a file using `nano`           |`nano myfile.txt`
Delete a file                      |`rm myfile.txt`
Copy a file                        |`cp myfile.txt mycopy.txt`
Rename a file                      |`mv myfile.txt mycopy.txt`
Move a file to one folder up       |`mv myfile.txt ../`
Move a file to the home folder     |`mv myfile.txt ~`

 * Note: `nano` is one of many text editors. 
   It is the one recommended to beginners, 
   as its interface is closest to what one expects

### Create an executable script

Creating an executable script has two steps:

- 1. Create a script
- 2. Allow the script to execute

As an example, we create a script, called `do_it.sh`:

```
nano do_it.sh
```

!!! info "Why use a `.sh` file extension?"
    Using `.sh` as a file extension a social convention for how a Bash script is called,
    as (1) `sh` is short for 'shell', (2) Bash is short for 'Bourne Again Shell'.
    A 'shell' in this context is a program that allows working with an operating system. 

As an example, copy-paste this content into the script:

```
#!/bin/bash
echo "Hello!"
ls | rev
```

!!! info "What does this program do?"
     * The first line is called the [shebang](https://en.wikipedia.org/wiki/Shebang_(Unix)),
       and indicates this is a Bash script
     * The second line displays the text between double quotes
     * The third line displays the files in the folder reversed. 
       The `|` is called the [pipeline](https://en.wikipedia.org/wiki/Pipeline_(Unix)) operator

Save and close `nano`.

 * Use `CTRL-X` to start to exit, then press `y` to start saving the file, then
   press enter to use the current filename

Use [chmod](https://en.wikipedia.org/wiki/Chmod) to make the file executable:

```
chmod +x do_it.sh
```

 * `+x` can be read as: 'add the right to execute'

!!! info "Create read-only files"
    If you want to protect your data from being modified accidentally,
    `chmod` can create read-only files,
    by removing the writing rights using `chmod -w`.

## Exercises

???- tip "Video with solutions"

    There is a video that shows the solution of all these exercises: 
    [YouTube](https://youtu.be/7_LPeQbcmAo), [download (.ogv)](https://richelbilderbeek.nl/bianca_cli.ogv)

???- question "1. View the help of the command `cd`"

    Use `man` to view the help of any command, in this case `cd`:

    ```
    man cd
    ```
    This will fail, because Bianca has (close to) no internet access.

???- question "2a. Navigate to the project folder, e.g. `/proj/sens2023598`"

    ```
    cd /proj/sens2023598
    ```

    Don't forget the `/` at the start.

???- question "2b. Navigate to your home folder"

    The syntax to move to your home folder is:

    ```
    cd /home/[username]
    ```

    where `[username]` is your UPPMAX username, for example:

    ```
    cd /home/richel
    ```

    The squiggle/tilde (`~`) is a shorter notation, 
    that does exactly the same:

    ```
    cd ~
    ```

???- question "2c. Navigate to the wharf, e.g. `/proj/sens2023598/nobackup/wharf`"

    ```
    cd /proj/sens2023598/nobackup/wharf
    ```

???- question "3a. Create a folder `/proj/sens2023598/workshop/[your_login_name]`, for example, `/proj/sens2023598/workshop/richel`"

    ```
    mkdir /proj/sens2023598/workshop/richel
    ```

    Or navigate there first:

    ```
    cd /proj/sens2023598/workshop/
    mkdir richel
    ```


???- question "4a. Create a file, e.g. `richel.txt`"

    ```
    touch richel.txt
    ```


???- question " 4b. Copy the file (e.g. to `richel_again.txt`). "

    ```
    cp richel.txt richel_again.txt
    ```

???- question "4c. Move the copied file (e.g. move it one folder up to `../richel_again.txt`)"

    ```
    mv richel_again.txt ../
    ```

???- question "4d. Delete the copied file"


    ```
    rm ../richel_again.txt
    ```

    or:

    ```
    cd ..
    rm richel_again.txt
    ```


???- question "5. Create an executable script called `/proj/sens2023598/workshop/[your_login_name]/do_it.sh`, which, upon running, displays a welcome message in text (e.g. `Hello!`) and does something (e.g. show the files in reverse order)"

    Edit the script:

    ```
    nano do_it.sh
    ```

    Change the text to:

    ```
    #!/bin/bash
    echo "Hello!"
    ls | rev
    ```

    Run the script:

    ```
    ./do_it.sh
    ```

## Extra material

### Other useful commands

These are some commands that we enjoy,
but are not part of the learning objectives.

Command name|Purpose
------------|---------------------------------------------
`scp`       |Copy file between Bianca and your local computer
`cat`       |Show the content of a file
`less`      |Navigate through the content of a file
`head`      |Show the first lines of a file
`tail`      |Show the last lines of a file
`less`      |Show the content of a file
`wc`        |Count words, lines and/or characters
`|`         |[Pipe](https://en.wikipedia.org/wiki/Pipeline_(Unix)) the output of one command to serve as input for the next
`>`         |Write to file (removes existing content if any)
`>>`        |Append to file

With `ls /usr/bin | wc --lines` one can see that there are more than 1700
commands on Bianca.

### The terminal and the GUI are friends

If you are using the Bianca remote desktop environment,
you can see that its file browser and terminal are friends.

On a clean terminal, try typing `cd` 
and then drag a folder from the GUI to the terminal.

It types the absolute path for you!

### Commonly used symbolic links

These are some commonly used symbolic links, 
that will simplify navigation:

```
cd Desktop
ln -s /proj/sens2023598/ proj
ln -s /proj/sens2023598/nobackup nobackup
ln -s /proj/sens2023598/nobackup/wharf/richel/richel-sens2023598 wharf`
```

 * Replace `sens2023598` by your project
 * Replace `richel` by your username
