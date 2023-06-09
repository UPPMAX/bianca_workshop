# The command-line on Bianca

!!! info "Objectives"
    - We'll use the commands and investigate the Bianca environment
    - Tips and tricks for Bianca users

!!! warning

    - We assume that you have already covered the Command-line material and tested on Rackham
        - [LINUX](https://uppmax.github.io/uppmax_intro/linux.html)
        - [Basic toolkit](https://uppmax.github.io/uppmax_intro/linux_basics.html)

## Command-line intro

       
### Navigation and file management

1. `pwd`  &emsp; present directory
1. `ls`  &emsp;list content
1. `touch` &emsp;create an empty file
1. `cd`  &emsp;change directory
1. `mkdir`  &emsp;make directory
1. `cp`  &emsp;copy
1. `scp`  &emsp;securely remotely copy
1. `mv`  &emsp;move
1. `rm`  &emsp;remove
1. `rmdir`  &emsp;remove empty directory

### Read files and change file properties

10. `cat`  &emsp;print content on screen
11. `head`  &emsp;print first part
12. `tail`  &emsp;print last part
13. `less`  &emsp;browse content
14. `tar`  &emsp;compress or extract file
15. `chmod`  &emsp;change file permissions
16. `man`  &emsp;info about a command
17. `echo`  &emsp;print text on screen

## Type along

### Navigating Bianca

- Check the path to your $HOME folder

```
$ cd ~
$ pwd
$ pwd -P
```

??? answer
    ```
    /home/$USER
    /castor/project/home/bjornc
    ```
    The first `cd ~` gives no output, due to the 'Silence is golden' rule,
    which entails to display nothing if all went OK.

    It appears that your home folder, is some else's subfolder :-)

- Check the path to your projects

```
$ cd /proj
$ ls
$ pwd
$ pwd -P
```

??? answer
    ```
    /proj
    /proj
    ```
    Your `/proj` folder, is the machine's `/proj` folder.

```
$ cd /sensXXX
$ pwd
$ pwd -P
```
??? answer
    ```
    /proj/sensXXX
    /castor/project/proj
    ```
    Your `/proj/sensXXX` folder, however, is a subfolder of `/castor`

### Create empty file

Create an empty file, show it exists, and delete it:

```
$ touch my_file.txt
$ ls
$ rm my_file.txt
```

??? answer
    There is only output after `ls`:
    ```
    my_file.txt
    ```
    The exact output depends on all files present in your working directory

### Create simple script

Create a simple script, run it, and delete it:

```
$ echo '#!/bin/bash' > my_script.sh
$ echo 'echo "Hello world"' >> my_script.sh
$ chmod +x my_script.sh
$ ./my_script.sh
$ rm my_script.sh
```

Note that `>` means '(over)write the text to file',
and `>>` means 'append the text to file'. Also, single quotes (`'`)
for a cleaner syntax (i.e. avoiding `\"`).

??? answer
    There is only output after running the script:
    ```
    Hello world
    ```

## The terminal and the GUI are friends

On a clean terminal, try typing `cd` 
and then dragging a folder from the GUI to the terminal.

It types the absolute path for you!

## Use `chmod` to protect important files

Most projects rely on data that should in principle never be changed. 
Run `chmod -R -w data/` to remove write permissions 
on everything inside the `data` directory. 
Run `chmod -R +w data/` to undo this.

- Protect your files

```
$ mkdir dont_delete_me
$ touch dont_delete_me/please_let_me_live.txt
$ chmod -R -w dont_delete_me/
$ rm -rf dont_delete_me
$ chmod -R +w dont_delete_me/
$ rm -rf dont_delete_me
```

??? answer
    Only the first `rm` gives output:
    ```
    rm: cannot remove ‘dont_delete_me/please_let_me_live.txt’: Permission denied
    ```

## Aliases on your desktop

1. Be logged in with ThinLinc.
2. Open a terminal.
3. Run `cd Desktop\`
4. Make shortcuts to your heart's content:
  - `PROJ=sens2023531`
  - `ln -s /proj/$PROJ/ proj`
  - `ln -s /proj/$PROJ/nobackup nobackup`
  - `ln -s /proj/$PROJ/nobackup/wharf/<yourusername>/<yourusername-$PROJ> wharf`

You can also make aliases to executables, 
like the convenient `interactive`-job starting script in `proj/useful_files`.
