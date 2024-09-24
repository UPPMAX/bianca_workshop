# Command-line intro

!!! info "Objectives"
    - We'll use the commands and investigate the Bianca environment

!!! warning

    - We assume that you have already covered the Command-line material and tested on Rackham
        - [LINUX](https://uppmax.github.io/uppmax_intro/linux.html)
        - [Basic toolkit](https://uppmax.github.io/uppmax_intro/linux_basics.html)
        
### Navigation and file management

1. `pwd`  &emsp; present directory
1. `ls`  &emsp;list content
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

