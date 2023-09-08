# The command-line on Bianca

!!! info "Objectives"
    - Being able to navigate in/out folder
    - Being able to view/create/move/delete files
    - Create an executable bash script

## Exercises

Tips: 

 * You can Google all this
 * At the bottom of this page are the solutions.

 1. View the help of the command `cd`
 2. Go to the folder `/proj/sens2023598/workshop/`
 3. Create a folder `/proj/sens2023598/workshop/[your_login_name]`,
    for example, `/proj/sens2023598/workshop/richel`
 3. In that folder, create an executable script called `do_it.sh`.
 4. Upon running the script, it should:
   * display a welcome message in text, e.g. `Hello!`
   * do something, e.g. show the files in reverse order
 5. Copy the file (e.g. to `do_it_again.sh`). 
 6. Move the copied file (e.g. move it one folder up to `../do_it_again_here.sh`)
 7. Delete the copied file

![Using the command line on a computer cluster](./img/610803_a_woman_using_the_command_line_on_a_computer_cluster.png)

> Using the command line on a computer cluster ![./img/](public_domain_88x31.png)

## Command-line intro

### Help

1. `man`  &emsp;info about a command
      
### Navigation and file management

1. `pwd`  &emsp; present directory
1. `ls`  &emsp;list content
1. `mkdir`  &emsp;make directory
1. `cd`  &emsp;change directory
1. `rmdir`  &emsp;remove empty directory
1. `touch` &emsp;create an empty file
1. `cp`  &emsp;copy
1. `mv`  &emsp;move
1. `rm`  &emsp;remove
1. `scp`  &emsp;securely remotely copy

### Read files and change file properties

1. `cat`  &emsp;print content on screen
1. `head`  &emsp;print first part
1. `tail`  &emsp;print last part
1. `less`  &emsp;browse content
1. `chmod`  &emsp;change file permissions
1. `echo`  &emsp;print text on screen

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
$ cd /proj/sensXXX
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
  - `ln -s /proj/$PROJ/nobackup/wharf/${USER}/${USER}-$PROJ wharf`

## Solutions

### 1. View the help of the command `cd`

```
man cd
```

### 2. Go to the folder `/proj/sens2023598/workshop/`

```
cd /proj/sens2023598/workshop/
```

### 3. Create a folder `/proj/sens2023598/workshop/[your_login_name]`

```
cd /proj/sens2023598/workshop/
mkdir richel
```
### 4. In that folder, create an executable script called `do_it.sh`.

```
touch do_it.sh
chmod +x do_it.sh
```

It does nothing, but it is executable :-)

### 5. Upon running the script, it should do things

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

### 6. Copy the file (e.g. to `do_it_again.sh`)

```
copy do_it.sh do_it_again.sh
```

### 7. Move the copied file


```
mv do_it_again.sh ../do_it_again_richel.sh
```

### 8. Delete the copied file

```
rm ../do_it_again_richel.sh
```

or

```
cd ..
rm do_it_again_richel.sh
```



## Links

 * Video 'Using the command-line on the UPPMAX Bianca cluster': [YouTube](https://youtu.be/kjqLAx2bgJI), [download (.ogv)](https://richelbilderbeek.nl/uppmax_bianca_command_line.ogv)
 * [UPPMAX intro course materials on LINUX](https://uppmax.github.io/uppmax_intro/linux.html)
