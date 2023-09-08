# The command-line on Bianca

!!! info "Objectives"
    - Being able to navigate in/out folders
    - Being able to view/create/move/delete files
    - Create an executable bash script

## Exercises

 * 1a. View the help of the command `cd`
 * 2a. Navigate to the project folder, e.g. `/proj/sens2023598`
 * 2b. Navigate to your home folder
 * 2c. Navigate to the wharf, e.g. `/proj/sens2023598/nobackup/wharf`
 * 3a. Create a folder `/proj/sens2023598/workshop/[your_login_name]`, for example, `/proj/sens2023598/workshop/richel`
 * 4a. Create a file, e.g. `richel.txt`
 * 4b. Copy the file (e.g. to `richel_again.txt`). 
 * 4c. Move the copied file (e.g. move it one folder up to `../richel_again.txt`)
 * 4d. Delete the copied file
 * 5a. Create an executable script called `/proj/sens2023598/workshop/[your_login_name]/do_it.sh`,
   which, upon running, displays a welcome message in text (e.g. `Hello!`)
   and does something (e.g. show the files in reverse order)

![Using the command line on a computer cluster](./img/610803_a_woman_using_the_command_line_on_a_computer_cluster.png)

> Using the command line on a computer cluster ![Public domain](./img/public_domain_88x31.png)

## 1. Help

Use `man` to see the help pages about a command:

```
man cd
man man
man ls
```

Press CTRL-Z to exit `man`

## 2. Navigation

Use `cd` to change directory/folder:

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

 * Tip: use `ls` to see the content of a folder
 * Tip: use `pwd` to see your current location

!!! info "The Silence Is Golden Rule"
    When your command 'just works' there is no output
    (try, for example `cd ~`). 
    This is due to [The Silence Is Golden Rule](https://www.linfo.org/rule_of_silence.html)

## 3. Working with folders

Do what                            |Example command
-----------------------------------|---------------------
Create a folder                    |`mkdir myfolder`
Delete an empty folder             |`rmdir myfolder`

 * Tip: use `ls` to see the content of a folder
 * Tip: use `pwd` to see your current location
 * Tip for sysadmins: use `pwd -P` to see your real current location on the hardware

## 4. Working with files

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

## 5. Creating an executable script

Creating an executable script has two steps:

 * 1. Create a script
 * 2. Allow the script to execute

As an example, we create a script, called `do_it.sh`:

```
nano do_it.sh
```

Use `.sh` as a file extension a social convention for how a Bash script is called,
as (1) `sh` is short for 'shell', (2) Bash is short for 'Bourne Again Shell'.
A 'shell' in this context is a program that allows working with an operating system. 

As an example, copy-paste this content into the script:

```
#!/bin/bash
echo "Hello!"
ls | rev
```

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

## 6. Other useful commands

These are some commands that we enjoy.

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

## 7. The terminal and the GUI are friends

On a clean terminal, try typing `cd` 
and then dragging a folder from the GUI to the terminal.

It types the absolute path for you!

## 8. Commonly used links

These are some commonly used links:

```
cd Desktop
ln -s /proj/sens2023598/ proj
ln -s /proj/sens2023598/nobackup nobackup
ln -s /proj/sens2023598/nobackup/wharf/richel/richel-sens2023598 wharf`
```

 * Replace `sens2023598` by your project
 * Replace `richel` by your username

## Solutions

### 1a. View the help of the command `cd`

```
man cd
```

### 2a. Go to the project folder, e.g. `/proj/sens2023598`

```
cd /proj/sens2023598
```

Don't forget the `/` at the start.

### 2b. Go to your home folder

 * Your home folder

```
cd /home/richel
```


The squiggle/tilde (`~`) is a shorter notation:

```
cd ~
```

### 2c. Go to the the wharf, e.g. `/proj/sens2023598/nobackup/wharf`

```
cd /proj/sens2023598/nobackup/wharf
```

### 3a. Create a folder `/proj/sens2023598/workshop/[your_login_name]`

```
mkdir /proj/sens2023598/workshop/richel
```

Or navigate there first:

```
cd /proj/sens2023598/workshop/
mkdir richel
```

### 4a. Create a file, e.g. `richel.txt`

```
touch richel.txt
```

### 4b. Copy the file (e.g. to `richel_again.txt`). 

```
cp richel.txt richel_again.txt
```

### 4c. Move the copied file (e.g. move it one folder up to `../richel_again.txt`)

```
mv richel_again.txt ../
```
### 4d. Delete the copied file

```
rm ../richel_again.txt
```

or:

```
cd ..
rm richel_again.txt
```

### 5. In that folder, create an executable script called `do_it.sh`.

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


## Links

 * Video 'Using the command-line on the UPPMAX Bianca cluster': [YouTube](https://youtu.be/kjqLAx2bgJI), [download (.ogv)](https://richelbilderbeek.nl/uppmax_bianca_command_line.ogv)
 * [UPPMAX intro course materials on LINUX](https://uppmax.github.io/uppmax_intro/linux.html)
