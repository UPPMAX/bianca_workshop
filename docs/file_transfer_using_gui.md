# File transfer using a graphical tool

![](./img/filezilla_login_to_bianca_236_x_266.png)

> FileZilla connected to Bianca

## Overview

As a user, we need to transfer files
between our local computer and Bianca.
Here we learn how to transfer files
to Bianca using a graphical tool/program.

There are constraints on which programs
we can use, due to Bianca being an HPC cluster
for sensitive data.
Details are described in 'Bianca's constraints'.

When using such a graphical tool, 
one needs to be inside of SUNET.
FileZilla is a tool that is easy to setup.
The full procedure is described in 'Using FileZilla'.

The files you transfer will end up in your so-called `wharf` folder.
Where to find this folder is shown in the section 'Where do my files end up?'.

## Bianca's constraints

```mermaid
flowchart TD

    %% Give a white background to all nodes, instead of a transparent one
    classDef node fill:#fff,color:#000,stroke:#000

    %% Graph nodes for files and calculations
    classDef file_node fill:#faf,color:#000,stroke:#f0f
    classDef calculation_node fill:#aaf,color:#000,stroke:#00f

    subgraph sunet[Inside SUNET]
      subgraph bianca_outside[Bianca outside]
        login_node(login node):::calculation_node
        wharf(wharf):::file_node
        subgraph bianca_inside[Bianca]
          calculation_node(calculation/interative node):::calculation_node
          your_project_folder(Your project folder):::file_node
        end
      end
      user(User)
      sftp_client(SFTP client, e.g. FileZilla)
      sftp_server(SFTP server):::file_node
    end
    style bianca_inside fill:#afa,color:#000,stroke:#faa
    style bianca_outside fill:#ffa,color:#000,stroke:#faa
    style sunet fill:#faa,color:#000,stroke:#faa

    user --> |logs in |login_node
    user --> |uses| sftp_client
    sftp_client --> |transfer files|sftp_server
    sftp_server <--> |transfer files|wharf
    login_node --> |submit jobs|calculation_node
    login_node --> |can use|your_project_folder
    calculation_node --> |can use|your_project_folder
    wharf --> |transfer files| your_project_folder
```

Bianca is an HPC cluster for sensitive data.
To protect that sensitive data,
Bianca has no direct internet connection.
This means that files cannot be downloaded directly.

???- tip "What is an HPC cluster again?"

    What an HPC cluster is, is described in general terms [here](overview.md).

Instead, one needs to learn one of the many ways to do **secure** file transfer.

Here, we show how to transfer files using a graphical tool called FileZilla.

In general, one can pick any graphical tools with these constraints:

- the tool must support SFTP
- the tool must not store a password

Whatever tool one picks, it must do secure file transfer.
For secure file transfer, Bianca supports the SFTP protocol.
So, for secure file transfer to Bianca, one needs a tool
that supports SFTP.

???- warning "Use SFTP ... and why SCP will never work"

    You must use SFTP. 

    However, some users find tools that support
    another protocol called 'SCP'. 
    We understand the confusion, due
    to the many technical and abbreviated terms. 

    Only SFTP will work.

Whatever tool one picks, additionally, the tool must **not** store a password.
Due to security reasons, one needs to connect to Bianca using a password
**and** a two-factor authentication number (e.g. `VerySecret123456`).
If a tool stores a password, that password will be valid for only one session.

One tool that can be used for file transfer to Bianca
is FileZilla, which is described in detail below.
The extra materials at the bottom of this page contain
other tools.

## Using FileZilla

![FileZilla logo, from https://en.wikipedia.org/wiki/FileZilla#/media/File:FileZilla_logo.svg](./img/filezilla_logo_240_x_240.png)

> The FileZilla logo

FileZilla is a secure file transfer tool that works under Linux, Mac and Windows.

To transfer files to/from Bianca using FileZilla, do:

- Get inside SUNET

???- tip "Forgot how to get inside SUNET?"

    It is discussed [at the 'login to Bianca' page](login_bianca.md). 

- Start FileZilla
- From the menu, select 'File | Site manager'

???- tip "Where is that?"

    It is here:

    ![](./img/filezilla_file_site_manager.png)
    
    > The FileZilla 'File' menu contains the item 'Site manager'

- Click 'New site'

???- tip "Where is that?"

    It is here:

    ![](./img/filezilla_site_manager.png)

- Create a name for the site, e.g. `bianca-sens123456`.
- For that site, use all standards, except:
    - Set protocol to 'SFTP - SSH File Transfer Protocol'
    - Set host to `bianca-sftp.uppmax.uu.se`
    - Set user to `[username]-[project]`, e.g. `richel-sens123456`

???- tip "How does that look like?"

    It looks similar to these:

    ![](./img/filezilla_setup_bianca_pavlin.png)

    ![](./img/filezilla_setup_bianca_richel.png)

!!! tip "Storing a password is useless"

    Because Bianca holds sensitive data, 
    there is need to use the UPPMAX two-factor authentication
    code every time you login.
    Due to this, storing a password is hence useless

- Click 'Connect'
- You will be asked for your password with two-factor identification, hence
  type `[your password][2FA code]`, e.g. `VerySecret123456`

Now you can transfer files between your local computer and Bianca.

???- tip "How does that look like?"

    It looks like this:

    ![](./img/filezilla_login_to_bianca.png)

## Where do my files end up?

They end up in your personal `wharf` folder.

Its location is at `/home/[user_name]/[project_name]/nobackup/wharf/[user_name]/[user_name]-[project_name]`,
for example, at `/home/sven/sens123456/nobackup/wharf/sven/sven-sens123456`.

???- tip "How does that look like?"

    It looks like this:

    ![](./img/filezilla_file_on_bianca.png)

## Exercises

???- question "1. Your colleague uses an SCP tool for secure file transfer. It does not work on Bianca. What is the problem?"

    The assumptions of your colleague are already wrong:
    SCP is considered an **outdated** secure (`outdated + secure = insecure`) 
    file transfer protocol.

    Bianca does not support this outdated file transfer protocol.
    Instead, use a tool that allows the SFTP secure file transfer protocol.

???- question "2. Your colleague used a very user-friendly SFTP tool. However, only on Bianca, it fails after the first use. What is the problem?"

    This user-friendly program stores the password of the user.
    However, the password of Bianca changes, due to the need
    of adding a two-factor authentication number at the end
    of a password (e.g. `VerySecret123456`).

    Your colleague can either reset the password every time he/she
    uses the program, or use another tool (such as FileZilla) instead.

???- question "3a. Exercise: upload a file to Bianca"

    This is described in the text above.

    There is also a video [here](https://youtu.be/V-iPQLjvByc).

???- question "3b. Locate that file on Bianca using the graphical file explorer"

    Its location is at `/home/[user_name]/[project_name]/nobackup/wharf/[user_name]/[user_name]-[project_name]`,
    for example, at `/home/sven/sens123456/nobackup/wharf/sven/sven-sens123456`.

    Start the file explorer and navigate there.

    ![](./img/filezilla_file_on_bianca.png)

???- question "3c. Locate that file on Bianca using the terminal"

    Its location is at `/home/[user_name]/[project_name]/nobackup/wharf/[user_name]/[user_name]-[project_name]`,
    for example, at `/home/sven/sens123456/nobackup/wharf/sven/sven-sens123456`.

    Start a terminal and use `cd` to change to that directory, e.g., 
    from your home folder, run
    `cd sens123456/nobackup/wharf/sven/sven-sens123456`.

    In that folder, type `ls` to list the files present in that folder.

    ![](./img/filezilla_file_on_bianca.png)

## Extra material

### WinSCP

WinSCP is a secure file transfer tool that works under Windows.

To transfer files to/from Bianca using WinSCP, do:

- Get inside SUNET

???- tip "Forgot how to get inside SUNET?"

    It is discussed [at the 'login to Bianca' page](login_bianca.md). 

- Start WinSCP
- Create a new site
- For that site, use all standards, except:
    - Set file protocol to 'SFTP'
    - Set host name to `bianca-sftp.uppmax.uu.se`
    - Set user name to `[username]-[project]`, e.g. `richel-sens123456`

???- tip "How does that look like?"

    It looks like this:

    ![WinSCP](./img/winscp-snaphot1.png)

### File transfer overview

```mermaid
flowchart TD

    %% Give a white background to all nodes, instead of a transparent one
    classDef node fill:#fff,color:#000,stroke:#000

    %% Graph nodes for files and calculations
    classDef file_node fill:#faf,color:#000,stroke:#f0f
    classDef calculation_node fill:#aaf,color:#000,stroke:#00f

    subgraph sunet[Inside SUNET]
      subgraph bianca_outside[Bianca outside]
        login_node(login node):::calculation_node
        wharf(wharf):::file_node
        subgraph bianca_inside[Bianca]
          calculation_node(calculation node):::calculation_node
          interactive_node(interactive node):::calculation_node
          your_project_folder(Your project folder):::file_node
        end
      end
      user(User)
      sftp_client(SFTP client, e.g. FileZilla)
      transit(transit):::file_node
      sftp_server(SFTP server):::file_node
      terminal(terminal)
    end
    style bianca_inside fill:#afa,color:#000,stroke:#faa
    style bianca_outside fill:#ffa,color:#000,stroke:#faa
    style sunet fill:#faa,color:#000,stroke:#faa

    user --> |logs in |login_node
    user --> |uses| sftp_client
    user --> |uses| terminal
    sftp_client --> |transfer files|sftp_server
    sftp_server <--> |transfer files|wharf
    terminal --> |transfer files| transit
    transit <--> |transfer files|wharf
    login_node --> |submit jobs|calculation_node
    login_node --> |login to|interactive_node
    interactive_node --> |submit jobs|calculation_node
    login_node --> |can use|your_project_folder
    interactive_node --> |can use|your_project_folder
    calculation_node --> |can use|your_project_folder
    wharf --> |transfer files| your_project_folder
```
