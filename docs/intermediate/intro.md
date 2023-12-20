# Bianca In-Depth Hackathon: Improve Your Handling of Sensitive Research Data

Are you already working with sensitive data in your research 
and feel your workflow can be improved? 
If yes, welcome to a full day of learning of smarter ways 
to work on the Bianca UPPMAX cluster. 

???- tip "I am new to Bianca, where do I start?"

    As a beginner, you are encouraged to start with the
    basic Bianca course, at [here](../intro.md).

You will learn how to login from a terminal (bypassing ThinLinc), 
do file transfer from a terminal (allowing scripts), 
advanced SLURM, using IDEs (i.e. RStudio and/or VSCode), 
and installing custom software and packages.

Tentative schedule

Time |Topic
-----|--------------------------
9.00| [Introduction](https://uppmax.github.io/bianca_workshop/intermediate/intro/)
9.10| [NAISS-SENS summary](https://uppmax.github.io/bianca_workshop/sens_project_short/)
9.20| [Transferring files to and from Bianca](https://uppmax.github.io/bianca_workshop/intermediate/transfer/)
10.00| Coffee break
10.15| [Transferring files p2](https://uppmax.github.io/bianca_workshop/intermediate/transfer/)
10.35| [Slurm jobs p1](https://uppmax.github.io/bianca_workshop/intermediate/slurm_intermed/)
11.00| Break
11.15| [Slurm jobs p2](https://uppmax.github.io/bianca_workshop/intermediate/slurm_intermed/)
12.00| LUNCH
13.00| [Software and packages installation](https://uppmax.github.io/bianca_workshop/intermediate/install/)
13.50| break
14.05| [IDE:s on Bianca](https://uppmax.github.io/bianca_workshop/intermediate/ides/)
15.00| Coffee break
15.15| [Summary](https://uppmax.github.io/bianca_workshop/intermediate/summary/)
15.20| Q/A and [extra material](https://uppmax.github.io/bianca_workshop/extra/slurm/)
15.55| Closing words
16.00| END

!!! info "Q/A collaboration document"

    - Use the Q/A page for the workshop with your questions.

          - [https://hackmd.io/@bclaremar/bianca_hack_dec_2023?both](https://hackmd.io/@bclaremar/bianca_dec_2023?both)

    - Depending on how many helpers there are we’ll see how fast there are answers.

        - Some answers may come after the workshop.

    - Create a new line for new questions. Take care if others are editing at the same time.

## Overview of courses

```mermaid
%% Direction is top->down
flowchart TD

    %% Give a white background, instead of a transparent one
    classDef node fill:#fff,color:#000,stroke:#000
    classDef focus_node fill:#fff,color:#000,stroke:#000,stroke-width:4px
    
    subgraph sub_basic_use["Basic use of Bianca"]
      can_login_to_remove_desktop(Can login to remote deskop)
      can_login_to_console(Can login to console)
      can_use_command_line_1(Can use the command line 1)
      can_use_modules(Can use modules)
      can_use_interactive_node(Can use an interactive node):::focus_node
      can_manage_files_using_cli(Can manage files using CLI)
      can_tranfer_files_using_rsync(Can transfer files using rsync):::focus_node
      can_schedule_jobs(Can schedule jobs):::focus_node
      can_create_bash_script_using_cli(Can create a bash script using CLI)
    end
    style sub_basic_use fill:#fcc,color:#000,stroke:#fcc

    subgraph sub_intermediate_use["Intermediate use of Bianca"]
      can_use_command_line_2(Can use the command line 2)
      can_use_custom_software(Can use custom software):::focus_node
      can_monitor_jobs(Can monitor jobs):::focus_node
      can_use_ide(Can use an IDE)
    end
    style sub_intermediate_use fill:#ffc,color:#000,stroke:#ffc

    subgraph sub_non_goal["Not in course"]
      can_tranfer_files_using_gui(Can transfer files using GUI)
      can_create_bash_script_using_gui(Can create a bash script using GUI)
    end
    style sub_non_goal fill:#ccc,color:#000,stroke:#ccc


    %% Basic
    can_login_to_console --> can_tranfer_files_using_rsync
    can_login_to_console --> can_use_command_line_1
    can_login_to_remove_desktop ---> can_use_command_line_1
    can_use_command_line_1 --> can_use_modules
    can_use_command_line_1 --> can_use_interactive_node
    can_use_command_line_1 --> can_use_command_line_2
    can_use_command_line_1 --> can_create_bash_script_using_cli
    can_use_command_line_1 --> can_schedule_jobs
    can_use_command_line_1 --> can_manage_files_using_cli
    can_use_command_line_1 --> can_tranfer_files_using_rsync
    can_use_modules --> can_schedule_jobs
    can_create_bash_script_using_cli --> can_schedule_jobs

    %% Basic -> Intermediate
    can_schedule_jobs --> can_monitor_jobs
    can_use_interactive_node --> can_use_ide
    can_use_modules --> can_use_ide

    %% Make sure Intermediate is below Basic,
    %% using invisible nodes
    can_schedule_jobs ~~~ can_use_command_line_2

    %% Intermediate
    can_use_command_line_2 --> can_use_custom_software

    %% Basic -> None
    %% can_login_to_remove_desktop ---> can_tranfer_files_using_gui
    %% can_tranfer_files_using_gui --> can_use_custom_software
    %% can_login_to_remove_desktop ---> can_create_bash_script_using_gui
    %% can_create_bash_script_using_gui --> can_schedule_jobs

    %% Make sure Non-goals is below Intermediat,
    %% using invisible nodes
    can_use_custom_software ~~~ can_tranfer_files_using_gui
```
