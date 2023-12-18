# Introduction to Bianca: Handling Sensitive Research Data

Are you starting to work with your sensitive data in your research? 

If yes, welcome to a full day introduction to handling sensitive data on the UPPMAX cluster, Bianca!

You will learn about the national infrastructure Bianca is part of, how to login to Bianca, upload and download files, using pre-installed software and how to start your code.

Tentative schedule:

Time |Topic
-----|------------------------
9:00 |[Introduction](intro.md)
9:15 |[NAISS-Sens](sens_project_short.md)
9:35 |[Login](login_bianca.md)
9:55 |Break
10:10|[Command line](commandline.md)
10:55|Break
11:10|[Modules](modules1.md)
12:00|Lunch
13:00|[Transferring files to and from Bianca](transfer_basics.md)
13:55|Break
14:10|[Compute nodes and SLURM](slurm_intro.md)
14:55|Break
15:10|Summary
15:15|Q/A
15:55|Closing words
16:00|END

!!! info "Q/A collaboration document"

    - Use the Q/A page for the workshop with your questions.

          - [https://hackmd.io/@bclaremar/bianca_dec_2023?both](https://hackmd.io/@bclaremar/bianca_dec_2023?both)

    - Depending on how many helpers there are we’ll see how fast there are answers.

        - Some answers may come after the workshop.

    - Create a new line for new questions. Take care if others are editing at the same time.

## Overview of courses

%% Direction is top->down
flowchart TD

    %% Give a white background, instead of a transparent one
    classDef node fill:#fff,color:#000,stroke:#000
    classDef focus_node fill:#fff,color:#000,stroke:#000,stroke-width:4px
    
    subgraph sub_basic_use["Basic use of Bianca"]
      can_login_to_remove_desktop(Can login to remote deskop)
      can_use_command_line_1(Can use the command line 1)
      can_create_bash_script_using_gui(Can create a bash script using GUI)
      can_use_modules(Can use modules)
      can_use_interactive_node(Can use an interactive node)
      can_develop_code_interactively(Can develop code interactively)
      can_tranfer_files_using_gui(Can transfer files using GUI):::focus_node
      can_schedule_jobs(Can schedule jobs):::focus_node
      can_use_ide(Can use an IDE):::focus_node
    end
    style sub_basic_use fill:#ffa,color:#000,stroke:#ffa



    subgraph sub_intermediate_use["Intermediate use of Bianca"]
      can_login_to_console(Can login to console)
      can_use_command_line_2(Can use the command line 2)
      can_create_bash_script_using_cli(Can create a bash script using CLI)
      can_tranfer_files_using_cli(Can transfer files using a command-line tool):::focus_node
      can_use_custom_software(Can use custom software):::focus_node
      can_monitor_jobs(Can monitor jobs):::focus_node
    end
    style sub_intermediate_use fill:#afa,color:#000,stroke:#afa

    %% Basic
    can_login_to_remove_desktop ---> can_use_command_line_1
    can_login_to_remove_desktop ---> can_tranfer_files_using_gui
    can_login_to_remove_desktop ---> can_create_bash_script_using_gui
    can_use_command_line_1 --> can_use_modules
    can_use_command_line_1 --> can_use_interactive_node
    can_use_command_line_1 --> can_use_command_line_2
    can_use_command_line_1 --> can_schedule_jobs
    can_develop_code_interactively --> can_use_ide
    can_use_modules --> can_schedule_jobs
    can_use_interactive_node --> can_use_ide
    can_use_interactive_node --> can_develop_code_interactively
    can_create_bash_script_using_gui --> can_schedule_jobs

    %% Basic -> Intermediate
    can_tranfer_files_using_gui --> can_use_custom_software
    can_schedule_jobs --> can_monitor_jobs

    %% Make sure Intermediat in below Basic,
    %% using invisible nodes
    can_use_ide ~~~ can_use_command_line_2

    %% Intermediate
    can_login_to_console --> can_tranfer_files_using_cli
    can_use_command_line_2 --> can_create_bash_script_using_cli
    can_use_command_line_2 --> can_tranfer_files_using_cli
    can_tranfer_files_using_cli --> can_use_custom_software
```
    
