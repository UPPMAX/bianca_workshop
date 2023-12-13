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

    9.00 Introduction
    9.10 NAISS-SENS summary
    9.20 Transferring files to and from Bianca
    10.00 Coffee break
    10.15 Transferring files p2
    10.35 Slurm jobs p1
    11.00 Break
    11.15 Slurm jobs p2 
    12.00 LUNCH
    13.00 Software and packages installation on Bianca
    13.50 break
    14.05 IDE:s on Bianca
    15.00 Coffee break
    15.15 Summary
    15.20 Q/A and extra material
    15.55 Closing words
    16.00 END

!!! info "Q/A collaboration document"

    - Use the Q/A page for the workshop with your questions.

          - [https://hackmd.io/@bclaremar/bianca_hack_dec_2023?both](https://hackmd.io/@bclaremar/bianca_dec_2023?both)

    - Depending on how many helpers there are we’ll see how fast there are answers.

        - Some answers may come after the workshop.

    - Create a new line for new questions. Take care if others are editing at the same time.

## Overview of courses

```mermaid
flowchart TD

    %% Give a white background, instead of a transparent one
    classDef node fill:#fff,color:#000,stroke:#000

    subgraph sub_prerequisites["Preprequisites to use Bianca"]
      can_login(Can login)
      can_use_command_line(Can use the command line)
      can_create_bash_script(Can create a bash script)
      can_use_modules(Can use modules)
      can_use_interactive_node(Can use an interactive node)
    end
    style sub_prerequisites fill:#f00,color:#000,stroke:#faa

    subgraph sub_basic_use["Basic use of Bianca"]
      can_develop_code_interactively(Can develop code interactively)
      can_tranfer_files_using_gui(Can transfer files using graphical user interface)
      can_tranfer_files_using_cli(Can transfer files using a command-line tool)
      can_schedule_jobs(Can schedule jobs)
    end
    style sub_basic_use fill:#ff0,color:#000,stroke:#ffa

    subgraph sub_intermediate_use["Intermediate use of Bianca"]
      can_use_custom_software(Can use custom software)
      can_monitor_jobs(Can monitor jobs)
      %% Richel: I think this is basic use, as it is beginners that want this
      can_use_ide(Can use an IDE)
    end
    style sub_intermediate_use fill:#0f0,color:#000,stroke:#afa

    can_login ---> can_use_command_line
    can_login ---> can_tranfer_files_using_gui
    can_use_command_line --> can_create_bash_script
    can_use_command_line --> can_use_modules
    can_use_command_line --> can_use_interactive_node
    can_use_command_line --> can_tranfer_files_using_cli
    can_use_modules --> can_schedule_jobs
    can_create_bash_script --> can_schedule_jobs
    can_schedule_jobs --> can_monitor_jobs
    can_use_interactive_node --> can_use_ide
    can_use_interactive_node --> can_develop_code_interactively
    can_tranfer_files_using_cli --> can_use_custom_software
    can_tranfer_files_using_gui --> can_use_custom_software
```
    
