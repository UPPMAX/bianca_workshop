# Introduction to Bianca: Handling Sensitive Research Data

???- question "Prefer a video?"

    In case you prefer a video over reading,
    [watch this YouTube video (6 minutes long)](https://youtu.be/o0fRHoa8C5U)

Are you starting to work with your sensitive data in your research? 

If yes, welcome to a full day introduction to handling sensitive data on the UPPMAX cluster, Bianca!

You will learn about the national infrastructure Bianca is part of, how to login to Bianca, upload and download files, using pre-installed software and how to start your code.

Tentative schedule:

When  | Who  | What
------|------|-----------------------------
9:00  | R    | [Introduction](intro.md) and [logging in](login_bianca.md)
10:00 | .    | Break
10:15 | R    | [Command line](commandline.md)
11:00 | .    | Break
11:15 | R    | [Modules](modules1.md)
12:00 | .    | Lunch
13:00 | P    | [Transferring files to and from Bianca](transfer_basics.md)
13:45 | .    | Break
14:00 | B    | [Compute nodes and SLURM](slurm_intro.md), including interactive nodes
14:45 | .    | Break
15:00 | L    | [Legal aspects of sensitive data](sens_project_short.md)
15:30 | R    | Summary
15:35 | R    | Anonymous evaluation
15:45 | All  |Optional Q&A

 * Who: `B`: Björn, `L`: Lars, `P`: Pavlin, `R`: Richèl

!!! info "Q/A collaboration document"

    - Use the Q/A page for the workshop with your questions.

          - [https://hackmd.io/@UPPMAX/Bianca_Intro_QaA](https://hackmd.io/@UPPMAX/Bianca_Intro_QaA)

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
