# Introduction to Bianca: Handling Sensitive Research Data

![Bianca Castafiore](./img/bianca_castafiore_192_x_226.png)

> Bianca Castafiore, the Tintin character the cluster is named after.

???- question "Prefer a video?"

    In case you prefer a video over reading,
    [watch this YouTube video (6 minutes long)](https://youtu.be/o0fRHoa8C5U)

???- info "Notes for teachers"

    Teaching goals:

    - The learners have heard the topics of this course

    Schedule:

    ```mermaid
    gantt
      title Lesson plan Introduction and Logging in
      dateFormat X
      axisFormat %s
      section Introduction
      Prior knowledge: prior_1, 0, 5s
      Theory : theory_1, after prior_1, 5s
      section Logging In
      Prior knowledge: prior_2, after theory_1, 5s
      Theory: theory, after prior_2, 5s
      Exercises: crit, exercise, after theory, 30s
      Feedback: feedback, after exercise, 10s
    ```

## Introduction

Are you starting to work with your sensitive data in your research?

If yes, welcome to a full day introduction to
handling sensitive data on the UPPMAX cluster, Bianca!

You will learn how to login to Bianca, upload and download files,
using pre-installed software and how to run your code.

## Schedule

When  | Who  | What
------|------|-----------------------------
9:00  | R    | [Introduction](intro.md) (this page), [login](login.md), [log in to remote desktop](login_remote_desktop.md)
10:00 | .    | Break
10:15 | R    | [Use the remote desktop](use_remote_desktop.md), [log in to console environment](login_console.md), [command line](commandline.md)
11:00 | .    | Break
11:15 | R    | [Create an executable script](create_script.md), [use modules](modules.md)
12:00 | .    | Lunch
13:00 | P    | Transferring files to and from Bianca 55m
13:55 | .    | Break 10m
14:05 | B    | Compute nodes and SLURM, including interactive nodes 50m
14:55 | .    | Break
15:10 | L    | Legal aspects of sensitive data 20m
15:30 | R    | Summary
15:35 | R    | Anonymous evaluation
15:45 | All  | Optional Q&A

- Who: `B`: Björn, `L`: Lars, `P`: Pavlin, `R`: Richèl

!!! info "Q/A collaboration document"

    - Use the Q/A page for the workshop with your questions.
        - [https://hackmd.io/@UPPMAX/Bianca_Intro_QaA](https://hackmd.io/@UPPMAX/Bianca_Intro_QaA){target="_blank"}
    - Depending on how many helpers there are,
      we'll see how fast there are answers
        - Some answers may come after the workshop.
    - Create a new line for new questions.
      Take care if others are editing at the same time.

## Overview of courses

```mermaid
%% Direction is top->down
flowchart TD

    %% Give a white background, instead of a transparent one
    classDef node fill:#fff,color:#000,stroke:#000
    classDef focus_node fill:#fff,color:#000,stroke:#000,stroke-width:4px
    classDef basic_node fill:#fdd,color:#000,stroke:#f00
    classDef basic_extra_node fill:#ffd,color:#000,stroke:#ff0
    classDef intermediate_node fill:#dfd,color:#000,stroke:#0f0
  
    %% subgraph sub_basic_use[Basic use of Bianca]
      understand_login(Understand login, has 2FA):::basic_node
      can_login_to_remove_desktop(Can login to remote deskop):::basic_node
      can_login_to_console(Can login to console):::basic_extra_node
      can_navigate_filesystem_using_gui(Can navigate filesystem using GUI):::basic_node
      can_navigate_filesystem_using_cli(Can navigate filesystem using CLI):::basic_extra_node
      can_find_wharf(Can find the wharf):::basic_node
      can_use_command_line_1(Can use the command line 1):::basic_node
      can_use_modules(Can use modules):::basic_node
      can_use_interactive_node(Can use an interactive node):::basic_node
      can_manage_files_using_cli(Can manage files using CLI):::basic_extra_node
      can_schedule_jobs(Can schedule jobs):::basic_node
      can_create_bash_script_using_cli(Can create a bash script using CLI):::basic_node
      can_tranfer_files_using_gui(Can transfer files using GUI):::basic_node
    %% end
    %% style sub_basic_use fill:#fcc,color:#000,stroke:#fcc

    %% subgraph sub_intermediate_use[Intermediate use of Bianca]
      can_tranfer_files_using_rsync(Can transfer files using rsync):::intermediate_node
      can_use_command_line_2(Can use the command line 2):::intermediate_node
      can_use_custom_software(Can use custom software):::intermediate_node
      can_monitor_jobs(Can monitor jobs):::intermediate_node
      can_use_ide(Can use an IDE):::intermediate_node
    %% end
    %% style sub_intermediate_use fill:#ffc,color:#000,stroke:#ffc


    %% Basic
    understand_login --> can_login_to_remove_desktop
    can_login_to_remove_desktop --> can_login_to_console
    can_login_to_remove_desktop --> can_use_command_line_1
    can_login_to_remove_desktop --> can_navigate_filesystem_using_gui
    can_navigate_filesystem_using_gui --> can_find_wharf
    can_login_to_console --> can_navigate_filesystem_using_cli
    can_navigate_filesystem_using_cli --> can_find_wharf
    can_find_wharf --> can_tranfer_files_using_gui
    understand_login --> can_login_to_console
    can_login_to_console --> can_use_command_line_1
    can_use_command_line_1 --> can_use_modules
    can_use_command_line_1 --> can_use_interactive_node
    can_use_command_line_1 --> can_use_command_line_2
    can_use_command_line_1 --> can_create_bash_script_using_cli
    can_use_command_line_1 --> can_schedule_jobs
    can_navigate_filesystem_using_cli --> can_manage_files_using_cli
    can_use_modules --> can_schedule_jobs
    can_create_bash_script_using_cli --> can_schedule_jobs

    %% Basic -> Intermediate
    can_manage_files_using_cli --> can_tranfer_files_using_rsync
    can_schedule_jobs --> can_monitor_jobs
    can_use_interactive_node --> can_use_ide
    can_use_modules --> can_use_ide

    %% Make sure Intermediate is below Basic,
    %% using invisible nodes
    can_schedule_jobs ~~~ can_use_command_line_2

    %% Intermediate
    can_use_command_line_1 --> can_tranfer_files_using_rsync
    can_use_command_line_2 --> can_use_custom_software
```

> Overview of the courses.
> Red nodes: Intro to Bianca.
> Yellow node: Intro to Bianca extra material.
> Green node: Intermediate Bianca.
