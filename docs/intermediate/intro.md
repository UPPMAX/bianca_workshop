# Bianca In-Depth Hackathon: Improve Your Handling of Sensitive Research Data

!!!- info "Learning objectives"

    - Know where to find an overview of topics discussed today
    - Know where to find the link to the shared document
    - Know where to find the schedule

???- question "For teachers"

    Prerequisites are:

    - Link to shared document

    Preparations are:

    - Shared document must contain:
        - a proper header with specifics: Zoom room, Bianca project, etc.
        - already redirect to evaluation at bottom
        - evaluation at the bottom

    Teaching goals are:

    - Learners know where to find an overview of topics discussed today
    - Learners know where to find the link to the shared document
    - Learners know where to find the schedule

    Lesson plan:

    ```mermaid
    gantt
      title Introduction
      dateFormat X
      axisFormat %s
      Prior: prior, 0, 5s
      Introduction: after prior, 5s
    ```

    Prior questions:

    - What do you think about when you think about Bianca?
    - How would you like to improve your Bianca experience?

Are you already working with sensitive data in your research
and feel your workflow can be improved?
If yes, welcome to a full day of learning of smarter ways
to work on the Bianca UPPMAX cluster.

???- tip "Help! I am new to Bianca, where do I start instead?"

    As a beginner, you are encouraged to start with the
    basic Bianca course, at [here](../intro.md).

???- tip "I may be rusty on how to log in to Bianca, how do I login again?"

    You can refresh your mind at 
    [the Basic Bianca course's session on accessing Bianca](../login_bianca.md).

???- tip "I may be rusty on using Linux commands, is there any documentation?"

    Sure! It is at
    [the Basic Bianca course's session on using the command-line](../commandline.md).

???- tip "How do I copy/paste within a ThinLinc session?" 

     - Windows/MAC: Right-click and choose, or
     - Windows:
         - paste: `shift+insert`
         - copy: `ctrl+insert`

???- tip "How do I copy/paste between ThinLinc and locally?" 

      ![copy-paste](../img/copy_paste.PNG)

      - Like a "wharf" for clipboard

      - Copy in ThinLinc by the ThinLinc command (see above) and it ends up here in the ThinLinc clipboard
          - Mark and copy with Windows/Mac command
          - Paste locally with Windows/Mac command
       
      - Copy from locally 
          - paste in the ThinLinc clipboard with Windows/Mac command
          - paste to ThinLinc place by the ThinLinc command (see above)

!!! info "Course shared document"

    The course shared document can be found at: [https://hackmd.io/@UPPMAX/bianca-workshop-2024/edit](https://hackmd.io/@UPPMAX/bianca-workshop-2024/edit).

    It contains exercises, room for questions and room for feedback.

    To view the shared document in both text and rendered form, 
    click on the button in the top between the pen and eye.

    Regarding questions: 

    - create a new line for new questions
    - take care if others are editing at the same time
    - the more teaching assistants available, the quicker questions are answered
    - some questions may be answered after the workshop

## Schedule

When | Who  | What
-----|------|-------------------------
9:00 | R    | [Introduction](https://uppmax.github.io/bianca_workshop/intermediate/intro/)
9:10 | R    | 50 minutes [transferring files to and from Bianca](https://uppmax.github.io/bianca_workshop/intermediate/transfer/), which includes login
10:00| .    | Break
10:15| R    | 20 minutes [Transferring files p2](https://uppmax.github.io/bianca_workshop/intermediate/transfer/)
10:35| D    | 25 minutes [Slurm jobs p1](https://uppmax.github.io/bianca_workshop/intermediate/slurm_intermed/)
11:00| .    | Break
11:15| D    | 45 minutes [Slurm jobs p2](https://uppmax.github.io/bianca_workshop/intermediate/slurm_intermed/)
12:00| .    | Lunch
13:00| B    | [Software and packages installation](https://uppmax.github.io/bianca_workshop/intermediate/install/)
13:45| .    | Break
14:00| BR   | [IDEs on Bianca](https://uppmax.github.io/bianca_workshop/intermediate/ides/)
14:45| .    | Break
15:00| L    | [NAISS-SENS summary](https://uppmax.github.io/bianca_workshop/sens_project_short/)
15:30| R    | [Summary](https://uppmax.github.io/bianca_workshop/intermediate/summary/)
15:35| R    | [Evaluation](https://uppmax.github.io/bianca_workshop/intermediate/summary/)
15:45| All  | Optional Q/A and [extra material](https://uppmax.github.io/bianca_workshop/extra/slurm/)
16:00| .    | END

* Who: `B`: Björn, `L`: Lars, `R`: Richèl

## Overview of courses

```mermaid
%% Direction is top->down
flowchart TD

    subgraph sub_basic_use["Basic use of Bianca"]
      can_login_to_remove_desktop(Can login to remote deskop)
      can_login_to_console(Can login to console)
      can_use_command_line_1(Can use the command line 1)
      can_use_modules(Can use modules)
      can_use_interactive_node(Can use an interactive node)
      can_manage_files_using_cli(Can manage files using CLI)
      can_tranfer_files_using_rsync(Can transfer files using rsync)
      can_schedule_jobs(Can schedule jobs)
      can_create_bash_script_using_cli(Can create a bash script using CLI)
    end

    subgraph sub_intermediate_use["Intermediate use of Bianca"]
      can_use_command_line_2(Can use the command line 2)
      can_use_custom_software(Can use custom software)
      can_monitor_jobs(Can monitor jobs)
      can_use_ide(Can use an IDE)
      can_tranfer_files_using_rsync2(Can transfer files using rsync)
    end

    subgraph sub_non_goal["Not in courses"]
      can_tranfer_files_using_gui(Can transfer files using GUI)
      can_create_bash_script_using_gui(Can create a bash script using GUI)
    end


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
    can_tranfer_files_using_rsync --> can_tranfer_files_using_rsync2
    can_schedule_jobs --> can_monitor_jobs
    can_use_interactive_node --> can_use_ide
    can_use_modules --> can_use_ide
    can_tranfer_files_using_rsync2 -.-> |extra| can_tranfer_files_using_gui

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
