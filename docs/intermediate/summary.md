# Summary

!!!- info "Learning objectives"

    - Repeat seeing the overview of topics discussed today
    - Share confidence on having learned the topics discussed today

???- question "For teachers"

    Teacher goals are:

    - Learners give feedback on how well topics were taught today
    - Learners give feedback on the course anonymously

    Teaching goals are:

    - Learners have again seen an overview of topics discussed today

    Lesson plan:

    ```mermaid
    gantt
      title Summary
      dateFormat X
      axisFormat %s
      Monologue: 0, 5s
    ```

## Overview of today

Copied from [Introduction](intro.md):

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

    %% Make sure Non-goals is below Intermediate,
    %% using invisible nodes
    can_use_custom_software ~~~ can_tranfer_files_using_gui
```

## Final steps

```mermaid
flowchart TD
  learning_objectives[Fill in learning objectives\nin breakout room]
  evaluation[Fill in evaluation\nin breakout room]
  questions(Questions?)
  done[Enjoy the rest of your day!\nThanks!]
  q_and_a[Go to main room\nQ & A after evaluation]
  learning_objectives --> evaluation
  evaluation --> questions
  questions --> |yes| q_and_a
  questions --> |no| done
  q_and_a --> done
```

## Learning objectives

Most are copied from their respective pages:

- [Introduction](https://uppmax.github.io/bianca_workshop/intermediate/intro/)
    - [ ] See an overview of topics discussed today
    - [ ] See the link to the shared document
    - [ ] See the schedule
- [Transferring files to and from Bianca](https://uppmax.github.io/bianca_workshop/intermediate/transfer/)
	- [ ] Explore the UPPMAX documentation
	- [ ] Understand what the wharf is
	- [ ] Understand what the Transit server allows
	- [ ] Mount the wharf on Transit
	- [ ] Transfer files to/from Bianca using rsync
	- [ ] Transfer files to/from Bianca using FileZilla
- [Slurm jobs](https://uppmax.github.io/bianca_workshop/intermediate/slurm_intermed/)
    - [ ] Understand what Slurm is
    - [ ] Understand some Slurm parameters
    - [ ] Understand what sbatch is
    - [ ] Understand what the job queue is
    - [ ] Start jobs
    - [ ] See job CPU and memory usage
    - [ ] Can start an interactive job
    - [ ] Understand how compute nodes are moved between project clusters
    - [ ] Have used other Slurm tools
- [Software and packages installation](https://uppmax.github.io/bianca_workshop/intermediate/install/)
    - [ ] Understand how to install software yourself
    - Understand how to use Packages and libraries for scripts
        - [ ] ... using Conda
        - [ ] ... using Python packages with pip
        - [ ] ... using R packages
        - [ ] ... using Julia packages
    - [ ] Understand what containers are
    - [ ] Understand what Singularity is
    - [ ] Understand what Docker is
    - [ ] Understand how to build from source
- [IDEs on Bianca](https://uppmax.github.io/bianca_workshop/intermediate/ides/)
    - [ ] Understand what an IDE is
    - [ ] Have heard that RStudio, Jupyter, VSCodium are IDEs
    - [ ] Understand that there are IDEs that can run on Bianca
    - [ ] Have run the voted-for IDE on Bianca
- [NAISS-SENS](https://uppmax.github.io/bianca_workshop/sens_project_short/)
    - [ ] Understand what sensitive personal data is
    - [ ] Understand the difference between pseudonymisation and anonymisation
    - [ ] Know where to apply for project


### Exercise

Share your confidence on having learned the topics discussed today,
by going thought the list on the shared document.

Although it will be messy, between `[ ]`, add a number for confidence:

Grade|Description
-----|------------------------------------
`0`  |I have no idea what this is about
`1`  |I have no confidence I can do this
`2`  |I have low confidence I can do this
`3`  |I have some confidence I can do this
`4`  |I have good confidence I can do this
`5`  |I absolutely can do this!

This may result in a measurement like this:

- `[00101000111201]`: most learners have low confidence
- `[44345454545454]`: most learners have high confidence
