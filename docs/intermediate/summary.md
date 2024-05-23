# Summary and evaluation

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

## Exercises

Goal of these exercises are:

- You reflect on today
- We learn from your feedback

This is an anonymous evaluation
and will ideally be published in raw form at 
[https://github.com/UPPMAX/bianca_workshop/tree/main/evaluations/20240524](https://github.com/UPPMAX/bianca_workshop/tree/main/evaluations/20240524).
To do so, please do not share sensitive data here!

This evaluation uses a shared document, 
because its advantages outweigh the disadvantages.

The drawbacks of using a shared document are:

- one needs to be careful when multiple people are editing at the same time. This may result in minor data loss possible while editing

The advantages that outweigh these are: 

- it is already there, hence there is no new link needed
- you already have worked with it, hence there will be no new technical problems
- it guarantees integrity: a literal copy-paste of the data will perfectly preserve the data
- there is no owner of this data: all learners and teachers can access, verify and upload the data

You will work on these exercises in a shared Breakout room,
as having teachers present gives a bias in your feedback.

### Exercise 1: what did I learn today?

The goal of this exercise is for you to reflect on what you learned today:

- how many things?
- how confident are you now?

To do so, the learning objectives are collected in the shared document.

There, you can share your confidence on having learned the topics discussed today,
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

Thanks for your time in helping us understand where we can improve
our teaching of these topics! 

### Exercise 2: evaluation

The goal of this exercise is for you to reflect on the course.
How did you enjoy the course and how would you improve it?

To do so, the evaluation questions are at the bottom of the share document.

Thanks for your time  in helping us improve our course as a whole!
