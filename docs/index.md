---
template: home.html
---

<center>

<br/><br/>

<img src="assets/UU_logo_color.svg" alt="drawing" width="200"/>

<br/><br/>


# Welcome to the Bianca workshop!
    
</center>

!!! info "Introduction to Bianca: Handling Sensitive Research Data"
    
    - Are you *starting to* work with your sensitive data in your research? 
    - If yes, welcome to [a full day introduction to handling sensitive data on the UPPMAX cluster, Bianca](intro.md). 

    - Do you want to *deepen your existing knowledge* to work 
      with your sensitive data in new, more flexible ways? 
    - If yes, welcome to [a full day intermediate-level workshop/hackathon to handling sensitive data on the UPPMAX cluster, Bianca](intermediate/intro.md). 

    - You will learn about 
    
        - NAISS-SENS, i.e. the legal and administrative aspects
        - how to login to Bianca
        - transferring files
        - the job scheduler system
        - using pre-installed software
        - installing your own software

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
      can_use_ide(Can use an IDE)
    end
    style sub_basic_use fill:#ff0,color:#000,stroke:#ffa

    subgraph sub_intermediate_use["Intermediate use of Bianca"]
      can_use_custom_software(Can use custom software)
      can_monitor_jobs(Can monitor jobs)
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

<center>
<br>
    
<br/><br/>

[Practicalities](practicalities.md){ .md-button .md-button--primary }
[Bianca intro ](intro.md){ .md-button .md-button--primary }
[Intermediate workshop](intermediate/intro.md){ .md-button .md-button--primary }

<br/><br/>


</center>
