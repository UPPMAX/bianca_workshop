# Overview

!!! info "Objectives"
    - We'll get an overview of UPPMAX and SNIC/NAISS and how a computer cluster works


## The bigger picture

Here we place UPPMAX within the bigger, national, picture,
starting from the biggest source of money for research in Sweden.

![Vetenskapsrådet logo](./img/vr_logo_128_x_154.png)

[Vetenskapsrådet](https://www.vr.se) ('Science counsel', VR) is biggest funder
of research in Sweden and funds the national HPC infrastructure. 

![NAISS logo](./img/naiss_logo_416_x_68.png)

The [National Academic Infrastructure for Supercomputing in Sweden](https://www.naiss.se/) (NAISS) provides such HPC infrastructure: computing power, storage and data services. Applications for these resources starts at
[this NAISS page](https://www.naiss.se//#application-rounds-for-compute-and-storage-resources). These resources are physically located in multiple places in Sweden,
among other Uppsala.

![UPPMAX logo](./img/uppmax_logo.png)

[Uppsala Multidisciplinary Center for Advanced Computational Science](https://www.uppmax.uu.se/) (**UPPMAX = UppMACS**) provides the HPC infrastructure 
that is physically located in Uppsala

```mermaid
flowchart TD
    HPC_Sweden(HPC in Sweden)
    HPC_others(HPC in other cities)
    HPC_Uppsala(HPC in Uppsala)
    NAISS(NAISS)
    UPPMAX(UPPMAX)
    UU(Uppsala University)
    Users(Users)
    VR(Vetenskapsrådet)

    VR --> |money| HPC_Sweden
    HPC_Sweden -->|done by| NAISS
    NAISS --> |money| HPC_others
    NAISS --> |money| HPC_Uppsala
    HPC_Uppsala -->|done by| UPPMAX
    UU -->|money| HPC_Uppsala
    Users -->|apply for HPC|NAISS
```

## UPPMAX systems

Here we place Bianca between the other UPPMAX systems.

- For computing power, the [UPPMAX clusters](https://www.uppmax.uu.se/resources/systems/) are:
    - Rackham: regular data, general purpose
    - Snowy: regular data, long runs and GPU:s
    - **Bianca: for sensitive data, general use**
    - Miarka: for sensitive data, SciLifeLab-only

A technical summary can be found below.

```mermaid
flowchart TD
    UPPMAX(Which UPPMAX cluster?)
    Bianca
    Rackham
    Miarka
    Snowy
    is_sensitive{Do you use sensitive data?}
    is_scilifelab{Do you work at SciLifeLab?}
    is_long{Do you use long runs and/or GPUs?}

    UPPMAX --> is_sensitive 
    is_sensitive --> |yes|is_scilifelab
    is_scilifelab --> |yes|Miarka
    is_scilifelab --> |no|Bianca
    is_sensitive --> |no|is_long
    is_long --> |no|Rackham
    is_long --> |yes|Snowy
```

- For storage, the [UPPMAX storage systems](https://www.uppmax.uu.se/resources/systems/storage-systems/) are:
    - On-load, active use: Castor of Bianca, Crex for Rackham
    - Off-load, archive: Lutra for Rackham

```mermaid
flowchart TD
    UPPMAX(Which UPPMAX storage system?)
    which_cluster{Which UPPMAX cluster?}
    Castor
    Lutra
    usage_type{Type of use?}

    UPPMAX-->which_cluster
    which_cluster-->|Rackham|usage_type
    which_cluster-->|Bianca|Castor
    usage_type-->|active|Crex
    usage_type-->|archive|Lutra
```

- For cloud servies: use the [UPPMAX cloud](https://www.uppmax.uu.se/resources/systems/the-uppmax-cloud/). It is called 'Dis' (the Swedish word for 'haze') and
the `EAST-1` region of the SNIC science cloud. 

## Bianca is a computer cluster

Bianca is a computer cluster and not a supercomputer.
A computer cluster is a group of computers that can run
many calculations, as requested by multiple people, at the same time.

### Difference between supercomputer and computer cluster

![A supercomputer, from https://en.wikipedia.org/wiki/File:IBM_Blue_Gene_P_supercomputer.jpg](IBM_Blue_Gene_P_supercomputer.jpg)

A supercomputer is a machine that is optimized for doing calculations
quickly. For example, to predict the weather for tomorrow, the calculation
may not take a week. The image above is a supercomputer.

![A computer cluster using a Raspberry Pi ](./img/IMG_5111.jpeg){ width="400" }

![The Rackham computer cluster](./img/uppmax-light2.jpg)

A computer cluster is a machine that is optimized for doing a lot of calculations.
The images above shows a home-made computer cluster above,
and Rackham below, another UPPMAX computer cluster.

Bianca is a computer cluster.

### What is a computer cluster?

A computer cluster is a group of computers that can run
many calculations, as requested by multiple people, at the same time.

Each computer is called a **node**.

There are three types of nodes:

- **login nodes**: nodes where a user enters and interacts with the system
- **calculation nodes**: nodes that do the calculations
- **interactive nodes**: a type of calculation node, where a user can do calculations directly

Each node contains several CPU/GPU cores, RAM and local storage space.

A computer cluster typically is used by multiple people.
To ensure fair use of this shared resource, regular users
are restricted in some ways:
- Users cannot run calculations directly. 
  Instead, users need to request either (1) a calculation to be run,
  or (2) an interactive node
- Users cannot install software directly. 
  Instead, users need to use pre-installed software or learn
  techniques how to run custom software anyway

A user logs in to a login node via the Internet through:

- [http://bianca.uppmax.uu.se/](http://bianca.uppmax.uu.se/)
  for a remote desktop environment
- ThinLinc, for a remote desktop environment
- SSH, for a terminal environment

Note that [http://bianca.uppmax.uu.se/](http://bianca.uppmax.uu.se/) uses
ThinLinc.



![Node](./img/node.png)


  - Here the file management and lighter data analysis can be performed.

![RaspBerry](./img/nodes.png)

![RaspBerry](./img/Bild1.png){ width="400" }

- The **calculation nodes** have to be used for intense computing. 


## Questions

???+ question "At which organization does one start to request access to any Swedish computer cluster?"

    NAISS

???+ question "Which UPPMAX cluster is best for general purpose regular data?"

    Rackham

## Summary

!!! abstract "keypoints"
    - NAISS makes available large-scale high-performance computing resources, storage capacity, and advanced user support, for Swedish research. 
    - UPPMAX runs the local resources placed at Uppsala University
    - A cluster consists of several inter-connected computers that can work individually or together.

## Extra material

### UPPMAX clusters technical summary

|                        |Rackham        |Snowy                     |Bianca                                      |
|------------------------|---------------|--------------------------|--------------------------------------------|
|**Purpose**             |General-purpose|General-purpose           |Sensitive                                   |
|**# Intel CPU Nodes**   |486+144        |228                       |288                                         |
|**# GPU Nodes**         |-              |50, Nvidia T4             |10, 2x Nvidia A100 each                     |
|**Cores per node**      |20/16          |16                        |16/64                                       |
|**Memory per node**     |128 GB         |128 GB                    |128 GB                                      |
|**Fat nodes**           |256 GB & 1 TB  |256, 512 GB & 4 TB        |256 & 512 GB                                |
|**Local disk (scratch)**|2/3 TB         |4 TB                      |4 TB                                        |
|**Login nodes**         |Yes            |No (reached from Rackham) |Yes (2 cores and 15 GB)                     |
|**"Home" storage**      |Domus          |Domus                     |Castor                                      |
|**"Project" Storage**   |Crex, Lutra    |Crex, Lutra               |Castor                                      |

## Detailed overview of the UPPMAX systems

```mermaid

  graph TB

  Node1 -- interactive --> SubGraph2Flow
  Node1 -- sbatch --> SubGraph2Flow
  subgraph "Snowy"
  SubGraph2Flow(calculation nodes) 
        end

        thinlinc -- usr-sensXXX + 2FA + VPN ----> SubGraph1Flow
        terminal -- usr --> Node1
        terminal -- usr-sensXXX + 2FA + VPN ----> SubGraph1Flow
        Node1 -- usr-sensXXX + 2FA + no VPN ----> SubGraph1Flow
        
        subgraph "Bianca"
        SubGraph1Flow(Bianca login) -- usr+passwd --> private(private cluster)
        private -- interactive --> calcB(calculation nodes)
        private -- sbatch --> calcB
        end

        subgraph "Rackham"
        Node1[Login] -- interactive --> Node2[calculation nodes]
        Node1 -- sbatch --> Node2
        end
```
