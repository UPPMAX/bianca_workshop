# Overview

!!! info "Objectives"
    - We'll get an overview of UPPMAX and SNIC/NAISS and how a computer cluster works


## The bigger picture

![Vetenskapsrådet logo](./img/vr_logo_128_x_154.png)

[Vetenskapsrådet](https://www.vr.se) ('Science counsel', VR) is biggest funder
of research in Sweden and funds HPC infrastructure. 

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

- [UPPMAX clusters](https://www.uppmax.uu.se/resources/systems/) are:
    - Rackham: regular data, general purpose
    - Snowy: regular data, long runs and GPU:s
    - **Bianca: for sensitive data, general use**
    - Miarka: for sensitive data, SciLifeLab-only

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

- [UPPMAX storage](https://www.uppmax.uu.se/resources/systems/storage-systems/)
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

- [UPPMAX cloud](https://www.uppmax.uu.se/resources/systems/the-uppmax-cloud/),
      called 'Dis' (the Swedish word for 'haze'), it is
      the `EAST-1` region of the SNIC science cloud. 

## High Performance Computing — HPC

### What is a cluster?

- A network of computers, each computer working as a **node**.

- From small scale RaspberryPi cluster... 
     
![RaspBerry](./img/IMG_5111.jpeg){ width="400" }

- To supercomputers like Rackham.

![Rackham](./img/uppmax-light2.jpg)

- Each node contains several processor cores and RAM and a local disk called scratch.

![Node](./img/node.png)

- The user logs in to **login nodes** via Internet through ssh or ThinLinc.

  - Here the file management and lighter data analysis can be performed.

![RaspBerry](./img/nodes.png)

![RaspBerry](./img/Bild1.png){ width="400" }

- The **calculation nodes** have to be used for intense computing. 

## Summary about the three "common" UPPMAX clusters

| |Rackham|Snowy|Bianca|
|-------|-----|------|---|
|**Purpose**|General-purpose|General-purpose|Sensitive|
|**#  Nodes (Intel)**|486+144|228+ <br>50 Nvidia T4 GPUs|288 +  <br>10 nodes á 2 <br>NVIDIA A100 GPUs|
|**Cores per node**|20/16|16|16/64|
|**Memory per node**|128 GB|128 GB|128 GB
|**Fat nodes**|256 GB & 1 TB| 256, 512 GB & 4 TB| 256 & 512 GB|
|**Local disk (scratch)**|2/3 TB| 4 TB| 4 TB |
|**Login nodes**|Yes| No (reached from Rackham)|Yes (2 cores and 15 GB)|
|**"Home" storage**|Domus|Domus|Castor|
|**"Project" Storage**|Crex, Lutra|Crex, Lutra|Castor|

## Overview of the UPPMAX systems

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

!!! info "Next session"

    We will try the different forms to log in to Bianca!


!!! abstract "keypoints"
    - NAISS makes available large-scale high-performance computing resources, storage capacity, and advanced user support, for Swedish research. 
    - UPPMAX runs the local resources placed at Uppsala University
    - A cluster consists of several inter-connected computers that can work individually or together.


