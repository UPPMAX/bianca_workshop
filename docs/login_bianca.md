# Log in to Bianca

!!! info "Objectives" 

    - First step in understanding why Bianca login is the way it is
    - Log in via ThinLinc
    - Log in via SSH
    - First step in understanding what a login node is 

## Exercises

 1. Discuss: what is the purpose of Bianca? What kind of consequences will this have for its design?
 2. Login via ThinLinc
 3. Login via SSH
 4. Start an interactive session

![The voyage from outside the university network to a cluster login node](971_the_voyage_from_outside_the_university_network_to_a_cluster_login_node.png)

> The voyage from outside the university network to a cluster login node

## Bianca's design

Bianca was designed to, among other:

 * Protect the sensitive data
   * Accidental data leaks should be difficult
   * Law: if data is leaked, the person doing so should be possibly identified
 * Emulate a standard HPC cluster environment
   * Use the hardware as efficient as possible, 
     by using a queuing system.
     See [UPPMAX usage for the current usage of Bianca](https://status.uppmax.uu.se/usage/)
   * Distributes shared resources (CPU, memory) in a fair way,
     by using a queuing system
   * make correct data management as easy as possible
    
### Bianca and the Internet

![Bianca](./img/biancaorganisation-01.png)

> The relation between Bianca and the Internet

Bianca and the Internet have this relation:

 * Bianca has no internet [1], to prevent accidental data leaks. 
 * Bianca is only accessible from within SUNET (i.e. from university networks),
   to protect the sensitive data better.

Login depends on where you are and what you need:

Where you are|What you need   |What to do
-------------|----------------|------------------------------
Inside SUNET |A terminal      |`ssh` into Bianca
Inside SUNET |A remote desktop|Login at [https://bianca.uppmax.uu.se](https://bianca.uppmax.uu.se)
Outside SUNET|A terminal      |`ssh` into a VPN/Rackham first
Outside SUNET|A remote desktop|Login at [https://bianca.uppmax.uu.se](https://bianca.uppmax.uu.se) via a VPN

The be able to use VPN:

 * For Uppsala Univerity: [go to this page](https://mp.uu.se/en/web/info/stod/it-telefoni/anvandarguider/network/vpn-service)
 * For other Swedish universities, search their websites to get a VPN setup

Data can be transferred to/from the `wharf`, 
which is a special folder that is visible from the Internet.

## Logging in

You can log in either through ThinLinc or via SSH, `ssh`:

 * ThinLinc: provides a remote desktop, needed for using graphical tools
 * SSH: provides only the command line
    - ``ssh`` from home terminal
    - ``ssh`` from a session on Rackham 

### Log in to the Bianca remote desktop environment

!!! warning

    You need to be within SUNET or use a VPN

Bianca offers a remote desktop environment (which uses ThinLinc to establish
the connection). Here is how to login:

 1. In your web browser, go to [https://bianca.uppmax.uu.se](https://bianca.uppmax.uu.se)

 2. Fill in the first dialog. Do use the `UPPMAX` [2-factor authentication](https://www.uppmax.uu.se/support/user-guides/setting-up-two-factor-authentication/) (i.e. not SUPR!)

![Bianca login, first dialog](./img/bianca_gui_login_1st.png)

 3. Fill in the second dialog, using your regular password (i.e. no need for two-factor authentication)

![Bianca login, second dialog](./img/bianca_gui_login_2nd.png)

> The second Bianca remote desktop login dialog. 
> Note that it uses ThinLinc to establish this connection

 4. Enjoy! You are in!

![The Bianca remote desktop](./img/bianca_remote_desktop.png)

> The Bianca remote desktop

### Log in to the Bianca command-line environment

You can use your favorite terminal to login (see <https://uppmax.github.io/uppmax_intro/login2.html#terminals> for an overview of many)
to the Bianca command-line environment.

  1. From a terminal, within SUNET, use `ssh` to log in:

```bash
ssh [user]-[project name]@bianca.uppmax.uu.se
```

For example:

```bash
ssh richel-sens2023531@bianca.uppmax.uu.se
```

 2. Type your UPPMAX password, 
    directly followed by the UPPMAX 2-factor authentication number,
    for example `verysecret678123`, then press enter

 3. Type your UPPMAX password,
    for example `verysecret`

 4. Enjoy! You are in!

## Login node

When you are logged in, you are on a login node.
There are two types of nodes:

Type        |Purpose
------------|--------------------------
Login node  |Start jobs for worker nodes, do easy things
Worker node |Do hard calculations, either from scripts of an interactive session

Bianca contains hundreds of nodes, each of which is isolated from each other and the Internet.

As Bianca is a shared resources, there are rules to use it together in fair way:

 * The login node is only for easy things, such as moving files,
   starting jobs or starting an interactive session
 * The worker nodes are for harder things, such as
   running a script or running an interactive session.

To start an interactive session [2], type:


```bash
interactive -A [project name] -p core -n 2 -t 8:0:0
```

For example:

```bash
interactive -A sens2023531 -p core -n 2 -t 8:0:0
```

## Conclusions

 * Bianca makes it hard to leak data
 * Login differs from where you are and what you need
 * Only do light things on login nodes

## Footnotes

 * [1] 'no internet' meaning 'no direct way to download or upload data from/to
   the internet
 * [2] In this case, 8 hour long, with 2 cores

## Videos

 * Login from inside SUNET: [YouTube](https://youtu.be/upBozh2BI5c), [download (.ogv)](https://richelbilderbeek.nl/login_bianca_inside_sunet.ogv)
 * Login to Bianca from outside SUNET, no VPN: [YouTube](https://youtu.be/W-PMTyNcbYI), [download (.mp4)](https://richelbilderbeek.nl/login_bianca_outside_sunet.mp4)
 * Login to Bianca remote desktop from outside SUNET, with VPN: [YouTube](https://youtu.be/AIJKbJeu0MI), [download (.ogv)](https://richelbilderbeek.nl/login_bianca_outside_sunet_vpn.ogv)

## Links

 * [The Bianca remote desktop login](https://bianca.uppmax.uu.se)
 * [How get a VPN for UU](https://mp.uu.se/en/web/info/stod/it-telefoni/anvandarguider/network/vpn-service)
