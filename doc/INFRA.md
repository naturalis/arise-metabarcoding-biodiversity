## Introduction

This document records how the infrastructure needs to be configured to 
allow users to run dada2 analyses on Naturalis's HPC infra. The 
[background](#background) section is intended for people setting up a new
server. As this is currently taken care of this information is provided 
only as background context, it does not need to be repeated. The section 
below it describes how to [add users](#user-level-configuration) and must 
be performed for every new user.

## Background

### HPC Hardware

The analyses that are described in this repo require fairly generous
resources in storage, RAM and CPU. Hence, we perform the calculations
on Naturalis's [Metal as a Service](https://maas.netdc.naturalis.io/MAAS/#/machines)
solution for O&O. The MaaS dashboard through which the chosen machine 
is administered is accessible to members of the user group `nmri`,
which includes selected members from ICT as well as the bioinformaticians.

Here we choose machine `netdc-bms-c11h.maas`. This machine has 56 cores,
384GB RAM and 20TB storage space.

### Basic installs

To prepare the server for customization by the user to do the analyses, several 
basic installations need to be performed. Here we setup git (for cloning our
code onto the machine), curl (for fetching things from the web), python3 with
its package manager pip, and [the R runtime](https://linuxize.com/post/how-to-install-r-on-ubuntu-20-04/):

    sudo apt install git libcurl4-openssl-dev python3 python3-pip
    sudo apt install dirmngr gnupg apt-transport-https ca-certificates software-properties-common
    sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
    sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
    sudo apt install r-base

Installation instructions for the RStudio server, which also needs to be installed, are [here](https://www.how2shout.com/linux/install-rstudio-server-open-source-on-ubuntu-20-04-lts/) (link included as reference only).  
To prevent incompatibility issues with the R graphics engine use version 2022.12.0-353 or higher:

    sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
    sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/"
    sudo apt update
    sudo apt-get install gdebi-core
    wget https://download2.rstudio.org/server/bionic/amd64/rstudio-server-2022.12.0-353-amd64.deb
    sudo gdebi rstudio-server-2022.12.0-353-amd64.deb

## User-level configuration

### SSH Access

The newly deployed machine runs Ubuntu 20.04LTS and contains the SSH public
keys of the members of the `nmri` user group. Any others that are going to use
the machine as well will have their public SSH key added to 
`/home/ubuntu/.ssh/authorized_keys` on the server. The machine can then be 
accessed as:

    ssh -i id_rsa ubuntu@145.136.253.38

...where `id_rsa` specifies the location of the private key on the local machine 
that corresponds with the public key previously injected on the server. **Note that 
connection attempts only succeed from behind an EduVPN connection authenticated as 
a Naturalis member.**

### Adding users

_RStudio server_ follows a multi-user model. When accessing the server through
the web browser, the user encounters a login screen. The credentials correspond
with linux users and their passwords. Hence, users must be created:

    sudo su
    adduser --force-badname firstname.lastname
    
...where firstname.lastname follows the format of Naturalis accounts and email
addresses. During this step, the user is prompted to enter a password. This
password corresponds with the login through the server process.

### Accessing _RStudio server_'s port

_RStudio server_ is a webserver process that (by default) listens on
port 8787. This is an unusual port number that is normally blocked by the
EduVPN configuration. This can be circumvented by 
[SSH tunneling](https://www.tecmint.com/create-ssh-tunneling-port-forwarding-in-linux/).
With the following command we enact tunneling as a background process that maps port 
8080 on the user end (one of the usual ports for HTTP traffic) to 8787 on the server:

    ssh -i id_rsa -f -N ubuntu@145.136.253.38 -L 8080:145.136.253.38:8787

When using an EduVPN connection as Naturalis and with this tunneling set up,
it should be possible to access the running webserver at http://localhost:8080/
