# HPC Hardware

The analyses that are described in this repo require fairly generous
resources in storage, RAM and CPU. Hence, we perform the calculations
on Naturalis's [Metal as a Service](https://maas.netdc.naturalis.io/MAAS/#/machines)
solution for O&O. The MaaS dashboard through which the chosen machine 
is administered is accessible to members of the user group `nmri`,
which includes selected members from ICT as well as the bioinformaticians.

Here we choose machine `netdc-bms-c11h.maas`. This machine has 56 cores,
384GB RAM and 20TB storage space.

# SSH Access

The newly deployed machine runs Ubuntu 20.04LTS and contains the SSH public
keys of the members of the `nmri` user group. If others are going to use
the machine, their public SSH key needs to be added to `/home/ubuntu/.ssh/authorized_keys`.
The machine can then be accessed as:

    ssh -i id_rsa ubuntu@145.136.253.38

...where `id_rsa` specifies the location of the private key that corresponds
with the public key previously injected in the machine. **Note that connection
attempts only succeed from behind an EduVPN connection authenticated as a Naturalis
member.**

# Accessing _RStudio server_'s port

_RStudio server_ is a webserver process that (by default) listens on
port 8787. This is an unusual port number that is currently blocked by the
EduVPN configuration. This can circumvented by SSH tunneling. With the following
command we enact tunneling as a background process that maps port 8080 on the
user end (one of the usual ports for HTTP traffic) to 8787 on the server:

    ssh -i id_rsa -f -N ubuntu@145.136.253.38 -L 8080:145.136.253.38:8787

When using an EduVPN connection as Naturalis and with the tunneling set up,
it should be possible to access a running webserver at http://localhost:8080/

# Adding users

_RStudio server_ follows a multi-user model. When accessing the server through
the web browser, the user encounters a login screen. The credentials correspond
with linux users and their passwords. Hence, users must be created:

    sudo su
    adduser --force-badname firstname.lastname
    
...where firstname.lastname follows the format of Naturalis accounts and email
addresses. During this step, the user is prompted to enter a password. This
password corresponds with the login through the server process.
