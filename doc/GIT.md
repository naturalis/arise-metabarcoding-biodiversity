# Introduction

This section describes how to clone a project from github
into the RStudio webserver.

## Generating SSH keys

To be able to clone from, and push to, github, RStudio needs
to use keypair. The simplest way to do this is to go to
'Tools > Global options', go to the 'git/svn' tab, and create
a new keypair.

## Copy pubkey to github

Using the RStudio terminal, copy the pub key from
~/.ssh/id_rsa.pub and add it to github: https://github.com/settings/keys
