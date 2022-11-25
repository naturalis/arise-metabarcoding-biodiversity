# Introduction

This section describes how to clone a project from github
into the RStudio webserver.

## Generating SSH keys

To be able to clone from, and push to, github, RStudio needs
to use an SSH keypair. The simplest way to do this is to go to
'Tools > Global options' in the server interface, go to the 
'git/svn' tab, and create a new keypair.

## Copy pubkey to github

Using the RStudio terminal, copy the pub key from
~/.ssh/id_rsa.pub and add it to github: https://github.com/settings/keys

## Cloning the repo as a project

The code from the git repo can now be cloned as a project,
using 'File > New project', specifying the special URL
`git@github.com:naturalis/arise-metabarcoding-biodiversity.git`
