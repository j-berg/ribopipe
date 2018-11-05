#!/bin/sh

#For local install of conda on an HPC, etc
cd $HOME

# Start with a clean env
module purge

# Get the latest anaconda3 distro
wget https://repo.anaconda.com/archive/Anaconda3-5.3.0-Linux-x86_64.sh

# Change the permissions -> make it executable
chmod 700 Anaconda3-5.3.0-Linux-x86_64.sh

# Run install
./Anaconda3-5.3.0-Linux-x86_64.sh

#This will launch an interactive script. Follow along, add conda to $PATh via .bashrc
#Exit out of terminal and log in again 
