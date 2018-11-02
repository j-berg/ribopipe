#!/bin/bash
#Requires installations of git, wget, and python3
#Usage:
#bash local_install.sh <HISAT2>,<STAR> <yeast>,<mouse>,<human> clone
#or
#bash local_install.sh <HISAT2>,<STAR> <yeast>,<mouse>,<human> version github_tag

echo "Do you want to clone the git repository or download a version of RiboPipe?"
read action

if [$action == 'clone']; then
  git clone https://github.com/j-berg/ribopipe.git
  cd ribopipe/ribopipe/references

if [$action == 'version']; then
  echo "What version tag do you want to fetch?"
  read tag
  wget https://github.com/j-berg/ribopipe/archive/$tag.zip
  unzip ribopipe-${tag:1}.zip
  cd ribopipe-${tag:1}/ribopipe/references

echo "What program do you want to use to align? Options: 'hisat2' or 'star' (case sensitive)"
read program
echo "What model organism do you want to align to? Options: 'yeast' or 'human' or 'mouse' (case sensitive)"
read model

wget https://sourceforge.net/projects/ribopipe/files/${program}_references/${model}_reference_${program}.zip
unzip ${model}_reference_${program}.zip
rm ${model}_reference_${program}.zip
cd ../../
module load python3
python setup.py install --prefix ~/.local
cd ../
