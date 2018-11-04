#!/bin/bash
#RiboPipe
#An assembly and analysis pipeline for sequencing data
#alias: ripopipe

#Copyright (C) 2018  Jordan A. Berg
#jordan <dot> berg <at> biochem <dot> utah <dot> edu

#This program is free software: you can redistribute it and/or modify it under
#the terms of the GNU General Public License as published by the Free Software
#Foundation, either version 3 of the License, or (at your option) any later
#version.

#This program is distributed in the hope that it will be useful, but WITHOUT ANY
#WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
#PARTICULAR PURPOSE. See the GNU General Public License for more details.

#You should have received a copy of the GNU General Public License along with
#this program.  If not, see <https://www.gnu.org/licenses/>.

#Requires installations of git, wget, and python3
#Usage:
#bash local_install.sh <HISAT2>,<STAR> <yeast>,<mouse>,<human> clone
#or
#bash local_install.sh <HISAT2>,<STAR> <yeast>,<mouse>,<human> version github_tag

echo "Do you want to <clone> the git repository or download a <version> of RiboPipe?"
read action

if [$action == 'clone']; then
  git clone https://github.com/j-berg/ribopipe.git
  cd ribopipe/ribopipe/references

if [$action == 'version']; then
  echo "What version tag do you want to fetch?"
  read tag
  wget https://github.com/j-berg/ribopipe/archive/$tag.zip
  unzip ribopipe-${tag:1}.zip
  rm ribopipe-${tag:1}.zip
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

#Install htseq
echo "Do you want to <clone> the git repository or download <version> 0.11.0 of HTSeq?"
read action2

if [$action == 'clone']; then
  git clone https://github.com/simon-anders/htseq.git
  cd htseq
  python setup.py install --prefix ~/.local
  cd ../

if [$action == 'version']; then
  echo "Fetching HTSeq v0.11.0"
  wget https://github.com/simon-anders/htseq/archive/release_0.11.0.zip
  unzip htseq-release_0.11.0.zip
  rm htseq-release_0.11.0.zip
  cd htseq-release_0.11.0
  python setup.py install --prefix ~/.local
  cd ../
