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

sudo apt-get update && sudo apt-get install -y \
    build-essential \
    libssl-dev \
    uuid-dev \
    libgpgme11-dev \
    mono-devel

sudo apt-get install golang-go

echo 'export GOPATH=${HOME}/go' >> ~/.bashrc
echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc
source ~/.bashrc

sudo mkdir -p $GOPATH/src/github.com/sylabs
cd $GOPATH/src/github.com/sylabs
sudo git clone https://github.com/sylabs/singularity.git
cd singularity

sudo go get -u -v github.com/golang/dep/cmd/dep

cd $GOPATH/src/github.com/sylabs/singularity
sudo ./mconfig
sudo make -C builddir
sudo make -C builddir install
