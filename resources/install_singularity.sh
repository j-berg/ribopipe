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
