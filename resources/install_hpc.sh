#For installation on an HPC without root access
git clone https://github.com/j-berg/ribopipe.git
cd ribopipe
module load python3
python setup.py install --prefix ~/.local/
#if an error comes up saying the site-packages folder doesn't exist for python3, create it and re-run install

#Add ribopipe to path
#At the end of the previous output, a line something like:
Installing ribopipe script to /uufs/chpc.utah.edu/common/home/USER/.local/bin
#will appear. Copy this path to your PATH in your bash_profile
echo 'export PATH="/uufs/chpc.utah.edu/common/home/USER/.local/bin:$PATH"' >> ~/.bash_profile
source ~/.bash_profile
#Now you should be able to run ribopipe by simply calling the 'ribopipe' command
