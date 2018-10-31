#For installation on an HPC without root access
git clone https://github.com/j-berg/ribopipe.git
cd ribopipe
module load python3
python setup.py install --prefix ~/.local/
#if an error comes up saying the site-packages folder doesn't exist for python3, create it and re-run install
