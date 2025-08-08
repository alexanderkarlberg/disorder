#!/bin/bash
# This script contains all the commands executed by the CI of github

# Get sem
sudo apt install parallel

# Clone Hoppet
git clone https://github.com/hoppet-code/hoppet.git
cd hoppet
# Compile, check and install Hoppet
./configure
make -j
make check
sudo make install
cd ..

# Do the same for fastjet
wget https://fastjet.fr/repo/fastjet-3.4.2.tar.gz
tar -xzvf fastjet-3.4.2.tar.gz
cd fastjet-3.4.2
./configure
make -j
make check
sudo make install
cd ..

# And LHAPDF
wget https://lhapdf.hepforge.org/downloads/LHAPDF-6.5.5.tar.gz
tar -xzvf LHAPDF-6.5.5.tar.gz
cd LHAPDF-6.5.5
./configure --disable-python
make -j
sudo make install
wget http://lhapdfsets.web.cern.ch/lhapdfsets/current/MSHT20an3lo_as118.tar.gz
tar -xzvf MSHT20an3lo_as118.tar.gz
sudo cp -r MSHT20an3lo_as118 /usr/local/share/LHAPDF/.

# Set dynamic library path
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
export LD_RUN_PATH=$LD_RUN_PATH:/usr/local/lib
sudo ldconfig

ldconfig -v 2>/dev/null | grep -v ^$'\t'
