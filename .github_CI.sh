#!/bin/bash
# This script contains all the commands executed by the CI of github


# Clone Hoppet
git clone --branch 2024-01-n3lo-splittings-functions https://github.com/hoppet-code/hoppet.git
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
wget https://lhapdf.hepforge.org/downloads/LHAPDF-6.5.4.tar.gz
tar -xzvf LHAPDF-6.5.4.tar.gz
cd LHAPDF-6.5.4
./configure
make -j
sudo make install
