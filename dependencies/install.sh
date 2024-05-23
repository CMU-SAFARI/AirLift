#!/bin/bash

#If you have already installed the environment and packages, comment out the following line and just execute the conda activate command below.
conda create -n airlift-env samtools=1.20 seqtk=1.4 bedops=2.4.41 bedtools=2.31.1 bbmap=39.06 bwa=0.7.18 crossmap=0.7.0 fastremap-bio=1.0.0 python=3.10.14

conda activate airlift-env

# mkdir -p bin;
# tar -xf bwa-0.7.18.tar.gz; cd bwa-0.7.18/; make; cp ./bwa ../bin; make clean; cd ..; rm -rf bwa-0.7.18;
# tar -xf samtools-1.20.tar.bz2; cd samtools-1.20/; ./configure --without-curses --prefix=$PWD/..; make; make install; make clean; cd ..; rm -rf samtools-1.20/; rm -rf share/;
# tar -xf seqtk-v1.3.tar.gz; cd seqtk-1.3/; make; cp seqtk ../bin/; make clean; cd ..; rm -rf seqtk-1.3/;
# tar -xf bedops-v2.4.39.tar.gz; cd bedops-2.4.39/; make; make install; cp bin/* ../bin/; make clean; cd ..; rm -rf bedops-2.4.39/;
# tar -xf bedtools-2.29.2.tar.gz; cd bedtools2/; make; cp bin/* ../bin/; make clean; cd ..; rm -rf bedtools2/;
# tar -xf BBMap_38.87.tar.gz; cp -r bbmap/* bin/; rm -rf bbmap/;
# cp -r utils/* bin/
# pip3 install --no-cache-dir --prefix $PWD --force-reinstall --ignore-installed git+https://github.com/canfirtina/CrossMap.git
# export PATH="${PWD}/bin":$PATH
# export PYTHONPATH=$(echo $PWD/lib/python*/site-packages):$PYTHONPATH
# echo Installation automatically ran the following command to add the bin directory to your PATH. You may need the copy and paste this command to your "~/.bash_profile" to make sure the environment variables are set automatically everytime you login. Command: $'\n' export PATH="${PWD}/bin":'$PATH'
# echo Installation automatically ran the following command to add the bin directory to your PYTHONPATH. You may need the copy and paste this command to your "~/.bash_profile" to make sure the environment variables are set automatically everytime you login. Command: $'\n' export PYTHONPATH='$(echo '$PWD'/lib/python*/site-packages):$PYTHONPATH'

