#!/bin/bash

mkdir -p bin;
tar -xf bwa-0.7.17.tar.bz2; cd bwa-0.7.17/; make; cp ./bwa ../bin; make clean; cd ..; rm -rf bwa-0.7.17;
tar -xf samtools-1.11.tar.bz2; cd samtools-1.11/; ./configure --prefix=$PWD/..; make; make install; make clean; cd ..; rm -rf samtools-1.11/; rm -rf share/;
tar -xf seqtk-v1.3.tar.gz; cd seqtk-1.3/; make; cp seqtk ../bin/; make clean; cd ..; rm -rf seqtk-1.3/;
tar -xf bedops-v2.4.39.tar.gz; cd bedops-2.4.39/; make; make install; cp bin/* ../bin/; make clean; cd ..; rm -rf bedops-2.4.39/;
tar -xf bedtools-2.29.2.tar.gz; cd bedtools2/; make; cp bin/* ../bin/; make clean; cd ..; rm -rf bedtools2/;
tar -xf BBMap_38.87.tar.gz; cp -r bbmap/* bin/; rm -rf bbmap/;
cp -r utils/* bin/
pip3 install --no-cache-dir --prefix $PWD --force-reinstall --ignore-installed git+https://github.com/canfirtina/CrossMap.git
export PATH="${PWD}/bin":$PATH
export PYTHONPATH=$(echo $PWD/lib/python*/site-packages):$PYTHONPATH
echo Installation automatically ran the following command to add the bin directory to your PATH. You may need the copy and paste this command to your "~/.bash_profile" to make sure the environment variables are set automatically everytime you login. Command: $'\n' export PATH="${PWD}/bin":'$PATH'
echo Installation automatically ran the following command to add the bin directory to your PYTHONPATH. You may need the copy and paste this command to your "~/.bash_profile" to make sure the environment variables are set automatically everytime you login. Command: $'\n' export PYTHONPATH='$(echo '$PWD'/lib/python*/site-packages):$PYTHONPATH'
