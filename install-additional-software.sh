#! /usr/bin/env sh
git clone https://github.com/lh3/bwa
cd bwa
make
cd ..

git clone -b master https://github.com/samtools/samtools
cd samtools
make
cd ..

git clone https://github.com/docopt/docopt.cpp.git src/docopt
