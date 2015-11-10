#!/usr/bin/env bash

# Put these lines in .bashrc
PREFIX=$HOME/.local
export PATH=$PREFIX:$PREFIX/lib:$PREFIX/bin:$PATH

# Install prerequisites
mkdir -p $PREFIX/bin
sudo apt-get install -y vim
sudo apt-get install -y git
sudo apt-get install -y dos2unix
sudo apt-get install -y parallel
sudo apt-get install -y build-essential g++
sudo apt-get install -y gdebi-core
sudo apt-get install -y curl

# Misc libraries required by prerequisites
sudo apt-get install -y zlib1g-dev
sudo apt-get install -y libreadline-dev
sudo apt-get install -y libpcre3 libpcre3-dev
sudo apt-get install -y liblzma-dev
sudo apt-get install -y libbz2-dev
sudo apt-get install -y libatlas-base-dev
sudo apt-get install -y gfortran

# Required by R packages
sudo apt-get install -y libcurl4-openssl-dev
sudo apt-get install -y libssl-dev
sudo apt-get install -y libxml2-dev

# Install python-related stuff
sudo apt-get install -y python2.7
sudo apt-get install -y python2.7-dev
sudo apt-get install -y python-dev
sudo apt-get install -y python-pip
sudo apt-get install -y python-scipy
sudo apt-get install -y python-numpy
sudo apt-get install -y python-rpy2
sudo apt-get install -y python-beautifulsoup
sudo apt-get install -y python-bs4
sudo apt-get install -y python-mako
sudo apt-get install -y python-simplejson
sudo apt-get install -y python-cherrypy3
sudo apt-get install -y python-celery
sudo apt-get install -y python-redis
sudo apt-get install -y python-singledispatch

# Additional prerequisites
sudo apt-get -y upgrade gcc 
sudo pip install -U cython
sudo apt-get -y install redis-server
sudo pip install celery
sudo pip install flower # Tool to monitor Celery jobs
sudo pip install -U Celery # Helps to solve issue #12, 'module' object has no attribute 'celeryconfiguration'
sudo apt-get -y install kyotocabinet-utils 
sudo apt-get install realpath
sudo apt-get install bedtools

# install R version 3.2.1 for Ubuntu 14.04
wget https://cran.rstudio.com/bin/linux/ubuntu/trusty/r-base-core_3.2.1-4trusty0_amd64.deb
sudo gdebi -n r-base-core_3.2.1-4trusty0_amd64.deb 
wget https://cran.rstudio.com/bin/linux/ubuntu/trusty/r-recommended_3.2.1-4trusty0_all.deb
sudo gdebi -n r-recommended_3.2.1-4trusty0_all.deb 
wget https://cran.rstudio.com/bin/linux/ubuntu/trusty/r-doc-html_3.2.1-4trusty0_all.deb
sudo gdebi -n r-doc-html_3.2.1-4trusty0_all.deb 
wget https://cran.rstudio.com/bin/linux/ubuntu/trusty/r-base_3.2.1-4trusty0_all.deb
sudo gdebi -n r-base_3.2.1-4trusty0_all.deb 

# Manual download and installation of required binaries
mkdir downloads
cd downloads
curl -L https://github.com/bedops/bedops/releases/download/v2.4.14/bedops_linux_x86_64-v2.4.14.tar.bz2 | tar -xjvf - -C $PREFIX

sudo wget -np -R -A "bedToBigBed" -A "bedGraphToBigWig" -A "bigWig*" -A "bigBed*" -N -e robots=off -r -P $PREFIX/bin -nd "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/"
sudo wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/rowsToCols
sudo mv rowsToCols $PREFIX/bin
# Install tabix v0.2.5
git clone https://github.com/mdozmorov/tabix.git
cd tabix/
make
sudo mv tabix bgzip $PREFIX/bin/
cd ../..
sudo rm -r downloads
sudo chmod a+x $PREFIX/bin/*

# Install R packages
sudo Rscript -e 'install.packages(c("Hmisc", "RColorBrewer", "gplots", "xml2", "curl", "httr", "RCurl", "rversions", "git2r", "devtools", "shiny", "shinyBS", "DT", "dendextendRcpp", "colorRamps", "dplyr", "scales"), repos="http://cran.revolutionanalytics.com")'
sudo Rscript -e 'devtools::install_github("mdozmorov/d3heatmap")'
sudo Rscript -e 'devtools::install_github("mdozmorov/MDmisc")'
# Check if all installed correctly, in R
# require("Hmisc"); require("RColorBrewer"); require("gplots"); require("xml2"); require("curl"); require("httr"); require("RCurl"); require("rversions"); require("git2r"); require("devtools"); require("shiny"); require("shinyBS"); require("DT"); require("dendextendRcpp"); require("colorRamps"); require("dplyr"); require("scales"); require("d3heatmap")

# Installing Genomic Region Tool Kit
cd grtk
sudo python setup.py install
cd ..
# Finally, installing GenomeRunner itself. Keep uncommented only one type of installation
# Developmental mode. Changes made in github-cloned folder are immediately active
sudo python setup.py install develop -d $PREFIX/lib/python2.7/dist-packages/
# Standard mode, default. Changes made in github-cloned folder require reinstallation to be active
#sudo python setup.py install
