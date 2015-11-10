GRTK (Genomic Region Tool Kit) is a package containing various programs for analyzing continuous or binary genomic track data using basic statistics and machine learning. It builds on several popular community packages which perform low-level region manipulation.

# Installation

## Dependencies 

These packages are required:

- GNU parallel
- bedops
- bedtools
- tabix
- various Kent source utilities
- Kyoto Cabinet

### Installing bedops

```
curl http://bedops.googlecode.com/files/bedops_linux_x86_64-v2.2.0.tar.bz2 | sudo tar xvj -C /usr/local
sudo chmod a+x /usr/local/bin/*
```

### Installing the Kent utilities

```
sudo wget -np -R -A "bedToBigBed" -A "bedGraphToBigWig" -A "bigWig*" -A "bigBed*" -N -e robots=off -r -P /usr/local/bin -nd "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/"
sudo curl -o /usr/local/bin/rowsToCols http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/rowsToCols
sudo chmod a+x /usr/local/bin/*
```

### Ubuntu

```
sudo apt-get install -y parallel bedtools tabix kyotocabinet-utils realpath
```

## Building and installing grtk

```
git clone git@bitbucket.org:wrenlab/grtk.git
cd grtk
automake --add-missing
autoreconf
./configure --prefix=/usr/local (or whatever prefix you want)
make install
```

You will be yelled at during the configure step if you're missing dependencies.

# Usage

TODO

# License

AGPL v3
