# CryptoMMU
The artifact of CryptoMMU in [Zenodo](https://doi.org/10.5281/zenodo.8287142) provides ready-to-run virtual machine (VM),
however, that repo requires users to download huge data at once with low portability

This repo describes how to set up CryptoMMU in Ubuntu 22.04 modularly without interferencing other users and library tools:
+ Prerequisite: Automatically install tools
+ Prerequisite: Manually install remaining tools
+ Build & Run


## Prerequisite: Automatically install tools
Assuming that our Ubuntu machine is in vanilla state

1. Update apt list and reboot

        $ sudo apt update
        $ sudo apt upgrade
        $ reboot now 

2. Install other essential tools (some might be already installed)

        $ sudo apt install build-essential python2 python2-dev vim git m4 zlib1g zlib1g-dev libtool libltdl-dev libfabric-dev cmake kmod uuid-dev bison doxygen perl valgrind htop iptables libhdf5-dev patch gettext libtool-bin ca-certificates gnupg


## Prerequisite: Manually install remaining tools
Let's assume we start from $HOME. If the wget-server provider does not work, you can download tools below through this [link](https://drive.google.com/drive/folders/1HdZvY4MGNwxd-EdZA5oMysj2W7hX6SSH?usp=drive_link)

1. Download repo and disk image for full system simulation

        $ git clone https://github.com/harrylee365/cryptommu.git
        $ cd $HOME/cryptommu
        $ wget http://www.sfu.ca/~zhenman/files/software/disk-binary.tar.gz
        $ tar zxvf disk-binary.tar.gz
        $ cp -r disks $HOME/cryptommu/CryptoMMU && cp -r binaries $HOME/cryptommu/CryptoMMU
        $ cp -r disks $HOME/cryptommu/CryptoMMU_ReadAcc && cp -r binaries $HOME/cryptommu/CryptoMMU_ReadAcc
        $ cp -r disks $HOME/cryptommu/Full_IOMMU && cp -r binaries $HOME/cryptommu/Full_IOMMU
        $ cp -r disks $HOME/cryptommu/Border_Control && cp -r binaries $HOME/cryptommu/Border_Control
        $ cp -r disks $HOME/cryptommu/ATS-only-IOMMU && cp -r binaries $HOME/cryptommu/ATS-only-IOMMU

 2. Download tools not available through usual `apt install`

        $ mkdir $HOME/toolset 
        $ cd $HOME/toolset
        $ wget  https://sourceforge.net/projects/scons/files/scons/1.3.1/scons-1.3.1.tar.gz && tar zxvf scons-1.3.1.tar.gz
        $ wget https://sourceforge.net/projects/swig/files/swig/swig-2.0.9/swig-2.0.9.tar.gz && tar zxvf swig-2.0.9.tar.gz
        $ wget https://download.mono-project.com/sources/mono/mono-6.12.0.199.tar.xz && tar xvf mono-6.12.0.199.tar.xz
        $ wget https://github.com/google/protobuf/releases/download/v2.5.0/protobuf-2.5.0.tar.gz && tar zxvf protobuf-2.5.0.tar.gz

3. In the following procedure, we will install all these library tools in the isolated directory `volatile` without annoying other users on the server.

        $ mkdir ~/volatile

### Install Scons

        $ cd scons-1.3.1
        $ python2 setup.py install --prefix=$HOME/volatile

### Install Mono 
Assuming that ca-certificates and gnupg are installed

Method 1: install using dowanload file
1. Automake

        $ cd $HOME/toolset/mono-6.12.0.199
        $ ./configure --prefix=$HOME/volatile

2. Install

        $ make 
        $ make install

Method 2: install using apt (not recommended)

1. Add Mono repository to our system using apt

        $ sudo gpg --homedir /tmp --no-default-keyring --keyring /usr/share/keyrings/mono-official-archive-keyring.gpg --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF
        $ echo "deb [signed-by=/usr/share/keyrings/mono-official-archive-keyring.gpg] https://download.mono-project.com/repo/ubuntu stable-focal main" | sudo tee /etc/apt/sources.list.d/mono-official-stable.list
        $ sudo apt update

2. Install

        $ sudo apt install mono-devel

### Install Swig
The procedure is similar to Mono

        $ cd $HOME/toolset/swig-2.0.9
        $ ./configure --prefix=$HOME/volatile
        $ make 
        $ make install

### Install protobuf
Installing protobuf is similar to Mono and Swig. 
But let's skip this to avoid linking protobuf in gem5 as it still has some version issues

### Install older gcc and g++ (4.8)
Since CryptoMMU project is based on old gem5, it needs gcc & g++ higher than v4.7. 
However, Ubuntu 22.04 does not provide it, so we need to add the old repo to our apt system. 
Please note that this is recommended to be done as the last procedure, cause building other tools need new gxx

1. Add old Ubuntu repo to our system using apt (you can exclude kr to use US server)

        $ sudo gpg --homedir /tmp --no-default-keyring --keyring /usr/share/keyrings/kr-archive-ubuntu-keyring.gpg --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3B4FE6ACC0B21F32 40976EAF437D05B5
        $ echo "deb [signed-by=/usr/share/keyrings/kr-archive-ubuntu-keyring.gpg] http://kr.archive.ubuntu.com/ trusty universe main" | sudo tee /etc/apt/sources.list.d/kr-archive-ubuntu.list
        $ sudo apt update

2. Install 

        $ sudo apt install gcc-4.8
        $ sudo apt install g++-4.8

3. This repository requires gcc/g++>=4.7 and python>=2.5. However, directly modifying versions might annoy other users on the server. Hence, we can make symbolic links and temporarily use them for each terminal open. 

        $ ln -s /usr/bin/python2.7 ~/volatile/bin/python
        $ ln -s /usr/bin/gcc-4.8 ~/volatile/bin/gcc
        $ ln -s /usr/bin/g++-4.8 ~/volatile/bin/g++

4. Clean up added server. First, commentize the list in **kr-archive-ubuntu.list** and **mono-official-stable.list**. Then, clean up cached list as follows:

        $ sudo rm -rf /var/lib/apt/lists/*
        $ sudo apt update


## Build & Run

1. Set enviorment variables to use gcc, g++, scons, mono, and swig installed by us in *current terminal*. Please note that the directory order should be followed in a strict manner. **The variables should be set everytime we open the terminal**

        $ export PATH=~/volatile/bin:$PATH
        $ export LD_LIBRARY_PATH=~/volatile/lib:$LD_LIBRARY_PATH

2. The repository provides five models, as explained in the paper. Let's assume that we are interested in CryptoMMU model. **The variable should be set everytime we run differnt model**

        $ export M5_PATH=$HOME/cryptommu/CryptoMMU
        $ cd $M5_PATH

3. Build 

        $ ./build.gem5.sh

4. Run a benchmark, Denoise. Option `-c` is the number L2$ banks, whereas `-a` is the number of accelerators

        $ ./run_bench.sh Denoise -c X -a Y


## Contributors
+ Faiz Alam                 falam3@ncsu.edu
+ Hyokeun Lee               ~~hlee48@ncsu.edu~~ hyokeunlee@ajou.ac.kr
+ Abhishek Bhattacharjee    abhishek@cs.yale.edu
+ Amro Awad                 ajawad@ncsu.edu


## Citation
```
@inproceedings{cryptommu-micro2023, 
author = {Alam, Faiz and Lee, Hyokeun and Bhattacharjee, Abhishek and Awad, Amro},
title = {CryptoMMU: Enabling Scalable and Secure Access Control of Third-Party Accelerators},
booktitle = {IEEE/ACM International Symposium on Microarchitecture (MICRO)},
year = {2023}
}
```

