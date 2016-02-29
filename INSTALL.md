# Installation
traitar is available for Linux via the python packaging index. We tested the software under Ubuntu 14.04 LTS (trusty). We didn't test older versions but these might work as well. 
Prior to installation with pip make sure the following packages are installed by running

``sudo apt-get install python-scipy python-matplotlib python-pip python-pandas``
and optionally
``sudo apt-get install python-virtualenv``

These package might not be available for your Linux distribution. Please let us know, if you have troubles installing this software, so that we can assist you with the installation.

Once finished we strongly encourage you to install traitar with pip rather than by manually cloning the repository. Install locally by

``pip install traitar-<version>.tar.gz --user``

You need to add the following line to your ~/.bashrc to adjust the PATH variable so that the bash finds the executables needed for running traitar. 

``PATH=$PATH:~/.local/bin/``

You may need to run ``source ~/.bashrc`` once in your current session.

You can also install globally with

``sudo pip install traitar-<version>.tar.gz``
### Creating a virtual environment
You may want to use virtualenv to create a clean environment for traitar i.e. run

```
virtualenv <environment_name> --system-site-packages
source <environment_path>/bin/activate

pip install -U traitar-<version>.tar.gz
PATH=$PATH:<environment_path>/bin/
source ~/.bashrc
```

### Additional requirements
traitar further needs prodigal and hmmsearch available on the command line. For parallel execution it further requires GNU parallel.
All three are available as preconfigured package for many Linux installation e.g. for Debian / Ubuntu. For Ubuntu 14.04 it is available via the trusty-backports.

```
#include trusty backports source in the sources.list
"sudo deb http://archive.ubuntu.com/ubuntu trusty-backports main restricted universe multiverse ">> /etc/apt/sources.list

#update sources
apt-get install update

#install packages
sudo apt-get install parallel prodigal hmmer
```

``traitar -h`` will provide help regarding the different options of traitar.
traitar requires the Pfam 27.0 HMM models. These are not distributed with this package but need to be downloaded

``traitar pfam <path to download folder>``

will download and extract the Pfam models to the specified folder. The download can take a while depending on you internet connection. The Pfam models are extracted automatically to the same folder. This can slow down some machines. 
You may also download and extract the Pfam models manually from ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam27.0/Pfam-A.hmm.gz and run 

``traitar pfam --local <path to Pfam folder>``

to let traitar know where.
