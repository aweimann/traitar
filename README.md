
# traitAR - the microbial trait analyzer
traitAR is a software for characterizing microbial samples from nucleotide or protein sequences. It can accurately phenotype 67 diverse traits related to growth, shape, carbon source utilized etc.

# Installation
traitAR is available for Linux via the python packaging index. We strongly encourage to install it with pip rather than manually by cloning the repository. Install locally by

``pip install -i https://testpypi.python.org/pypi  --user traitar``

You might need further packages that are not available via pip e.g. freetype, which can be installed via your software managment system. Please let us know, if you have troubles installing this software, so that we can assist you with the installation. 
You need to add the following line to your ~/.bashrc to adjust the PATH variable so that the bash finds the executables provided by packages installed with pip. 

``PATH=$PATH:~/.local/bin/``

You may need to run ``source ~/.bashrc`` once in your current session.

You can also install globally with

``sudo pip install -i https://testpypi.python.org/pypi traitar``

traitAR further needs Prodigal and HMMer available on the command line. For parallel execution it further requires GNU parallel.
All three are available as preconfigured package for the major Linux installation e.g. for Debian / Ubuntu install by

``sudo apt-get install parallel prodigal HMMer``

``traitar -h`` will provide help regarding the differnt options of traitAR.
traitAR requires the Pfam 27.0 HMM models. These are not distributed with this package but need to be downloaded

``traitar config <path to folder>``

will download and extract the Pfam models to the specified folder. The download can take a while depending on you internet connection. The Pfam models are extracted automatically to the same folder. This can slow down some machines. 
You may also download and extract the Pfam models manually from ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam27.0/Pfam-A.hmm.gz and run 

``traitar config --local <path to Pfam folder>``

to let traitAR know where.
# Basic usage

``traitar phenotype <in dir>  <sample file> from_nucleotides <out_dir> `` 

will trigger the standard workflow of traitAR, which is to predict open reading frames with Prodigal, annotate the coding sequences with Pfam families using HMMer and finally predict phenotypes from the models for the 67 traits. 

``traitar phenotype <in dir>  <sample file> from_genes <out_dir> `` 
 
assumes that gene prediction has been conducted already externally. In this case analysis will start with the Pfam annotation. If the output directory already exists, traitAR will offer to recompute or resume the individual analysis steps. This option is only available if the process is run interactively.
## Parallel usage
traitAR can benefit from parallel execution. The ``-c`` parameter sets the number of processes used e.g. ``-c 2`` for using two processes

``traitar phenotype <in dir>  <sample file> from_nucleotides out_dir -c 2`` 

This requires installing GNU parallel as noted above.
##Run traitAR with packaged sample data.
``traitar phenotype <traitar_dir>/data/sample_data <traitar_dir>/data/sample_data/samples.txt from_nucleotides <out_dir> -c 2`` will trigger phenotyping of *Listeria grayi DSM_20601* and *Listeria ivanovii WSLC3009*. Computation should be done within 5 minutes. You can find out ``<traitar_dir>`` by running

``python``

``>>> import traitar``

``>>> traitar.__path__``
## Output
traitAR provides the gene prediction results in ``<out_dir>/gene_prediction``, the Pfam annotation in ``<out_dir>/pfam_annotation`` and the phenotype prediction in``<out_dir>/phenotype prediction``. The phenotype prediction is summarized in heatmaps individually for the phyletic pattern classifier in ``heatmap_phypat.png``, for the phylogeny-aware classifier in ``heatmap_phypat_ggl.png`` and for both classifiers combined in ```heatmap_comb.png```. These heatmaps are based on tab separated text files e.g. ``predictions_majority-votes_combined.txt``. A negative prediction is encoded as 0, a prediction made only by the pure phyletic classifier as 1, one made by the phylogeny-aware classifier by 2 and a prediction supported by both algorithms as 3. ``predictions_flat_majority-votes_combined.txt`` provides a flat version of this table with one prediction per row. The expert user might also want to access the individual results for each algorithm in the respective sub folders ``phypat`` and ``phypat+GGL``.

# Web service
traitAR is also offered as a web service at
http://algbio.cs.uni-duesseldorf.de/webapps/wa-webservice/pipe.php?pr=phenolyzer 
