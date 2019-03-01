# SPiP **S**plicing **P**ipeline **P**rediction

---

SPiP is a decisional tree running a cascade of bioinformatics tools. Briefly, SPiP uses SPiCE tool for the consensus splice sites (donor and acceptor sites), MES for polypyrimidine tract between -13 and -20, BPP for branch point area between -18 and -44, an homemade score to research cryptic/de novo activation and ΔtESRseq for exonic splicing regulatory element until to 120 nt in exon

## Repository contents

---

* SPiPv0.2: the SPiP scripts
* testCrypt.txt: an example of input data in text format
* testVar.vcf: an example of input data in vcf format
* *RefFiles*: folder where are the reference files use by SPiP

## Install SPiP

---

To get SPiP from this repository, you can enter in the linux consoles:

    git clone http://gitlab.baclesse.fr/LEMRAP/splicelauncherpipeline
    cd ./SPiP

SPiP needs also to install the Rcurl library, from the R console:

    install.packages(“Rcurl”)

In the main to optimize the run of SPiP you can also install samtools. SPiP will use samtools to get the DNA sequences.

### Install Samtools

    cd /path/to/SPiP
    wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
    tar xfvj samtools-1.9.tar.bz2
    cd ./samtools-1.9
    ./configure --prefix=/where/to/install
    make
    make install

For more information, please see the [samtools manual](http://www.htslib.org/doc/samtools.html "tittle")

### Load the human genome

Here we use the [GENCODE](https://www.gencodegenes.org/ "tittle") database, with the hg19 example:

    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz -O genomehg19.fa.gz
    gunzip genomehg19.fa.gz
    /path/to/samtools faidx genomehg19.fa

## Run SPiP

---

you can get the different argument of SPiP by `Rscript /path/to/SPiPv0.2.r --help`

An example of SPiP run with test file:

    cd /path/to/SPiP/
    Rscript ./SPiPv0.2.r -I ./testCrypt.txt -O ./outputTest.txt -s /path/to/samtools -f /path/to/fastaGenome

In this example SPiP will generate a text file "outputTest.txt" where the predictions will be save.
