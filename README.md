# SPiP, **S**plicing **P**ipeline **P**rediction
![SPiP logo](https://github.com/raphaelleman/SPiP/blob/master/RefFiles/logoPipeline.gif)

---

SPiP is a decisional tree running a cascade of bioinformatics tools. Briefly, SPiP uses SPiCE tool for the consensus splice sites (donor and acceptor sites), MES for polypyrimidine tract between -13 and -20, BPP for branch point area between -18 and -44, an homemade score to research cryptic/de novo activation and ΔtESRseq for exonic splicing regulatory element until to 120 nt in exon

SPiP is available for Windows OS at https://sourceforge.net/projects/splicing-prediction-pipeline/

**Table**

* [Repository contents](#1)
* [Install SPiP](#2)
    * [Install Samtools](#3)
    * [Load the human genome](#4)
* [Run SPiP](#5)
    * [SPiP options](#6)

## Repository contents<a id="1"></a>

---

* SPiPv0.6: the SPiP scripts
* testCrypt.txt: an example of input data in text format
* testVar.vcf: an example of input data in vcf format
* *RefFiles*: folder where are the reference files use by SPiP

## Install SPiP<a id="2"></a>

---

To get SPiP from this repository, you can enter in the linux consoles:
```shell
git clone https://github.com/raphaelleman/SPiP
cd ./SPiP
```

SPiP needs also to install 2 libraries, from the R console:
```R
install.packages(“RCurl”)
install.packages(“parallel”)
```

For some version of Linux the installation of RCurl package requiers instalation of libcurl dependencies:
```shell
sudo apt-get update
sudo apt-get install libcurl4-openssl-dev
sudo apt-get install libcurl4-gnutls-dev
sudo apt-get install libcurl4-nss-dev
```

### Install Samtools<a id="3"></a>
In the main to optimize the run of SPiP you can also install samtools. SPiP will use samtools to get the DNA sequences.

```shell
cd /path/to/SPiP/
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar xfvj samtools-1.9.tar.bz2
cd ./samtools-1.9
./configure --prefix=/where/to/install
make
make install
```
For more information, please see the [samtools manual](http://www.htslib.org/doc/samtools.html "tittle")

### Load the human genome<a id="4"></a>

Here we use the [GENCODE](https://www.gencodegenes.org/ "tittle") database, with the hg19 example:
```shell
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz -O genomehg19.fa.gz
gunzip genomehg19.fa.gz
/path/to/samtools faidx genomehg19.fa
```

## Run SPiP<a id="5"></a>

---

you can get the different argument of SPiP by `Rscript /path/to/SPiPv0.6.r --help`

An example of SPiP run with test file [testCrypt.txt](http://gitlab.baclesse.fr/LEMRAP/spip/blob/master/testCrypt.txt "tittle"):

```shell
cd /path/to/SPiP/
Rscript ./SPiPv0.6.r -I ./testCrypt.txt -O ./outputTest.txt -s /path/to/samtools -f /path/to/genomehg19.fa
```

In this example SPiP will generate a text file "outputTest.txt" where the predictions will be save. The scheme of this output is:

| Column names | Example | Description |
|------------:|:--------:|:------------:|
| varID | NM_007294:c.213-6T>G | The variant id (Transcript:mutation) |
| Interpretation | Alter by SPiCE | The overall prediction |
| InterConfident | 92.9 % +/- 2.1 % | The risk that the variant impact splicing<br/>Estimated from collection of variant with *in vitro* RNA studies and frequent variant |
| chr | chr17 | Chromosome number |
| strand | - | Strand of the junction ('+': forward;<br/> '-':reverse) |
| gNomen | 41256979 | Genomic position of variant |
| seqPhysio | ACGG...AGGA | (A, C, G, T)-sequence before the mutation |
| seqMutated | ACGG...AGGA | (A, C, G, T)-sequence after the mutation |
| NearestSS | acceptor | The nearest natural splice site to the variant |
| distSS | -6 | Distance between the nearest splice site and the mutation |
| RegType | IntronCons | The type of region where located the variant |
| SPiCEproba | 1 | The SPiCE probability for variant in consensus splice site |
| SPiCEinter_2thr | high | The SPiCE classes (high/medium/low) |
| deltaMES | 0 | MES variation for variant in the polypyrimidine tract |
| mutInPBarea | No | If the mutation is located in branch point predicted by BPP tool |
| deltaESRscore | NA | ESR score variation for exonic variant |
| posCryptMut | 41256978 | Genomic position of cryptic splice site after the mutation |
| sstypeCryptMut | Acc | Splice type of cryptic splice site after the mutation |
| probaCryptMut | 0.000710404942432828 | Score of cryptic splice site after the mutation |
| classProbaCryptMut | No | Prediction of cryptic splice site after the mutation (Yes: used, No: Not used) |
| nearestSStoCrypt | Acc | Splice type of the nearest natural splice site |
| nearestPosSStoCrypt | 41256973 | Genomic position of the nearest natural splice site |
| nearestDistSStoCrypt | -5 | Distance between the cryptic site and the natural site |
| posCryptWT | 41256970 | Genomic position of cryptic splice site before the mutation |
| probaCryptWT | 4.89918764104143e-07 | Score of cryptic splice site before the mutation |
| classProbaCryptWT | No | Prediction of cryptic splice site before the mutation (Yes: used, No: Not used) |
| posSSPhysio | 41256973 | Genomic position of natural splice site that same splice site type of the mutated cryptic |
| probaSSPhysio | 0.00408919066993282 | Score of natural splice site that same splice site type of the mutated cryptic |
| classProbaSSPhysio | Yes | Prediction of natural splice site that same splice site type of the mutated cryptic (Yes: used, No: Not used) |
| probaSSPhysioMut | 1.74991364794327e-06 | Score of natural splice site that same splice site type of the mutated cryptic after the mutation |
| classProbaSSPhysioMut | No | Score of natural splice site that same splice site type of the mutated cryptic after the mutation (Yes: used, No: Not used) |

### SPiP options <a id="6"></a>

**-I, --input** /path/to/inputFile

+ list of variants file (.txt or .vcf). SPiP supports VCF version 4.1 or later (see example [testVar.vcf](https://github.com/raphaelleman/SPiP/blob/master/testVar.vcf "tittle")). The txt file must be tab-delimated and the column with mutation, in format Transcript:mutation, is indicated by 'varID' column name (see example [testCrypt.txt](https://github.com/raphaelleman/SPiP/blob/master/testCrypt.txt "tittle")).

**-O, --output** /path/to/outputFile

+ Name of ouput file (.txt). Directory to the output file (in text format)

**-g, --GenomeAssenbly** hg19

+ Genome assembly version (hg19 or hg38) [default= hg19]

**-s, --SamPath** /path/to/samtools]

+ Path to samtools, if you want to use Ensembl api keep this argument empty

**-f, --fastaGenome** /path/to/fastaGenome

+ Fasta file of genome used by samtools, see [Load the human genome](#4) section

**-t, --threads** N

+ Number of threads used for the calculation [default= 1]

**-l, --maxLines** N

+ Number of lines read in each time [default= 1000]

**--header**

+ Add the meta-column information to the file, to explain the significance of each SPiP column

```R
    ### SPiP output v0.6
    ## varID \tThe name of variant (transcript:mutation)
    ## Interpretation \tOverall prediction of SPiP
    ##    Alter by SPiCE \tAlteration of consensus splice site predicted by SPiCE (corresponding to classes \"medium\" and \"high\" of SPiCE)
    ##    Alter by MES (Poly TC) \tAlteration of polypyrimidine tract by MES (threshold of -15 %)
    ##    Alter BP \tAlteration of branch point predicted by BPP (variant in 4-mer of BP)
    ##    Alter ESR \tAlteration of ESR motifs predicted by deltaESRseq (threshold of 1.10)
    ##    Alter by create Cryptic/Exon \tCreation of new splice site
    ##    NTR \tNothing To Report, no alteration predicted
    ## InterConfident \tProbability of splicing alteration with CI_95%, estimated from mutations 65,955 mutations
    ## chr \tChromosome number
    ## strand \tStrand of the transcripts
    ## varType \tType of variant
    ## ntChange \tNucleotides variation
    ## ExonInfo \tNumber and size of Exon/Intron
    ## transcript \tTranscript (RefSeq)
    ## gene \tGene symbol (RefSeq)
    ## gNomen \tGenomic coordinates
    ## seqPhysio \t(A, C, G, T)-sequence before the mutation
    ## seqMutated \t(A, C, G, T)-sequence after the mutation
    ## NearestSS \tNearest splice site to the mutation
    ## distSS \tDistance between the splice site and the mutation
    ## RegType \tType of region in the transcript, Exon/Intron
    ##    Cons \tConsensus splice site (5\': -3; +6 / 3\': -12; +2)
    ##    PolyTC \tPolypyrimidine tract (-20; -13)
    ##    BP \tBranch point area (-44; -18)
    ## SPiCEproba \tSPiCE score
    ## SPiCEinter_2thr \tClasses of SPiCE (low, medium, high)
    ## deltaMES \tDelta score of MES",
    ## mutInPBarea \tMutation in branch point
    ## deltaESRscore \tScore of deltaESRscore
    ## posCryptMut \tPostion of mutated cryptic splice site
    ## sstypeCryptMut \tSplice type of mutated cryptic splice site
    ## probaCryptMut \tScore of mutated cryptic splice site
    ## classProbaCryptMut \tUse of mutated cryptic splice site (Yes/No)
    ## nearestSStoCrypt \tSplice type of the nearest natural splice site to the mutated cryptic site
    ## nearestPosSStoCrypt \tPosition of the nearest natural splice site to the mutated cryptic site
    ## nearestDistSStoCrypt \tDistance of the nearest natural splice site to the mutated cryptic site
    ## posCryptWT \tPostion of wild-type cryptic splice site
    ## probaCryptWT \tScore of wild-type cryptic splice site
    ## classProbaCryptWT \tUse of wild-type cryptic splice site (Yes/No)
    ## posSSPhysio \tPosition of the natural splice site (same splice type of cryptic site)
    ## probaSSPhysio \tScore of the natural splice site (same splice type of cryptic site)
    ## classProbaSSPhysio \tUse of the natural splice site (same splice type of cryptic site) (Yes/No)
    ## probaSSPhysioMut \tScore of the natural splice site (same splice type of cryptic site) after the mutation
    ## classProbaSSPhysioMut \tUse of the natural splice site (same splice type of cryptic site) after the mutation (Yes/No)
```
