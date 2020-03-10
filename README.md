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
* [Authors](#7)
* [License](#8)

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

SPiP needs also to install 3 libraries, from the R console:
```R
install.packages("foreach")
install.packages("doParallel")
```

### Install Samtools<a id="3"></a>
SPiP will use samtools to get the DNA sequences.

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
| varType | substitution | Type of variant |
| ntChange | T>G | Nucleotides variation |
| ExonInfo | Intron 4 (1499) | Number and size of Exon/Intron |
| transcript | NM_007294 | Transcript (RefSeq) |
| gene | BRCA1 | Gene symbol (RefSeq) |
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

**--VCF**

+ Get the SPiP output in VCF format (v4.0)

```shell
##fileformat=VCFv4.0
##SPiP output v0.6
##SPiPCommand=/path/to/SPiPv0.6.r -I inputFile -O outputFile -s /path/to/samtools -f /path/to/genomeReference.fe --VCF
##assembly=GRCh37/hg19
##contig=<ID=chr1,length=249250621>
##contig=<ID=chr2,length=243199373>
##contig=<ID=chr3,length=198022430>
##contig=<ID=chr4,length=191154276>
##contig=<ID=chr5,length=180915260>
##contig=<ID=chr6,length=171115067>
##contig=<ID=chr7,length=159138663>
##contig=<ID=chr8,length=146364022>
##contig=<ID=chr9,length=141213431>
##contig=<ID=chr10,length=135534747>
##contig=<ID=chr11,length=135006516>
##contig=<ID=chr12,length=133851895>
##contig=<ID=chr13,length=115169878>
##contig=<ID=chr14,length=107349540>
##contig=<ID=chr15,length=102531392>
##contig=<ID=chr16,length=90354753>
##contig=<ID=chr17,length=81195210>
##contig=<ID=chr18,length=78077248>
##contig=<ID=chr19,length=59128983>
##contig=<ID=chr20,length=63025520>
##contig=<ID=chr21,length=48129895>
##contig=<ID=chrX,length=155270560>
##contig=<ID=chr22,length=51304566>
##contig=<ID=chrY,length=59373566>
##INFO=<ID=varID,Number=1,Type=String,Description="The name of variant (transcript:mutation)">
##INFO=<ID=Interpretation,Number=1,Type=String,Description="Overall prediction of SPiP">
##INFO=<ID=InterConfident,Number=1,Type=String,Description="Probability of splicing alteration with CI_95%, estimated from mutations 53,048 mutations">
##INFO=<ID=chr,Number=1,Type=String,Description="Chromosome number">
##INFO=<ID=strand,Number=1,Type=String,Description="Strand of the transcripts">
##INFO=<ID=varType,Number=1,Type=String,Description="Type of variant">
##INFO=<ID=ntChange,Number=1,Type=String,Description="Nucleotides variation">
##INFO=<ID=ExonInfo,Number=1,Type=String,Description="Number and size of Exon/Intron">
##INFO=<ID=transcript,Number=1,Type=String,Description="Transcript (RefSeq)">
##INFO=<ID=gene,Number=1,Type=String,Description="Gene symbol (RefSeq)">
##INFO=<ID=gNomen,Number=1,Type=String,Description="Genomic coordinates">
##INFO=<ID=seqPhysio,Number=1,Type=String,Description="(A, C, G, T)-sequence before the mutation">
##INFO=<ID=seqMutated,Number=1,Type=String,Description="(A, C, G, T)-sequence after the mutation">
##INFO=<ID=NearestSS,Number=1,Type=String,Description="Nearest splice site to the mutation">
##INFO=<ID=distSS,Number=1,Type=String,Description="Distance between the splice site and the mutation">
##INFO=<ID=RegType,Number=1,Type=String,Description="Type of region in the transcript, Exon/Intron">
##INFO=<ID=SPiCEproba,Number=1,Type=Float,Description="SPiCE score">
##INFO=<ID=SPiCEinter_2thr,Number=1,Type=String,Description="Classes of SPiCE (low, medium, high)">
##INFO=<ID=deltaMES,Number=1,Type=Float,Description="Delta score of MES">
##INFO=<ID=mutInPBarea,Number=1,Type=String,Description="Mutation in branch point">
##INFO=<ID=deltaESRscore,Number=1,Type=Float,Description="Score of deltaESRscore">
##INFO=<ID=posCryptMut,Number=1,Type=Integer,Description="Postion of mutated cryptic splice site">
##INFO=<ID=sstypeCryptMut,Number=1,Type=String,Description="Splice type of mutated cryptic splice site">
##INFO=<ID=probaCryptMut,Number=1,Type=Float,Description="Score of mutated cryptic splice site">
##INFO=<ID=classProbaCryptMut,Number=1,Type=String,Description="Use of mutated cryptic splice site (Yes/No)">
##INFO=<ID=nearestSStoCrypt,Number=1,Type=String,Description="Splice type of the nearest natural splice site to the mutated cryptic site">
##INFO=<ID=nearestPosSStoCrypt,Number=1,Type=Integer,Description="Position of the nearest natural splice site to the mutated cryptic site">
##INFO=<ID=nearestDistSStoCrypt,Number=1,Type=Integer,Description="Distance of the nearest natural splice site to the mutated cryptic site">
##INFO=<ID=posCryptWT,Number=1,Type=Integer,Description="Postion of wild-type cryptic splice site">
##INFO=<ID=probaCryptWT,Number=1,Type=Float,Description="Score of wild-type cryptic splice site">
##INFO=<ID=classProbaCryptWT,Number=1,Type=String,Description="Use of wild-type cryptic splice site (Yes/No)">
##INFO=<ID=posSSPhysio,Number=1,Type=Integer,Description="Position of the natural splice site (same splice type of cryptic site)">
##INFO=<ID=probaSSPhysio,Number=1,Type=Float,Description="Score of the natural splice site (same splice type of cryptic site)">
##INFO=<ID=classProbaSSPhysio,Number=1,Type=String,Description="Use of the natural splice site (same splice type of cryptic site) (Yes/No)">
##INFO=<ID=probaSSPhysioMut,Number=1,Type=Float,Description="Score of the natural splice site (same splice type of cryptic site) after the mutation">
##INFO=<ID=classProbaSSPhysioMut,Number=1,Type=String,Description="Use of the natural splice site (same splice type of cryptic site) after the mutation (Yes/No)">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	15765825	NM_007272:g.15765825:G>A	G	A	.	.	Interpretation="NTR"|InterConfident="00.04 % [00.02 % ; 00.08%]"|strand="+"|varType="substitution"|ntChange="G>A"|ExonInfo="Intron 1 (1795)"|transcript="NM_007272"|gene="CTRC"|NearestSS="donor"|DistSS="825"|RegType="DeepIntron"|SPiCEproba="0"|SPiCEinter_2thr="Outside SPiCE Interpretation"|deltaMES="0"|mutInPBarea="No"|deltaESRscore="NA"|posCryptMut="15765816"|sstypeCryptMut="Acc"|probaCryptMut="0.00206159394907144"|classProbaCryptMut="No"|nearestSStoCrypt="Don"|nearestPosSStoCrypt="15765000"|nearestDistSStoCrypt="816"|posCryptWT="15765816"|probaCryptWT="0.00161527498798199"|classProbaCryptWT="No"|posSSPhysio="15766795"|probaSSPhysio="0.0775463330795674"|classProbaSSPhysio="Yes"|probaSSPhysioMut="0.0775463330795674"|classProbaSSPhysioMut="Yes"
```

**--header**

+ Add the meta-column information to the file, to explain the significance of each SPiP column (only in text format output)

```shell
    ##ALT=<ID=*,Description=\"Represents allele(s) other than observed.\">
    ##INFO=<ID=Interpretation,Number=1,Type=String,Description=\"Overall prediction of SPiP\">
    ##INFO=<ID=InterConfident,Number=1,Type=String,Description=\"Probability of splicing alteration with CI_95%, estimated from mutations 53,048 mutations\">
    ##INFO=<ID=strand,Number=1,Type=String,Description=\"Strand of the transcripts\">
    ##INFO=<ID=varType,Number=1,Type=String,Description=\"Type of variant\">
    ##INFO=<ID=ntChange,Number=1,Type=String,Description=\"Nucleotides variation\">
    ##INFO=<ID=ExonInfo,Number=1,Type=String,Description=\"Number and size of Exon/Intron\">
    ##INFO=<ID=transcript,Number=1,Type=String,Description=\"Transcript (RefSeq)\">
    ##INFO=<ID=gene,Number=1,Type=String,Description=\"Gene symbol (RefSeq)\">
    ##INFO=<ID=NearestSS,Number=1,Type=String,Description=\"Nearest splice site to the mutation\">
    ##INFO=<ID=distSS,Number=1,Type=String,Description=\"Distance between the splice site and the mutation\">
    ##INFO=<ID=RegType,Number=1,Type=String,Description=\"Type of region in the transcript, Exon/Intron\">
    ##INFO=<ID=SPiCEproba,Number=1,Type=Float,Description=\"SPiCE score\">
    ##INFO=<ID=SPiCEinter_2thr,Number=1,Type=String,Description=\"Classes of SPiCE (low, medium, high)\">
    ##INFO=<ID=deltaMES,Number=1,Type=Float,Description=\"Delta score of MES\">
    ##INFO=<ID=mutInPBarea,Number=1,Type=String,Description=\"Mutation in branch point\">
    ##INFO=<ID=deltaESRscore,Number=1,Type=Float,Description=\"Score of deltaESRscore\">
    ##INFO=<ID=posCryptMut,Number=1,Type=Integer,Description=\"Postion of mutated cryptic splice site\">
    ##INFO=<ID=sstypeCryptMut,Number=1,Type=String,Description=\"Splice type of mutated cryptic splice site\">
    ##INFO=<ID=probaCryptMut,Number=1,Type=Float,Description=\"Score of mutated cryptic splice site\">
    ##INFO=<ID=classProbaCryptMut,Number=1,Type=String,Description=\"Use of mutated cryptic splice site (Yes/No)\">
    ##INFO=<ID=nearestSStoCrypt,Number=1,Type=String,Description=\"Splice type of the nearest natural splice site to the mutated cryptic site\">
    ##INFO=<ID=nearestPosSStoCrypt,Number=1,Type=Integer,Description=\"Position of the nearest natural splice site to the mutated cryptic site\">
    ##INFO=<ID=nearestDistSStoCrypt,Number=1,Type=Integer,Description=\"Distance of the nearest natural splice site to the mutated cryptic site\">
    ##INFO=<ID=posCryptWT,Number=1,Type=Integer,Description=\"Postion of wild-type cryptic splice site\">
    ##INFO=<ID=probaCryptWT,Number=1,Type=Float,Description=\"Score of wild-type cryptic splice site\">
    ##INFO=<ID=classProbaCryptWT,Number=1,Type=String,Description=\"Use of wild-type cryptic splice site (Yes/No)\">
    ##INFO=<ID=posSSPhysio,Number=1,Type=Integer,Description=\"Position of the natural splice site (same splice type of cryptic site)\">
    ##INFO=<ID=probaSSPhysio,Number=1,Type=Float,Description=\"Score of the natural splice site (same splice type of cryptic site)\">
    ##INFO=<ID=classProbaSSPhysio,Number=1,Type=String,Description=\"Use of the natural splice site (same splice type of cryptic site) (Yes/No)\">
    ##INFO=<ID=probaSSPhysioMut,Number=1,Type=Float,Description=\"Score of the natural splice site (same splice type of cryptic site) after the mutation\">
    ##INFO=<ID=classProbaSSPhysioMut,Number=1,Type=String,Description=\"Use of the natural splice site (same splice type of cryptic site) after the mutation (Yes/No)\">
```
## Authors <a id="7"></a>


* Raphael Leman - [raphaelleman](https://github.com/raphaelleman/ "tittle")
    * You can contact me at: r.leman@baclesse.unicancer.fr or raphael.leman@orange.fr

> **Cite as:** SPiP: a Splicing Prediction Pipeline addressing the diversity of splice alterations validated on a diagnostic set of 3,048 exonic and intronic variants
Raphaël Leman, Béatrice Parfait, Dominique Vidaud, Emmanuelle Girodon, Laurence Pacot, Gérald Legac, Chandran Ka, Claude Ferec, Yann Fichou, Céline Quesnelle, Etienne Muller, Dominique Vaur, Laurent Castera, Agathe Ricou, Hélène Tubeuf, Omar Soukarieh, Pascaline Gaildrat, Florence Riant, Marine Guillaud-Bataille, Sandrine M. Caputo, Virginie Caux-Moncoutier, Nadia Boutry-Kryza, Françoise Bonnet-Dorion, Ines Schultz, Maria Rossing, Michael T. Parsons, Amanda B. Spurdle, Thierry Frebourg, Alexandra Martins, Claude Houdayer, Sophie Krieger, [in preparation](https://doi.org/10.1093/bioinformatics/btz784 "tittle")

## License <a id="8"></a>


This project is licensed under the MIT License - see the [LICENSE](https://github.com/raphaelleman/SPiP/blob/master/LICENSE "tittle") file for details
