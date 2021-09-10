# SPiP, **S**plicing **P**ipeline **P**rediction
![SPiP logo](https://github.com/raphaelleman/SPiP/blob/master/RefFiles/logoPipeline.gif)

---

SPiP is a randomForest model running a cascade of bioinformatics tools. Briefly, SPiP uses SPiCE tool for the consensus splice sites (donor and acceptor sites), MES for polypyrimidine tract between -13 and -20, BPP for branch point area between -18 and -44, an homemade score to research cryptic/de novo activation and ΔtESRseq for exonic splicing regulatory element until to 120 nt in exon

SPiP is available for Windows OS at https://sourceforge.net/projects/splicing-prediction-pipeline/

**Table**

* [Repository contents](#1)
* [Install SPiP](#2)
    * [Load the transcriptome files](#3)
* [Run SPiP](#4)
    * [SPiP options](#5)
* [Authors](#6)
* [License](#7)

## Repository contents<a id="1"></a>

---

* SPiPv2.1_main.r: the SPiP script
* testCrypt.txt: an example of input data in text format
* testVar.vcf: an example of input data in vcf format
* *RefFiles*: folder where are the reference files used by SPiP

## Install SPiP<a id="2"></a>

---

To get SPiP from this repository, you can enter in the linux consoles:
```shell
git clone https://github.com/raphaelleman/SPiP
cd ./SPiP
```

SPiP needs also to install 2 libraries, from the R console:
```R
install.packages("foreach")
install.packages("doParallel")
install.packages("randomForest")
```

### Load the transcriptome files<a id="3"></a>

you have to download frome sourcforge the RData files containing the transcripts sequences.
**hg19** assembly : [transcriptome_hg19.RData](https://sourceforge.net/projects/splicing-prediction-pipeline/files/transcriptome/transcriptome_hg19.RData/download "tittle")
**hg38** assembly : [transcriptome_hg38.RData](https://sourceforge.net/projects/splicing-prediction-pipeline/files/transcriptome/transcriptome_hg38.RData/download "tittle")

Put these files in `/path/to/SPiP/RefFiles/` or you can define it manually by the option `--transcriptome`.

NB: commands to regenerate these files are available in [getGenomeSequenceFromBSgenome.r](https://github.com/raphaelleman/SPiP/blob/master/RefFiles/getGenomeSequenceFromBSgenome.r "tittle")

## Run SPiP<a id="4"></a>

---

you can get the different argument of SPiP by `Rscript /path/to/SPiPv2.1_main.r --help`

An example of SPiP run with test file [testCrypt.txt](https://github.com/raphaelleman/SPiP/blob/master/testCrypt.txt "tittle"):

```shell
cd /path/to/SPiP/
Rscript ./SPiPv2.1_main.r -I ./testCrypt.txt -O ./outputTest.txt
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

### SPiP options <a id="5"></a>

**-I, --input** /path/to/inputFile

+ list of variants file (.txt or .vcf). SPiP supports VCF version 4.1 or later (see example [testVar.vcf](https://github.com/raphaelleman/SPiP/blob/master/testVar.vcf "tittle")). The txt file must be tab-delimated and the column with mutation, in format Transcript:mutation, is indicated by 'varID' column name (see example [testCrypt.txt](https://github.com/raphaelleman/SPiP/blob/master/testCrypt.txt "tittle")).

**-O, --output** /path/to/outputFile

+ Name of ouput file (.txt). Directory to the output file (in text format)

**-g, --GenomeAssenbly** hg19

+ Genome assembly version (hg19 or hg38) [default= hg19]

**-t, --threads** N

+ Number of threads used for the calculation [default= 1]

**-l, --maxLines** N

+ Number of lines read in each time [default= 1000]

**--verbose**

+ Show run process, *i.e.* displays progression bar tool

**--geneList** /path/to/geneList.txt

+ You can process analysis exclusively on a gene list, available only if VCF input

**--transcriptList** /path/to/transcriptList.txt

+ You can process analysis exclusively on a transcript list, available only if VCF input

**--transcriptome** /path/to/transcriptome_hgXX.RData

+ You can define where you have installed the file transcriptome_hgXX.RData if your file is not in /path/to/SPiP/RefFiles/

**--VCF**

+ Get the SPiP output in VCF format (v4.0)

```shell
# dynamic line modified in script : paste0("##SPiP output v",version)
# dynamic line modified in script : paste0("##SPiPCommand=",CMD)
## SPiP=altUsed|varID|Interpretation|InterConfident|SPiPscore|strand|gNomen|varType|ntChange|ExonInfo|exonSize|transcript|gene|NearestSS|DistSS|RegType|SPiCEproba|SPiCEinter_2thr|deltaMES|BP|mutInPBarea|deltaESRscore|posCryptMut|sstypeCryptMut|probaCryptMut|classProbaCryptMut|nearestSStoCrypt|nearestPosSStoCrypt|nearestDistSStoCrypt|posCryptWT|probaCryptWT|classProbaCryptWT|posSSPhysio|probaSSPhysio|classProbaSSPhysio|probaSSPhysioMut|classProbaSSPhysioMut
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	15765825	NM_007272:g.15765825:G>A	G	A	.	.	SPiP=A|NTR|00.04 % [00.02 % ; 00.08%]|+|substitution|G>A|Intron 1 (1795)|NM_007272|CTRC|donor|825|DeepIntron|0|Outside SPiCE Interpretation|0|No|NA|15765816|Acc|0.00206159394907144|No|Don|15765000|816|15765816|0.00161527498798199|No|15766795|0.0775463330795674|Yes|0.0775463330795674|Yes
```

## Authors <a id="6"></a>


* Raphael Leman - [raphaelleman](https://github.com/raphaelleman/ "tittle")
    * You can contact me at: r.leman@baclesse.unicancer.fr or raphael.leman@orange.fr

> **Cite as:** SPiP, a comprehensive Splicing Prediction Pipeline for massive detection of exonic and intronic variant effect on mRNA splicing.
**Raphaël Leman**, Béatrice Parfait, Dominique Vidaud, Emmanuelle Girodon, Laurence Pacot, Gérald Le Gac, Chandran Ka, Claude Ferec, Yann Fichou, Céline Quesnelle, MEtienne Muller, Dominique Vaur, Laurent Castera, Agathe Ricou, Hélène Tubeuf, Omar Soukarieh, Pascaline Gaildrat, Florence Riant, Marine Guillaud-Bataille, Sandrine M. Caputo, Virginie Caux-Moncoutier, Nadia Boutry-Kryza, Françoise Bonnet-Dorion, Ines Schultz, Maria Rossing, Louis Goldenberg, Olivier Quenez, Valentin Harter, Michael T. Parsons, Amanda B. Spurdle, Thierry Frébourg, Alexandra Martins, Claude Houdayer, Sophie Krieger, [in preparation](https://www.researchgate.net/publication/339817193_SPiP_a_Splicing_Prediction_Pipeline_addressing_the_diversity_of_splice_alterations_validated_on_a_diagnostic_set_of_3048_exonic_and_intronic_variants "tittle")

## License <a id="7"></a>


This project is licensed under the MIT License - see the [LICENSE](https://github.com/raphaelleman/SPiP/blob/master/LICENSE "tittle") file for details
