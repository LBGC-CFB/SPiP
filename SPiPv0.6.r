#!/usr/bin/env Rscript

#######################
# SPiP software
#######################
#import librairy
tryCatch({
library(parallel)
},
	error=function(cond) {
		message("Here's the original error message:")
		message(cond)
		message("*****You need to install \'parallel\' library\nInstall it by: install.pakages(\'parallel\')")
})
tryCatch({
library(foreach)
},
	error=function(cond) {
		message("Here's the original error message:")
		message(cond)
		message("*****You need to install \'foreach\' library\nInstall it by: install.pakages(\'foreach\')")
})
tryCatch({
library(doParallel)
},
	error=function(cond) {
		message("Here's the original error message:")
		message(cond)
		message("*****You need to install \'doParallel\' library\nInstall it by: install.pakages(\'doParallel\')")
})

message("
      _.-'''-,
    .'        `\\
   /           /
  /      .--^_^
  |     /  C ,,\\
  |    |   \\  _.)
   \\   |   /  \\
    '-, \\./ \\)\\)
       `-/   );/
_________''--'-'________________
")

# Env variables
options(scipen=50)
samPath=NULL
fastaFile=NULL
threads = 1
genome="hg19"
maxLines = 1000
printHead = FALSE
printVCF = FALSE

#SPiP arguments
helpMessage=paste("Usage: SPiPv0.6.r\n
    Mandatory \n
        -I, --input /path/to/inputFile\t\tlist of variants file (.txt or .vcf)
        -O, --output /path/to/outputFile\t\tName of ouput file (.txt)
        -s, --SamPath /path/to/samtools]\t\tPath to samtools executable
        -f, --fastaGenome /path/to/fastaGenome\t\tfasta file of genome used by samtools\n
    Genome options \n
        -g, --GenomeAssenbly hg19\t\tGenome assembly version (hg19 or hg38) [default= ",genome,"] \n
    Parallel options \n
        -t, --threads N\t\tNumber of threads used for the calculation [default= ",threads,"]
        -l, --maxLines N\t\tNumber of lines read in each time [default= ",maxLines,"]\n
    Other options\n
        --VCF\t\tPrint output in vcf format
        --header \t\tPrint meta-header info
    -h, --help\t\tPrint this help message and exit\n
   You could : Rscript SPiPv0.6.r -I ./testCrypt.txt -O ./outTestCrypt.txt",sep="")

#get script argument
argsFull <- commandArgs()

Rscript <- argsFull[1]

scriptPath=dirname(normalizePath(sub("--file=","",argsFull[substr(argsFull,1,7)=="--file="])))
if (length(which(argsFull=="--args"))==0){message(helpMessage);q(save = "no")}

args = argsFull[(which(argsFull=="--args")+1):length(argsFull)]

if (length(args)<4){message(helpMessage);stop("Not enought arguments")}

i=1
while (i <= length(args)){
    if(args[i]=="-I"|args[i]=="--input"){
        inputFile=args[i+1];i = i+2
    }else if(args[i]=="-O"|args[i]=="--output"){
        outputFile=args[i+1];i = i+2
    }else if(args[i]=="-g"|args[i]=="--GenomeAssenbly"){
        genome=args[i+1];i = i+2
    }else if(args[i]=="-s"|args[i]=="--SamPath"){
        samPath=args[i+1];i = i+2
    }else if(args[i]=="-f"|args[i]=="--fastaGenome"){
        fastaFile=args[i+1];i = i+2
    }else if(args[i]=="-t"|args[i]=="--threads"){
        threads= as.numeric(args[i+1]);i = i+2
    }else if(args[i]=="-l"|args[i]=="--maxLines"){
        maxLines=as.numeric(args[i+1]);i = i+2
    }else if(args[i]=="--header"){
        printHead=TRUE;i = i+1
    }else if(args[i]=="--VCF"){
        printVCF=TRUE;i = i+1
    }else if(args[i]=="-h"|args[i]=="--help"){
        message(helpMessage);stop()
    }else{
        message(paste("********Unknown option:",args[i],"\n"));message(helpMessage);stop()
    }
}

#Other argument
if(genome!="hg19" & genome!="hg38"){
	message("###########################")
	message("#Define the assembly genome version (hg19 or hg38)")
	message("###########################")
	message(helpMessage)
	stop()
}

if(is.null(fastaFile)|is.null(samPath)){
	message("###########################")
	message("#No fasta file for samtools")
	message("###########################")
    message(helpMessage)
    stop()
}

fastaIndex = paste0(fastaFile,".fai")
if(!file.exists(fastaIndex)){
    message("###########################")
    message("#No index file for samtools")
    message("#Use: /path/to/samtools faidx /path/to/genome.fa")
    message("###########################")
    message(helpMessage)
    stop()
}

registerDoParallel(threads)
CMD = paste(normalizePath(sub("--file=","",argsFull[substr(argsFull,1,7)=="--file="])),
        " --input ", normalizePath(inputFile),
        " --output ", outputFile,
        " --SamPath ", normalizePath(samPath),
        " --fastaGenome ", normalizePath(fastaFile),
        " --GenomeAssenbly ", genome,
        " --threads ", threads,
        " --maxLines ", maxLines,
        if(printVCF){" --VCF "},
        if(printHead){" --header "},sep="")

headerHelp_hg19 = c("##fileformat=VCFv4.0",
        "##SPiP output v0.6",
        paste("##SPiPCommand=",CMD,sep=""),
        "##assembly=GRCh37/hg19",
        "##contig=<ID=chr1,length=249250621>",
        "##contig=<ID=chr2,length=243199373>",
        "##contig=<ID=chr3,length=198022430>",
        "##contig=<ID=chr4,length=191154276>",
        "##contig=<ID=chr5,length=180915260>",
        "##contig=<ID=chr6,length=171115067>",
        "##contig=<ID=chr7,length=159138663>",
        "##contig=<ID=chr8,length=146364022>",
        "##contig=<ID=chr9,length=141213431>",
        "##contig=<ID=chr10,length=135534747>",
        "##contig=<ID=chr11,length=135006516>",
        "##contig=<ID=chr12,length=133851895>",
        "##contig=<ID=chr13,length=115169878>",
        "##contig=<ID=chr14,length=107349540>",
        "##contig=<ID=chr15,length=102531392>",
        "##contig=<ID=chr16,length=90354753>",
        "##contig=<ID=chr17,length=81195210>",
        "##contig=<ID=chr18,length=78077248>",
        "##contig=<ID=chr19,length=59128983>",
        "##contig=<ID=chr20,length=63025520>",
        "##contig=<ID=chr21,length=48129895>",
        "##contig=<ID=chrX,length=155270560>",
        "##contig=<ID=chr22,length=51304566>",
        "##contig=<ID=chrY,length=59373566>")

headerHelp_hg38 = c("##fileformat=VCFv4.0",
        "##SPiP output v0.6",
        paste("##SPiPCommand=",CMD,sep=""),
        "##assembly=GRCh38/hg38",
        "##contig=<ID=chr1,length=248956422>",
        "##contig=<ID=chr2,length=242193529>",
        "##contig=<ID=chr3,length=198295559>",
        "##contig=<ID=chr4,length=190214555>",
        "##contig=<ID=chr5,length=181538259>",
        "##contig=<ID=chr6,length=170805979>",
        "##contig=<ID=chr7,length=159345973>",
        "##contig=<ID=chr8,length=145138636>",
        "##contig=<ID=chr9,length=138394717>",
        "##contig=<ID=chr10,length=133797422>",
        "##contig=<ID=chr11,length=135086622>",
        "##contig=<ID=chr13,length=114364328>",
        "##contig=<ID=chr14,length=107043718>",
        "##contig=<ID=chr12,length=133275309>",
        "##contig=<ID=chr15,length=101991189>",
        "##contig=<ID=chr16,length=90338345>",
        "##contig=<ID=chr17,length=83257441>",
        "##contig=<ID=chr18,length=80373285>",
        "##contig=<ID=chr19,length=58617616>",
        "##contig=<ID=chr20,length=64444167>",
        "##contig=<ID=chr21,length=46709983>",
        "##contig=<ID=chr22,length=50818468>",
        "##contig=<ID=chrX,length=156040895>",
        "##contig=<ID=chrY,length=57227415>")

headerHelp = c("##ALT=<ID=*,Description=\"Represents allele(s) other than observed.\">",
                "##INFO=<ID=Interpretation,Number=1,Type=String,Description=\"Overall prediction of SPiP\">",
                "##INFO=<ID=InterConfident,Number=1,Type=String,Description=\"Probability of splicing alteration with CI_95%, estimated from mutations 53,048 mutations\">",
                "##INFO=<ID=strand,Number=1,Type=String,Description=\"Strand of the transcripts\">",
                "##INFO=<ID=varType,Number=1,Type=String,Description=\"Type of variant\">",
                "##INFO=<ID=ntChange,Number=1,Type=String,Description=\"Nucleotides variation\">",
                "##INFO=<ID=ExonInfo,Number=1,Type=String,Description=\"Number and size of Exon/Intron\">",
                "##INFO=<ID=transcript,Number=1,Type=String,Description=\"Transcript (RefSeq)\">",
                "##INFO=<ID=gene,Number=1,Type=String,Description=\"Gene symbol (RefSeq)\">",
                "##INFO=<ID=NearestSS,Number=1,Type=String,Description=\"Nearest splice site to the mutation\">",
                "##INFO=<ID=distSS,Number=1,Type=String,Description=\"Distance between the splice site and the mutation\">",
                "##INFO=<ID=RegType,Number=1,Type=String,Description=\"Type of region in the transcript, Exon/Intron\">",
                "##INFO=<ID=SPiCEproba,Number=1,Type=Float,Description=\"SPiCE score\">",
                "##INFO=<ID=SPiCEinter_2thr,Number=1,Type=String,Description=\"Classes of SPiCE (low, medium, high)\">",
                "##INFO=<ID=deltaMES,Number=1,Type=Float,Description=\"Delta score of MES\">",
                "##INFO=<ID=mutInPBarea,Number=1,Type=String,Description=\"Mutation in branch point\">",
                "##INFO=<ID=deltaESRscore,Number=1,Type=Float,Description=\"Score of deltaESRscore\">",
                "##INFO=<ID=posCryptMut,Number=1,Type=Integer,Description=\"Postion of mutated cryptic splice site\">",
                "##INFO=<ID=sstypeCryptMut,Number=1,Type=String,Description=\"Splice type of mutated cryptic splice site\">",
                "##INFO=<ID=probaCryptMut,Number=1,Type=Float,Description=\"Score of mutated cryptic splice site\">",
                "##INFO=<ID=classProbaCryptMut,Number=1,Type=String,Description=\"Use of mutated cryptic splice site (Yes/No)\">",
                "##INFO=<ID=nearestSStoCrypt,Number=1,Type=String,Description=\"Splice type of the nearest natural splice site to the mutated cryptic site\">",
                "##INFO=<ID=nearestPosSStoCrypt,Number=1,Type=Integer,Description=\"Position of the nearest natural splice site to the mutated cryptic site\">",
                "##INFO=<ID=nearestDistSStoCrypt,Number=1,Type=Integer,Description=\"Distance of the nearest natural splice site to the mutated cryptic site\">",
                "##INFO=<ID=posCryptWT,Number=1,Type=Integer,Description=\"Postion of wild-type cryptic splice site\">",
                "##INFO=<ID=probaCryptWT,Number=1,Type=Float,Description=\"Score of wild-type cryptic splice site\">",
                "##INFO=<ID=classProbaCryptWT,Number=1,Type=String,Description=\"Use of wild-type cryptic splice site (Yes/No)\">",
                "##INFO=<ID=posSSPhysio,Number=1,Type=Integer,Description=\"Position of the natural splice site (same splice type of cryptic site)\">",
                "##INFO=<ID=probaSSPhysio,Number=1,Type=Float,Description=\"Score of the natural splice site (same splice type of cryptic site)\">",
                "##INFO=<ID=classProbaSSPhysio,Number=1,Type=String,Description=\"Use of the natural splice site (same splice type of cryptic site) (Yes/No)\">",
                "##INFO=<ID=probaSSPhysioMut,Number=1,Type=Float,Description=\"Score of the natural splice site (same splice type of cryptic site) after the mutation\">",
                "##INFO=<ID=classProbaSSPhysioMut,Number=1,Type=String,Description=\"Use of the natural splice site (same splice type of cryptic site) after the mutation (Yes/No)\">"
)

fileFormat = tolower(substr(basename(inputFile),nchar(basename(inputFile))-2,nchar(basename(inputFile))))
fileFormatOut = tolower(substr(outputFile,nchar(outputFile)-2,nchar(outputFile)))

if(fileFormat!="txt" & fileFormat!="vcf"){
    message("###########################")
    message("#Incorrect format of input, please try again with a txt or vcf file")
    message("###########################")
    message(helpMessage)
    stop()
}

if(fileFormatOut=="vcf"){printVCF = TRUE}

output<-file(outputFile,"w")

if(printHead & !printVCF){
    writeLines(headerHelp,con = output,sep="\n")
}else if(printVCF){
    writeLines(c(if(genome=="hg19"){headerHelp_hg19}else{headerHelp_hg38},headerHelp),con = output,sep="\n")
}

flush(output)
close(output)

message("##################")
message("#Your command:")
message("##################")
message(CMD)

#Get Ref files
inputref = paste(scriptPath, "/RefFiles",sep="")
message("Check RefSeq database...")
if(!file.exists(paste(inputref,"/dataRefSeqhg19.RData",sep="")) & !file.exists(paste(inputref,"/dataRefSeqhg38.RData",sep=""))){
    currentWD = getwd()
    setwd(scriptPath)
	message("Create RefSeq database...")
	source(paste(inputref,"/getRefSeqDatabase.r",sep=""),local =TRUE)
    setwd(currentWD)
}
message("Load RefSeq database...")
load(paste(inputref, "/dataRefSeq",genome,".RData",sep=""))
load(paste(inputref, "/RefFiles.RData",sep=""))

mint_GT=sum(as.numeric(as.vector(sub("Min.   :","",summary(ref_score_GT)[1,]))))
maxt_GT=sum(as.numeric(as.vector(sub("Max.   :","",summary(ref_score_GT)[6,]))))
mint_GC=sum(as.numeric(as.vector(sub("Min.   :","",summary(ref_score_GC)[1,]))))
maxt_GC=sum(as.numeric(as.vector(sub("Max.   :","",summary(ref_score_GC)[6,]))))
mint1=as.numeric(as.vector(sub("Min.   :","",summary(ref_score_AG)[1,1:10])))
maxt1=as.numeric(as.vector(sub("Max.   :","",summary(ref_score_AG)[6,1:10])))

maxt1=maxt1[order(maxt1,decreasing=T)]
maxt1=maxt1[1:8]
maxt1=sum(maxt1)

mint1=mint1[order(mint1,decreasing=F)]
mint1=mint1[1:8]
mint1=sum(mint1)

mint2=sum(as.numeric(as.vector(sub("Min.   :","",summary(ref_score_AG)[1,12:15]))))
maxt2=sum(as.numeric(as.vector(sub("Max.   :","",summary(ref_score_AG)[6,12:15]))))

i_score=NULL
i_score1=NULL
i_score2=NULL

dataESR$hexamer = as.character(dataESR$hexamer)
ESRmotif = dataESR$hexamer[dataESR$Assignment!="N"]
LEIsc_valuesWA = dataESR$LEIsc_valuesWA
LEIsc_valuesHA = dataESR$LEIsc_valuesHA
LEIsc_valuesHM = dataESR$LEIsc_valuesHM
LEIsc_valuesWD = dataESR$LEIsc_valuesWD
LEIsc_valuesHD = dataESR$LEIsc_valuesHD
ESRlistScore = dataESR$ESEseq_or_ESSseqscore

names(LEIsc_valuesWA) <- dataESR$hexamer
names(LEIsc_valuesHA) <- dataESR$hexamer
names(LEIsc_valuesHM) <- dataESR$hexamer
names(LEIsc_valuesWD) <- dataESR$hexamer
names(LEIsc_valuesHD) <- dataESR$hexamer
names(ESRlistScore) <- dataESR$hexamer
indAcc = c(1:11,17:62,68:95)
indDon = c(1:28,34:74,80:95)

me2x5 = ME2x5$V1
names(me2x5) <- as.character(ME2x5$V1.1)

inverseDic <- data.frame(V1 = c('T','G','C','A','N'),row.names = c('A','C','G','T','N'))
inverseDic$V1 <- as.character(inverseDic$V1)

contigToChr <- function(text){
    text1=unlist(strsplit(text, ".", fixed = TRUE))[1]
    chr=text
    if(substr(text1,1,2)=="NC"){
        chr=paste('chr',as.numeric(substr(text1,4,nchar(text1))),sep="")
        if(chr=="chr23"){chr="chrX"}else if(chr=="chr24"){chr="chrY"}
    }
    return(chr)
}

chrIndex = read.table(fastaIndex,sep="\t",header=FALSE,colClasses = "character")
chrIndex = chrIndex[,1]
if(substr(chrIndex[1],1,3)=="chr"){
    names(chrIndex) = chrIndex
}else if(nchar(chrIndex[1])<=2){
    names(chrIndex) = paste0("chr",chrIndex)
}else if(substr(chrIndex[1],1,3)=="NC_"){
    chrConvert = NULL
    for(i in chrIndex){chrConvert = c(chrConvert,contigToChr(i))}
    names(chrIndex) = chrConvert
}

#functions used in this pipeline

getPosSSphysio <- function(transcrit){

	if(dim(dataRefSeq[dataRefSeq$V4==transcrit,])[1]==0){
		paste("I don't find the transcript:",transcrit,"in the Refseq database")
	}else{

		chr <<- as.character(dataRefSeq[dataRefSeq$V4==transcrit,1])
		sens <<- as.character(dataRefSeq[dataRefSeq$V4==transcrit,6])

	if(sens=="+"){
		posStart=dataRefSeq[dataRefSeq$V4==transcrit,2]
		tailleCum=dataRefSeq[dataRefSeq$V4==transcrit,12]
		tailleCum=strsplit(as.character(tailleCum),split=",")
		tailleCum=as.numeric(unlist(tailleCum))
		posAcc <<- posStart+tailleCum
		taille=dataRefSeq[dataRefSeq$V4==transcrit,11]
		taille=strsplit(as.character(taille),split=",")
		taille=as.numeric(unlist(taille))
		posDon <<- posAcc+taille
		if(length(posDon)>1 & length(posAcc)>1){
			posDon <<- posDon[-length(posDon)]
			posAcc <<- posAcc[-1]
		}
	}else if(sens=="-"){
		posEnd=dataRefSeq[dataRefSeq$V4==transcrit,2]
		tailleCum=dataRefSeq[dataRefSeq$V4==transcrit,12]
		tailleCum=strsplit(as.character(tailleCum),split=",")
		tailleCum=as.numeric(unlist(tailleCum))
		posDon <<- posEnd+tailleCum
		taille=dataRefSeq[dataRefSeq$V4==transcrit,11]
		taille=strsplit(as.character(taille),split=",")
		taille=as.numeric(unlist(taille))
		posAcc <<- posDon+taille
		if(length(posDon)>1 & length(posAcc)>1){
			posDon <<- posDon[-1]
			posAcc <<- posAcc[-length(posAcc)]
		}
	}
	}
}

getNearestPos <- function(sens, varPos ,posDon, posAcc){
    varPos1 = varPos[1]
	distSS2 = NULL
	varPos2 = NULL
	if(length(varPos)==2){
		varPos2 = varPos[2]
	}
	posDon = posDon[order(posDon)]
	posAcc = posAcc[order(posAcc)]

	minPosDon = min(abs(posDon-varPos1 ))
	minPosAcc = min(abs(posAcc-varPos1 ))

	if(minPosDon<minPosAcc){
		SstypePhy <<- "donor"
		if(length( posDon[posDon==(varPos1 + min(abs(posDon-varPos1 )))])==0){
			nearestPosDon <-posDon[posDon==(varPos1 - min(abs(posDon-varPos1 )))]
		}else{
			nearestPosDon <-posDon[posDon==(varPos1 + min(abs(posDon-varPos1 )))]
		}

		if(sens=="+"){
			distSS1 <<- varPos1 - nearestPosDon
			if(length(varPos)==2){
				distSS2 = varPos2 - nearestPosDon
			}

			if(varPos1 <= nearestPosDon){
				nearestPosAcc = posAcc[which(posDon==nearestPosDon)-1]
				distSS1 <<- distSS1 -1
                if(abs(distSS1)<=120){
                    RegType <<- "ExonESR"
                }else{
                    RegType <<- "Exon"
                }
				if(length(varPos)==2){
					distSS2 = distSS2 -1
				}
			}else{
				nearestPosAcc = posAcc[which(posDon==nearestPosDon)]
                if(abs(distSS1)<=150){
                    RegType <<- "Intron"
                }else{
                    RegType <<- "DeepIntron"
                }
			}
		}else{
			distSS1 <<- nearestPosDon - varPos1+1
			if(length(varPos)==2){
				distSS2 = nearestPosDon - varPos2+1
			}
			if(varPos1 <= nearestPosDon){
				nearestPosAcc = posAcc[which(posDon==nearestPosDon)]
                if(abs(distSS1)<=150){
                    RegType <<- "Intron"
                }else{
                    RegType <<- "DeepIntron"
                }
			}else{
				nearestPosAcc = posAcc[which(posDon==nearestPosDon)+1]
				distSS1 <<- distSS1 -1
                if(abs(distSS1)<=120){
                    RegType <<- "ExonESR"
                }else{
                    RegType <<- "Exon"
                }
				if(length(varPos)==2){
					distSS2 = distSS2 -1
				}
			}
		}
		if(is.null(distSS2)){
			if(distSS1>=(-3) & distSS1<=(6)){
				RegType <<- paste(RegType,"Cons",sep="")
			}
		}else{
			if((distSS1>=(-3) & distSS1<=6) | (distSS2>=(-3) & distSS2<=6) | (distSS1<(-3) & distSS2>6)){
				RegType <<- paste(RegType,"Cons",sep="")
			}
		}
	}else{
		SstypePhy <<- "acceptor"
		if(length( posAcc[posAcc==(varPos1 + min(abs(posAcc-varPos1 )))])==0){
			nearestPosAcc <-posAcc[posAcc==(varPos1 - min(abs(posAcc-varPos1 )))]
		}else{
			nearestPosAcc <-posAcc[posAcc==(varPos1 + min(abs(posAcc-varPos1 )))]
		}
		if(sens=="+"){
			distSS1 <<- varPos1 - nearestPosAcc-1
			if(length(varPos)==2){
				distSS2 = varPos2 - nearestPosAcc-1
			}
			if(varPos1 <= nearestPosAcc){
				nearestPosDon = posDon[which(posAcc==nearestPosAcc)]
                if(abs(distSS1)<=150){
                    RegType <<- "Intron"
                }else{
                    RegType <<- "DeepIntron"
                }
			}else{
				nearestPosDon = posDon[which(posAcc==nearestPosAcc)+1]
                distSS1 <<- distSS1 +1
                if(abs(distSS1)<=120){
                    RegType <<- "ExonESR"
                }else{
                    RegType <<- "Exon"
                }
				if(length(varPos)==2){
					distSS2 = distSS2 + 1
				}
			}
		}else{
			distSS1 <<- nearestPosAcc - varPos1
			if(length(varPos)==2){
				distSS2 = nearestPosAcc - varPos2
			}
			if(varPos1 <= nearestPosAcc){
				nearestPosDon = posDon[which(posAcc==nearestPosAcc)-1]
                distSS1 <<- distSS1 +1
                if(abs(distSS1)<=120){
                    RegType <<- "ExonESR"
                }else{
                    RegType <<- "Exon"
                }
				if(length(varPos)==2){
					distSS2 = distSS2 + 1
				}
			}else{
				nearestPosDon = posDon[which(posAcc==nearestPosAcc)]
                if(abs(distSS1)<=150){
                    RegType <<- "Intron"
                }else{
                    RegType <<- "DeepIntron"
                }
			}
		}
		if(is.null(distSS2)){
			if(distSS1>=(-12) & distSS1<=2){
				RegType <<- paste(RegType,"Cons",sep="")
			}
			if (distSS1>=(-20) & distSS1<=(-13)){
				RegType <<- paste(RegType,"PolyTC",sep="")
			}
			if (distSS1>=(-44) & distSS1<=(-18)){
				RegType <<- paste(RegType,"BP",sep="")
			}
		}else{
			if((distSS1>=(-12) & distSS1<=2) | (distSS2>=(-12) & distSS2<=2) | (distSS1<(-12) & distSS2>2)){
				RegType <<- paste(RegType,"Cons",sep="")
			}
			if (distSS1>=(-20) & distSS1<=(-13) | (distSS2>=(-20) & distSS2<=(-13)) | (distSS1<(-20) & distSS2>(-13))){
				RegType <<- paste(RegType,"PolyTC",sep="")
			}
			if (distSS1>=(-44) & distSS1<=(-18) | distSS2>=(-44) & distSS2<=(-18) | (distSS1<(-44) & distSS2>(-18))){
				RegType <<- paste(RegType,"BP",sep="")
			}
		}
	}

	distSS <<- distSS1
	distSS2 <<- distSS2
	nearestPosAll <<- as.numeric(na.omit(c(nearestPosDon,nearestPosAcc)))
}

getNearestPosCrypt <- function(sens, cryptPos ,posDon, posAcc){
	nearestPosDon = 0
	nearestPosAcc = 0
	posDon = posDon[order(posDon)]
	posAcc = posAcc[order(posAcc)]
	minPosDon = min(abs(posDon-cryptPos))
	minPosAcc = min(abs(posAcc-cryptPos))

	if(minPosDon<minPosAcc){
		SstypePhyCrypt = "Don"
		if(length( posDon[posDon==(cryptPos + min(abs(posDon-cryptPos )))])==0){
			nearestPosDon <-posDon[posDon==(cryptPos - min(abs(posDon-cryptPos )))]
		}else{
			nearestPosDon <-posDon[posDon==(cryptPos + min(abs(posDon-cryptPos )))]
		}
		if(sens=="+"){
			distSSc = cryptPos - nearestPosDon
			if(cryptPos <= nearestPosDon){
				distSSc = distSSc -1
			}
		}else{
			distSSc = nearestPosDon - cryptPos+1
			if(cryptPos > nearestPosDon){
				distSSc = distSSc -1
			}
		}
	}else{
		SstypePhyCrypt = "Acc"
		if(length( posAcc[posAcc==(cryptPos + min(abs(posAcc-cryptPos )))])==0){
			nearestPosAcc <-posAcc[posAcc==(cryptPos - min(abs(posAcc-cryptPos )))]
		}else{
			nearestPosAcc <-posAcc[posAcc==(cryptPos + min(abs(posAcc-cryptPos )))]
		}
		if(sens=="+"){
			distSSc = cryptPos - nearestPosAcc-1
			if(cryptPos > nearestPosAcc){
				nearestPosDon = posDon[which(posAcc==nearestPosAcc)+1]
				distSSc = distSSc + 1
			}
		}else{
			distSSc = nearestPosAcc - cryptPos
			if(cryptPos <= nearestPosAcc){
				distSSc = distSSc + 1
			}
		}
	}
	distSScrypt <<- distSSc
	SstypePhyCrypt <<- SstypePhyCrypt
	if(SstypePhyCrypt=="Acc"){
		nearestPosPhyCrypt <<- nearestPosAcc
	}else{
		nearestPosPhyCrypt <<- nearestPosDon
	}
}

getExonInfo <- function(transcrit,posVar){

	sens = as.character(dataRefSeq[dataRefSeq$V4==transcrit,6])
	posStart=dataRefSeq[dataRefSeq$V4==transcrit,2]
	tailleExon = as.numeric(unlist(strsplit(as.character(dataRefSeq[dataRefSeq$V4==transcrit,11]),",")))
	tailleCum=dataRefSeq[dataRefSeq$V4==transcrit,12]
	tailleCum=strsplit(as.character(tailleCum),split=",")
	tailleCum=as.numeric(unlist(tailleCum))

	if(sens=="+"){

		gCDSstart = dataRefSeq[dataRefSeq$V4==transcrit,7]
		gCDSend = dataRefSeq[dataRefSeq$V4==transcrit,8]
		posAcc = posStart+tailleCum + 1
		posDon = posAcc+tailleExon - 1

		dataConvert=data.frame(idEx=c(1:length(tailleExon )),lenEx=tailleExon,gStart=posAcc,gEnd=posDon)

		if(nrow(dataConvert[dataConvert$gStart<=posVar & dataConvert$gEnd>=posVar,])>=1){
			ExonInfo = paste("Exon ",
					dataConvert$idEx[dataConvert$gStart<=posVar & dataConvert$gEnd>=posVar]," (",
					dataConvert$lenEx[dataConvert$gStart<=posVar & dataConvert$gEnd>=posVar],")",
					sep="")
		}else{
			ExonInfo = paste("Intron ",
							dataConvert$idEx[which(dataConvert$gEnd>=posVar)[1]-1]," (",
							abs(dataConvert$gEnd[which(dataConvert$gEnd>=posVar)[1]-1]-
								dataConvert$gStart[which(dataConvert$gEnd>=posVar)[1]])-1,")",
							sep="")
		}
	}else if(sens=="-"){

		gCDSstart = dataRefSeq[dataRefSeq$V4==transcrit,8]
		gCDSend = dataRefSeq[dataRefSeq$V4==transcrit,7]
		posDon = posStart+tailleCum + 1
		posAcc = posDon+tailleExon - 1

		dataConvert=data.frame(idEx=c(length(tailleExon):1),lenEx=tailleExon,gStart=posAcc,gEnd=posDon)

		if(nrow(dataConvert[dataConvert$gStart>=posVar & dataConvert$gEnd<=posVar,])>=1){
			ExonInfo = paste("Exon ",
					dataConvert$idEx[dataConvert$gStart>=posVar & dataConvert$gEnd<=posVar]," (",
					dataConvert$lenEx[dataConvert$gStart>=posVar & dataConvert$gEnd<=posVar],")",
					sep="")
		}else{
			ExonInfo = paste("Intron ",
							dataConvert$idEx[which(dataConvert$gEnd>=posVar)[1]]," (",
							abs(dataConvert$gEnd[which(dataConvert$gEnd>=posVar)[1]]-
								dataConvert$gStart[which(dataConvert$gEnd>=posVar)[1]-1])-1,")",
							sep="")
		}
	}
	return(ExonInfo)
}

convertcNomenIngNomen <- function(transcrit,posVar){

	sens = as.character(dataRefSeq[dataRefSeq$V4==transcrit,6])
	posStart=dataRefSeq[dataRefSeq$V4==transcrit,2]
	tailleExon = as.numeric(unlist(strsplit(as.character(dataRefSeq[dataRefSeq$V4==transcrit,11]),",")))
	tailleCum=dataRefSeq[dataRefSeq$V4==transcrit,12]
	tailleCum=strsplit(as.character(tailleCum),split=",")
	tailleCum=as.numeric(unlist(tailleCum))

	if(sens=="+"){

		gCDSstart = dataRefSeq[dataRefSeq$V4==transcrit,7]
		gCDSend = dataRefSeq[dataRefSeq$V4==transcrit,8]

		posAcc = posStart+tailleCum + 1
		posDon = posAcc+tailleExon - 1

		dataConvert=data.frame(idEx=c(1:length(tailleExon )),lenEx=tailleExon,gStart=posAcc,gEnd=posDon,cStart=0,cEnd=0 )

		ExCDSstart=dataConvert$idEx[dataConvert$gStart<=gCDSstart & dataConvert$gEnd>=gCDSstart]
		ExCDSend=dataConvert$idEx[dataConvert$gStart<=gCDSend & dataConvert$gEnd>=gCDSend ]

		if(length(ExCDSstart)==0){
			ExCDSstart = 1
		}
		if(length(ExCDSend)==0){
			ExCDSstart = length(tailleExon )
		}
		dataConvert$cStart[dataConvert$idEx==ExCDSstart]=dataConvert$gStart[dataConvert$idEx==ExCDSstart]-gCDSstart -1
		dataConvert$cEnd[dataConvert$idEx==ExCDSstart]=dataConvert$gEnd[dataConvert$idEx==ExCDSstart]-gCDSstart

		if(ExCDSstart < max(dataConvert$idEx)){
			for (i in seq(from=ExCDSstart+1,to=max(dataConvert$idEx),by=1)){
				dataConvert$cStart[dataConvert$idEx==i]=dataConvert$cEnd[dataConvert$idEx==(i-1)]+1
				dataConvert$cEnd[dataConvert$idEx==i]=dataConvert$cStart[dataConvert$idEx==i]+(dataConvert$lenEx[dataConvert$idEx==i]-1)
			}
		}
		if(ExCDSstart > 1){
			for (i in seq(from=ExCDSstart-1,to=1,by=-1)){
				dataConvert$cEnd[dataConvert$idEx==i]=dataConvert$cStart[dataConvert$idEx==(i+1)]-1
				dataConvert$cStart[dataConvert$idEx==i]=dataConvert$cEnd[dataConvert$idEx==i]-(dataConvert$lenEx[dataConvert$idEx==i]-1)
			}
		}
		cStop=dataConvert$cStart[dataConvert$idEx==ExCDSend]+(gCDSend -dataConvert$gStart[dataConvert$idEx==ExCDSend]-1)

		if(length(grep("-",substr(posVar,2,nchar(posVar)),fixed = T))>0){
			if(length(grep("-",substr(posVar,1,1),fixed = T))>0){
				posVarSplit=unlist(strsplit(posVar,"-",fixed = T))
				posVarSplit=c(paste("-",posVarSplit[2],sep=""),posVarSplit[3])
			}else{
				posVarSplit=unlist(strsplit(posVar,"-",fixed = T))
			}
				posVar1=posVarSplit[1]
				posVar2=as.numeric(posVarSplit[2])
			if(length(grep("*",posVar1,fixed = T))>0){
					posVar1=as.numeric(substr(posVar1,2,nchar(posVar1)))+cStop
			}else{
				posVar1=as.numeric(posVar1)
			}
			gVar <<- dataConvert$gStart[dataConvert$cStart==posVar1]-posVar2

		}else if(length(grep("+",substr(posVar,2,nchar(posVar)),fixed = T))>0){
			posVarSplit=unlist(strsplit(posVar,"+",fixed = T))
			posVar1=posVarSplit[1]
			posVar2=as.numeric(posVarSplit[2])
			if(length(grep("*",posVar1,fixed = T))>0){
				posVar1=as.numeric(substr(posVar1,2,nchar(posVar1)))+cStop
			}else{
				posVar1=as.numeric(posVar1)
			}
			gVar <<- dataConvert$gEnd[dataConvert$cEnd==posVar1]+posVar2

		}else if(length(grep("-",substr(posVar,1,1),fixed = T))>0){
			posVar = as.numeric(posVar)
			if(posVar <=min(dataConvert$cStart)){
				gVar <<- min(dataConvert$gStart)-(abs(posVar)-abs(min(dataConvert$cStart)))
			}else{
				if(abs(posVar)==abs(dataConvert$cEnd[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar])){
					gVar <<- dataConvert$gEnd[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar]
				}else{
					gVar <<- dataConvert$gStart[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar]+
					abs(abs(posVar)-abs(dataConvert$cStart[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar]))
				}
			}
		}else if(length(grep("*",substr(posVar,1,1),fixed = T))>0){
			posVar = as.numeric(substr(posVar,2,nchar(posVar)))
			posVar = posVar + cStop
			if(posVar >=max(dataConvert$cEnd)){
				gVar <<- max(dataConvert$gEnd) + (posVar - max(dataConvert$cEnd))
			}else{
				if(abs(posVar)==abs(dataConvert$cEnd[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar])){
					gVar <<- dataConvert$gEnd[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar]
				}else{
					gVar <<- dataConvert$gStart[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar]+
						(posVar-dataConvert$cStart[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar])
				}
			}
		}else{
			posVar = as.numeric(posVar)
			if(abs(posVar)==abs(dataConvert$cEnd[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar])){
				gVar <<- dataConvert$gEnd[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar]
			}else{
				if(dataConvert$cStart[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar]<0){
					gVar <<- dataConvert$gStart[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar]+
						(posVar-dataConvert$cStart[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar]-1)
				}else{
					gVar <<- dataConvert$gStart[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar]+
						(posVar-dataConvert$cStart[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar])
				}
			}
		}

	}else if(sens=="-"){

		gCDSstart = dataRefSeq[dataRefSeq$V4==transcrit,8]
		gCDSend = dataRefSeq[dataRefSeq$V4==transcrit,7]
		posDon = posStart+tailleCum + 1
		posAcc = posDon+tailleExon - 1

		dataConvert=data.frame(idEx=c(length(tailleExon):1),lenEx=tailleExon,gStart=posAcc,gEnd=posDon,cStart=0,cEnd=0 )

		ExCDSstart=dataConvert$idEx[dataConvert$gStart>=gCDSstart & dataConvert$gEnd<=gCDSstart]
		ExCDSend=dataConvert$idEx[dataConvert$gStart>=gCDSend & dataConvert$gEnd<=gCDSend ]
 		if(length(ExCDSstart)==0){
			ExCDSstart = 1
		}
		if(length(ExCDSend)==0){
			ExCDSstart = length(tailleExon )
		}

		dataConvert$cStart[dataConvert$idEx==ExCDSstart]=gCDSstart-dataConvert$gStart[dataConvert$idEx==ExCDSstart]
		dataConvert$cEnd[dataConvert$idEx==ExCDSstart]=gCDSstart-dataConvert$gEnd[dataConvert$idEx==ExCDSstart]+1

		if(ExCDSstart < max(dataConvert$idEx)){
			for (i in seq(from=ExCDSstart+1,to=max(dataConvert$idEx),by=1)){
				dataConvert$cStart[dataConvert$idEx==i]=dataConvert$cEnd[dataConvert$idEx==(i-1)]+1
				dataConvert$cEnd[dataConvert$idEx==i]=dataConvert$cStart[dataConvert$idEx==i]+(dataConvert$lenEx[dataConvert$idEx==i]-1)
			}
		}
		if(ExCDSstart > 1){
			for (i in seq(from=ExCDSstart-1,to=1,by=-1)){
				dataConvert$cEnd[dataConvert$idEx==i]=dataConvert$cStart[dataConvert$idEx==(i+1)]-1
				dataConvert$cStart[dataConvert$idEx==i]=dataConvert$cEnd[dataConvert$idEx==i]-(dataConvert$lenEx[dataConvert$idEx==i]-1)
			}
		}

		cStop=dataConvert$cStart[dataConvert$idEx==ExCDSend]+(dataConvert$gStart[dataConvert$idEx==ExCDSend] -gCDSend-1)

		if(length(grep("-",substr(posVar,2,nchar(posVar)),fixed = T))>0){
			if(length(grep("-",substr(posVar,1,1),fixed = T))>0){
				posVarSplit=unlist(strsplit(posVar,"-",fixed = T))
				posVarSplit=c(paste("-",posVarSplit[2],sep=""),posVarSplit[3])
			}else{
				posVarSplit=unlist(strsplit(posVar,"-",fixed = T))
			}
			posVar1=posVarSplit[1]
			posVar2=as.numeric(posVarSplit[2])
			if(length(grep("*",posVar1,fixed = T))>0){
				posVar1=as.numeric(substr(posVar1,2,nchar(posVar1)))+cStop
			}else{
				posVar1=as.numeric(posVar1)
			}
			gVar <<- dataConvert$gStart[dataConvert$cStart==posVar1]+(posVar2)

		}else if(length(grep("+",substr(posVar,2,nchar(posVar)),fixed = T))>0){
			posVarSplit=unlist(strsplit(posVar,"+",fixed = T))
			posVar1=posVarSplit[1]
			posVar2=as.numeric(posVarSplit[2])
			if(length(grep("*",posVar1,fixed = T))>0){
				posVar1=as.numeric(substr(posVar1,2,nchar(posVar1)))+cStop
			}else{
				posVar1=as.numeric(posVar1)
			}
				gVar <<- dataConvert$gEnd[dataConvert$cEnd==posVar1]-posVar2

		}else if(length(grep("-",substr(posVar,1,1),fixed = T))>0){
			posVar = as.numeric(posVar)
			if(posVar <=min(dataConvert$cStart)){
				gVar <<- max(dataConvert$gStart)+(abs(posVar)-abs(min(dataConvert$cStart)))
			}else{
				if(abs(posVar)==abs(dataConvert$cEnd[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar])){
					gVar <<- dataConvert$gEnd[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar]
				}else{
					gVar <<- dataConvert$gStart[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar]-
					abs(abs(posVar)-abs(dataConvert$cStart[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar]))
				}
			}
		}else if(length(grep("*",substr(posVar,1,1),fixed = T))>0){
			posVar = as.numeric(substr(posVar,2,nchar(posVar)))
			posVar = posVar + cStop
			if(posVar >=min(dataConvert$cEnd)){
				gVar <<- min(dataConvert$gEnd) - (posVar - max(dataConvert$cEnd))
			}else{
				if(abs(posVar)==abs(dataConvert$cEnd[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar])){
					gVar <<- dataConvert$gEnd[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar]
				}else{
					gVar <<- dataConvert$gStart[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar]-
						(posVar-dataConvert$cStart[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar])
				}
			}
		}else{
			posVar = as.numeric(posVar)
			if(abs(posVar)==abs(dataConvert$cEnd[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar])){
				gVar <<- dataConvert$gEnd[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar]
			}else{
				if(dataConvert$cStart[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar]<0){
					gVar <<- dataConvert$gStart[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar]-
						(posVar-dataConvert$cStart[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar]-1)
				}else{
					gVar <<- dataConvert$gStart[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar]-
						(posVar-dataConvert$cStart[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar])
				}
			}
		}
	}
}

getRevSeq <- function(sequence){
    splitRevSeq = rev(unlist(strsplit(sequence,'|')))
    seqRev = paste(unlist(lapply(list(splitRevSeq),function(x) inverseDic[x,1])),sep="",collapse="")
    return(seqRev)
}

getSeqFromSamtools <- function(chromosome,start,end,strand,command,fastaFile){
	posToAsk <- paste(chromosome,':',start,'-',end,sep="")
	path2script <- paste('faidx',fastaFile,posToAsk)
    samtoolsSeq = system2(command, args=path2script, stdout = TRUE, stderr=TRUE)
    samtoolsHead = samtoolsSeq[1]
	samtoolsSeq = samtoolsSeq[-1]
    if(substr(samtoolsHead,1,1)==">"){
        if(length(samtoolsSeq)>1){
            seqDNApool=""
            for(i in 1:length(samtoolsSeq)){
                seqDNApool=paste(seqDNApool,samtoolsSeq[i],sep="")
            }
            seqDNA = seqDNApool
        }else{
            seqDNA = samtoolsSeq
        }

        seqDNA = toupper(seqDNA)
        if(strand=="-"){
            seqDNA = getRevSeq(seqDNA)
        }
    }else{
        print(samtoolsSeq)
        seqDNA="ERROR samtools"
    }
	return(seqDNA)
}

getSequencePhysio <- function(genome,sens,chr,start,end){
	command = samPath
    chrCorr = chrIndex[chr]
	gPos=c(start,end)
	gPos = gPos[order(gPos)]
    seqDNA = getSeqFromSamtools(chrCorr, gPos[1], gPos[2], sens, command, fastaFile)
    if(seqDNA=="ERROR samtools"){
        message(paste("Request sequence caused a error:", chrCorr, gPos[1], gPos[2], sens, command, fastaFile,"\n"))
    }
    return(seqDNA)
}

checkCalculableScore <- function(seq){
    scoreCalculable<<-TRUE
    if(as.numeric(regexpr("N",seq))>0){
        print(paste("Unknown sequence:",seq))
        scoreCalculable<<-FALSE
    }
}

getSequenceMutated <- function(varPos, sens, seqPhysio, ntChange, varType, genome, chr){

	if(varType=="substitution"){
		ntMut <<- as.character(unlist(strsplit(ntChange,">")))[2]
		seqMutated <<- paste(substr(seqPhysio,1,150),ntMut,substr(seqPhysio,152,301),sep="")
	}else if (varType=="del"){
 		nbDel <<- abs(varPos[1]-varPos[2])+1
		if(sens=="+"){
			seqMutated1 = substr(seqPhysio,1,150)
			start = varPos[2]+1
			end = varPos[2]+151
			seqMutated2=getSequencePhysio (genome,sens,chr,start,end)
			seqMutated <<- paste(seqMutated1 , seqMutated2,sep="")
		}else if(sens=="-"){
			seqMutated1 = substr(seqPhysio,1,150)
			start=varPos[2]-1
			end=varPos[2]-151
			seqMutated2=getSequencePhysio (genome,sens,chr,end,start)
			seqMutated <<- paste(seqMutated1 , seqMutated2,sep="")
		}
	}else if (varType=="ins"){

		ntIns <<- gsub("ins", "", ntChange)
		nbIns <<- nchar(ntIns)
		seqMutated <<- paste(substr(seqPhysio,1,151),ntIns,substr(seqPhysio,152,301-nchar(ntIns)),sep="")

	}else if (varType=="dup"){

		nbIns <<- abs(varPos[1]-varPos[2])+1
		ntIns = substr (seqPhysio,151, 150+nbIns)
		seqMutated1 = substr(seqPhysio,1, 150)
		seqMutated2 = substr(seqPhysio,151+nbIns, 301)
		seqMutated <<- paste(seqMutated1, ntIns, ntIns, seqMutated2,sep="")

	}else if (varType=="delins"){

			nbDel <<- abs(varPos[1]-varPos[2])+1
			ntIns <<- gsub("delins", "", ntChange)
			nbIns <<- nchar(ntIns)

		if(sens=="+"){

			seqMutated1 = substr(seqPhysio,1,150)
			start = varPos[2]+1
			end = varPos[2]+151
			seqMutated2 = getSequencePhysio (genome,sens,chr,start,end)
			seqMutated <<- substr(paste(seqMutated1,ntIns ,seqMutated2,sep=""),1,301)

		}else if(sens=="-"){

			seqMutated1 = substr(seqPhysio,1,150)
			start = varPos[2]-1
			end = varPos[2]-151
			seqMutated2 = getSequencePhysio (genome,sens,chr,end,start)
			seqMutated <<- substr(paste(seqMutated1,ntIns ,seqMutated2,sep=""),1,301)

		}
	}
}

hashseq <- function( seq ){
	seqnum=gsub("A",0,seq)
	seqnum=gsub("C",1,seqnum)
	seqnum=gsub("G",2,seqnum)
	seqnum=gsub("T",3,seqnum)
	lenseq=length(seqnum)
	four=c(1,4,16,64,256,1024,4096,16384)
	i=1
	sum=1
	while(i<=lenseq){
		sum = sum + (as.numeric(seqnum[i])*four[lenseq-i+1])
		i=i+1
	}
		return(sum)
}

MESdonor <- function(seqDon){
    spliseq=unlist(strsplit(seqDon,split=""))
    valueCons1=cons1$score[cons1$N==spliseq[4]]
    valueCons2=cons2$score[cons2$N==spliseq[5]]
    valueN1=bgd$score[bgd$N==spliseq[4]]
    valueN2=bgd$score[bgd$N==spliseq[5]]
    valueCan=(valueCons1*valueCons2)/(valueN1*valueN2)
    seqNcan=paste(spliseq[1],spliseq[2],spliseq[3],spliseq[6],spliseq[7],spliseq[8],spliseq[9],sep="")
    valueNcan=me2x5[seqNcan]
    MESscore <<- log(valueNcan*valueCan)/log(2)
}

MESacceptor <- function(seqAcc){
    spliseqAcc=unlist(strsplit(seqAcc,split=""))
    valueCons1_acc=cons1_acc$score[cons1_acc$N==spliseqAcc[19]]
    valueCons2_acc=cons2_acc$score[cons2_acc$N==spliseqAcc[20]]
    valueN1_acc=bgd$score[bgd$N==spliseqAcc[19]]
    valueN2_acc=bgd$score[bgd$N==spliseqAcc[20]]
    valueCanAcc=(valueCons1_acc*valueCons2_acc)/(valueN1_acc*valueN2_acc)
    seqNcanAcc=c(spliseqAcc[1:18],spliseqAcc[21:23])
    sc=rep(0,9)
    sc[1]=me2x3acc1[hashseq(seqNcanAcc[1:7]),1]
    sc[2]=me2x3acc2[hashseq(seqNcanAcc[8:14]),1]
    sc[3]=me2x3acc3[hashseq(seqNcanAcc[15:21]),1]
    sc[4]=me2x3acc4[hashseq(seqNcanAcc[5:11]),1]
    sc[5]=me2x3acc5[hashseq(seqNcanAcc[12:18]),1]
    sc[6]=me2x3acc6[hashseq(seqNcanAcc[5:7]),1]
    sc[7]=me2x3acc7[hashseq(seqNcanAcc[8:11]),1]
    sc[8]=me2x3acc8[hashseq(seqNcanAcc[12:14]),1]
    sc[9]=me2x3acc9[hashseq(seqNcanAcc[15:18]),1]
    scoreNcan=(sc[1] * sc[2] * sc[3] * sc[4] * sc[5]) /
    	(sc[6] * sc[7] * sc[8] * sc[9])
    MESscoreAcc <<- log(scoreNcan*valueCanAcc)/log(2)
}

getMES <- function(sstype,seqCon){
    checkCalculableScore(seqCon)
    if(scoreCalculable){
    	if(sstype=="Acc"){
    		scoreMES = MESacceptor(seqCon)
    	}else{
    		if(substr(as.vector(seqCon),4,5)=="GT"){
    			scoreMES = MESdonor(seqCon)
    		}else{
    			scoreMES = 0
    		}
    	}
    }else{scoreMES = 0}
	return(scoreMES)
}

SSFdonGT <- function(SeqDonGT){
	for(j in 1:9){
		n=substr(as.vector(SeqDonGT),j,j)
		i_score[j] = ref_score_GT[n,j]
	}
	SSFdonScoreGT <<- ((sum(i_score)-mint_GT)/(maxt_GT-mint_GT))*100
}

SSFdonGC <- function(SeqDonGC){
	for(j in 1:9){
		n=substr(as.vector(SeqDonGC),j,j)
		i_score[j] = ref_score_GC[n,j]
	}
	SSFdonScoreGC <<- ((sum(i_score)-mint_GC)/(maxt_GC-mint_GC))*100
}

SSFacc <- function(SeqAcc1,SeqAcc2){
	if(regexpr(pattern="AG",SeqAcc1,fixed=T)>=1){
		SSFaccScore <<- 0
	}else{
		for(l in 1:10){
			n=substr(SeqAcc1,l,l)
			i_score1[l]=ref_score_AG[n,l]
		}

		for(m in 1:4){
			n=substr(SeqAcc2,m ,m )
			i_score2[m]=ref_score_AG[n,(m+11) ]
		}
		score1=i_score1[order(i_score1,decreasing=T)]
		score1=score1[1:8]
		SSFaccScore <<- ((((sum(score1)-mint1)/(maxt1-mint1)) + ((sum(i_score2)-mint2)/(maxt2-mint2)))/2)*100
	}
}

getSSF <- function(sstype,seqCon){
    checkCalculableScore(seqCon)
    if(scoreCalculable){
    	if (sstype=="Acc"){
    		seqSSFadj=substr(seqCon,7,21)
    		nEch1=substr(as.vector(seqSSFadj),1,10)
    		nEch2=substr(as.vector(seqSSFadj),12,15)
    		SSF = SSFacc(nEch1,nEch2)

    	}else{
    		if(substr(as.vector(seqCon),4,5)=="GC"){
    			SSF = SSFdonGT(seqCon)
    		}else{
    			SSF = SSFdonGT(seqCon)
    		}
    	}
    }else{SSF = 0}
	return(SSF)
}

getESRscore <- function(sstype,seq){
	ESRscore=NULL
	if(sstype=="Acc"){
		for (i in indAcc){
			motif = substr(seq,i,i+5)
			if(motif%in%ESRmotif & motif%in%ESRmotif){
				if(i<12){
					ESRscore = c(ESRscore,LEIsc_valuesWA[motif])
				}else if (i<63){
					ESRscore = c(ESRscore,LEIsc_valuesHA[motif])
				}else {
					ESRscore = c(ESRscore,LEIsc_valuesHM[motif])
				}
			}
		}
	}else{
		for (i in indDon){
			motif = substr(seq,i,i+5)
			if(motif%in%ESRmotif & motif%in%ESRmotif){
				if(i<29){
					ESRscore = c(ESRscore,LEIsc_valuesHM[motif])
				}else if (i<75){
					ESRscore = c(ESRscore,LEIsc_valuesHD[motif])
				}else {
					ESRscore = c(ESRscore,LEIsc_valuesWD[motif])
				}
			}
		}
	}
	if(is.null(ESRscore)){
		ESRscoreFinal = 0
	}else{
		ESRscoreFinal = mean(ESRscore)
	}
    return(ESRscoreFinal)
}

getScore <- function(sstype,seq,seqCon){
    MES = getMES(sstype,seqCon)
	SSF = getSSF(sstype,seqCon)
	ESR = getESRscore(sstype,seq)
	result <<-list(MES,SSF,ESR)
}

getProbaModel <- function(sstype, MES, SSF, ESR){
	if(sstype=="Don"){
		proba = exp(-1.170e+01 + SSF*(5.652e-02 - 1.770e-02) + MES*(4.547e-01 + 6.000e-02) + ESR*(7.005e+00 + 1.239e-01) + 1.756e+00)/
			(1+exp(-1.170e+01 + SSF*(5.652e-02 - 1.770e-02) + MES*(4.547e-01 + 6.000e-02) + ESR*(7.005e+00 + 1.239e-01) + 1.756e+00))
	}else{
		proba = exp(-1.170e+01 + SSF*(5.652e-02) + MES*(4.547e-01) + ESR*(7.005e+00))/
			(1+exp(-1.170e+01 + SSF*(5.652e-02) + MES*(4.547e-01) + ESR*(7.005e+00)))
	}
	return(proba)
}

getSeqToStudyWTAcc <- function(chr, sens, posAccPhy){
	if(sens=="+"){
		start=posAccPhy-19
		end=posAccPhy+100
		seqPhysioAcc <- getSequencePhysio (genome,sens,chr,start,end)
	}else if(sens=="-"){
		start=posAccPhy+20
		end=posAccPhy-99
		seqPhysioAcc <- getSequencePhysio (genome,sens,chr,end,start)
	}
	seqConsAccPhyNew <<- substr(seqPhysioAcc,1,23)
	seqExonAccPhyNew <<- substr(seqPhysioAcc,21,nchar(seqPhysioAcc))
}

getSeqToStudyWTDon <- function(chr, sens, posDonPhy){
	if(sens=="+"){
		start=posDonPhy-99
		end=posDonPhy+6
		seqPhysioDon <<- getSequencePhysio (genome,sens,chr,start,end)
	}else if(sens=="-"){
		start=posDonPhy+100
		end=posDonPhy-5
		seqPhysioDon <<- getSequencePhysio (genome,sens,chr,end,start)
	}
		seqConsDonPhyNew <<- substr(seqPhysioDon,98,nchar(seqPhysioDon))
		seqExonDonPhyNew <<- substr(seqPhysioDon,1,100)
}

getSplitTableSeq <- function(varName, chr, varPos, sens, seqPhysio, seqMutated, posDon, posAcc, nearestPosAll, varType){

	relPosAccPhy = as.numeric(gregexpr("AG",seqPhysio)[[1]])
	relPosDonPhy = as.numeric(gregexpr("GT",seqPhysio)[[1]])
	relPosAccMut = as.numeric(gregexpr("AG",seqMutated)[[1]])
	relPosDonMut = as.numeric(gregexpr("GT",seqMutated)[[1]])

	relPosAccPhyFilt = relPosAccPhy[relPosAccPhy >= 115 & relPosAccPhy <= 185]
	relPosDonPhyFilt = relPosDonPhy[relPosDonPhy >= 130 & relPosDonPhy <= 170]
	relPosAccMutFilt = relPosAccMut[relPosAccMut >= 115 & relPosAccMut <= 185]
	relPosDonMutFilt = relPosDonMut[relPosDonMut >= 130 & relPosDonMut <= 170]

	seqConsAccPhy = substr(rep(seqPhysio,length(relPosAccPhyFilt)),relPosAccPhyFilt-18,relPosAccPhyFilt+4)
	seqConsDonPhy = substr(rep(seqPhysio,length(relPosDonPhyFilt)),relPosDonPhyFilt-3,relPosDonPhyFilt+5)
	seqConsAccMut = substr(rep(seqMutated,length(relPosAccMutFilt)),relPosAccMutFilt-18,relPosAccMutFilt+4)
	seqConsDonMut = substr(rep(seqMutated,length(relPosDonMutFilt)),relPosDonMutFilt-3,relPosDonMutFilt+5)

	seqExonAccPhy = substr(rep(seqPhysio,length(relPosAccPhyFilt)),relPosAccPhyFilt+2,relPosAccPhyFilt+101)
	seqExonDonPhy = substr(rep(seqPhysio,length(relPosDonPhyFilt)),relPosDonPhyFilt-100,relPosDonPhyFilt-1)
	seqExonAccMut = substr(rep(seqMutated,length(relPosAccMutFilt)),relPosAccMutFilt+2,relPosAccMutFilt+101)
	seqExonDonMut = substr(rep(seqMutated,length(relPosDonMutFilt)),relPosDonMutFilt-100,relPosDonMutFilt-1)

	if(sens == "+"){
		PosAccPhy = varPos + (relPosAccPhyFilt-150)
		PosDonPhy = varPos + (relPosDonPhyFilt-152)
		PosAccMut = varPos + (relPosAccMutFilt-150)
		PosDonMut = varPos + (relPosDonMutFilt-152)
		if(varType!="substitution"){
			if(varType=="del"){
				PosAccMut[PosAccMut>varPos] = PosAccMut[PosAccMut>varPos] + nbDel
				PosDonMut[PosDonMut>varPos] = PosDonMut[PosDonMut>varPos] + nbDel
			}else if(varType=="ins" | varType=="dup"){
				PosAccMut[PosAccMut>varPos] = PosAccMut[PosAccMut>varPos] - nbIns
				PosDonMut[PosDonMut>varPos] = PosDonMut[PosDonMut>varPos] - nbIns
			}else if(varType=="delins"){
				PosAccMut[PosAccMut>varPos] = PosAccMut[PosAccMut>varPos] - (nbIns - nbDel)
				PosDonMut[PosDonMut>varPos] = PosDonMut[PosDonMut>varPos] - (nbIns - nbDel)
			}
		}
	}else{
		PosAccPhy = varPos - (relPosAccPhyFilt-149)
		PosDonPhy = varPos - (relPosDonPhyFilt-151)
		PosAccMut = varPos - (relPosAccMutFilt-149)
		PosDonMut = varPos - (relPosDonMutFilt-151)
		if(varType!="substitution"){
			if(varType=="del"){
				PosAccMut[PosAccMut<varPos] = PosAccMut[PosAccMut<varPos] - nbDel
				PosDonMut[PosDonMut<varPos] = PosDonMut[PosDonMut<varPos] - nbDel
			}else if(varType=="ins" | varType=="dup"){
				PosAccMut[PosAccMut<varPos] = PosAccMut[PosAccMut<varPos] + nbIns
				PosDonMut[PosDonMut<varPos] = PosDonMut[PosDonMut<varPos] + nbIns
			}else if(varType=="delins"){
				PosAccMut[PosAccMut<varPos] = PosAccMut[PosAccMut<varPos] + (nbIns - nbDel)
				PosDonMut[PosDonMut<varPos] = PosDonMut[PosDonMut<varPos] + (nbIns - nbDel)
			}
		}
	}
	tmpTableSeq = data.frame(var = rep(varName,length(c(PosAccPhy,PosDonPhy,PosAccMut,PosDonMut))),
					chr = rep(chr,length(c(PosAccPhy,PosDonPhy,PosAccMut,PosDonMut))),
					relPos = c(relPosAccPhyFilt,relPosDonPhyFilt,relPosAccMutFilt,relPosDonMutFilt),
					pos = c(PosAccPhy,PosDonPhy,PosAccMut,PosDonMut),
					seqCons = c(seqConsAccPhy,seqConsDonPhy,seqConsAccMut,seqConsDonMut),
					seqExon = c(seqExonAccPhy,seqExonDonPhy,seqExonAccMut,seqExonDonMut),
					seqType = c(rep("WT",length(c(PosAccPhy,PosDonPhy))),rep("Mut",length(c(PosAccMut,PosDonMut)))),
					sstype = c(rep("Acc",length(PosAccPhy)),rep("Don",length(PosDonPhy)),rep("Acc",length(PosAccMut)),rep("Don",length(PosDonMut)))
					)
    if(nrow(tmpTableSeq)>0){
    	tmpTableSeq$seqCons = as.character(tmpTableSeq$seqCons)
    	tmpTableSeq$seqExon = as.character(tmpTableSeq$seqExon)
    	tmpTableSeq$Physio = "No"
    	tmpTableSeq$Physio[which(tmpTableSeq$pos%in%posAcc & tmpTableSeq$sstype=="Acc")] = "Yes"
    	tmpTableSeq$Physio[which(tmpTableSeq$pos%in%posDon & tmpTableSeq$sstype=="Don")] = "Yes"
    	if(nrow(tmpTableSeq[tmpTableSeq$Physio=="Yes" & tmpTableSeq$seqType=="WT",])==1){
    		SSphy = tmpTableSeq$sstype[tmpTableSeq$Physio=="Yes" & tmpTableSeq$seqType=="WT"]
    		tmpTableSeq = tmpTableSeq[tmpTableSeq$sstype==SSphy,]
    	}

    	if(length(nearestPosAll)==2){
    		if(length(which(tmpTableSeq$sstype=="Don"))>0){
    			if(length(which(tmpTableSeq$pos==nearestPosAll[1] & tmpTableSeq$sstype=="Don"))==0){
    				getSeqToStudyWTDon(chr, sens, posDonPhy = nearestPosAll[1])
    				tmpTableSeq = rbind(tmpTableSeq,c(varName,chr,NA,nearestPosAll[1],as.character(seqConsDonPhyNew),as.character(seqExonDonPhyNew),"WT","Don","Yes"))
    			}
    		}
    		if(length(which(tmpTableSeq$sstype=="Acc"))>0){
    			if(length(which(tmpTableSeq$pos==nearestPosAll[2] & tmpTableSeq$sstype=="Acc"))==0){
    				getSeqToStudyWTAcc(chr, sens, posAccPhy = nearestPosAll[2])
    				tmpTableSeq = rbind(tmpTableSeq,c(varName[1],chr[1],NA,nearestPosAll[2],as.character(seqConsAccPhyNew),as.character(seqExonAccPhyNew),"WT","Acc","Yes"))
    			}
    		}
    	}
    	if(nrow(tmpTableSeq[tmpTableSeq$Physio=="Yes",])!=0){
    		tmpTableSeq = tmpTableSeq[!(tmpTableSeq$pos%in%tmpTableSeq$pos[tmpTableSeq$Physio=="Yes"] & tmpTableSeq$Physio!="Yes"),]
    	}
        tmpTableSeq = tmpTableSeq[!is.na(tmpTableSeq$sstype),]
    	tmpScore = mapply(getScore, tmpTableSeq$sstype, as.character(tmpTableSeq$seqExon), as.character(tmpTableSeq$seqCons))

    	tmpTableSeq$MES = unlist(tmpScore[1,])
    	tmpTableSeq$SSF = unlist(tmpScore[2,])
    	tmpTableSeq$ESR = unlist(tmpScore[3,])

    	tmpProba = mapply(getProbaModel, tmpTableSeq$sstype, tmpTableSeq$MES, tmpTableSeq$SSF , tmpTableSeq$ESR)
    	tmpTableSeq$proba = as.numeric(unlist(tmpProba))
    	tmpTableSeq$proba[tmpTableSeq$proba<0]=0
    	tmpTableSeq$proba[tmpTableSeq$proba>1]=1

    	tmpTableSeq$classProba = "No"
    	tmpTableSeq$classProba[tmpTableSeq$proba >= 0.00388264838136624] = "Yes"
    	tmp = merge(tmpTableSeq[tmpTableSeq$seqType=="WT",c("sstype","pos","proba","Physio")],tmpTableSeq[tmpTableSeq$seqType=="Mut",c("sstype","pos","proba","Physio")],
    					by.x = c("sstype","pos"),by.y = c("sstype","pos"),all.x=F,all.y=T)
    	tmp$proba.x[is.na(tmp$proba.x)]=0
    	tmp$Physio.x[is.na(tmp$Physio.x)] = "No"

    	tmp = tmp[order(tmp$proba.y,decreasing=T),]
    	tmp = tmp[!duplicated(tmp$pos),]
    	if(length(which(tmpTableSeq$seqType=="Mut" & tmpTableSeq$pos%in%tmp$pos[tmp$proba.x>=tmp$proba.y & tmp$Physio.x!="Yes"]))>0){
    		tmpTableSeq = tmpTableSeq[-which(tmpTableSeq$seqType=="Mut" & tmpTableSeq$pos%in%tmp$pos[tmp$proba.x>=tmp$proba.y & tmp$Physio.x!="Yes"]),]
    	}
    }
	return(tmpTableSeq)
}

getDeltaESRseq <- function(SstypePhy, distSS, seqPhysio, seqMutated){
	if(abs(distSS)>120){
		ESRscore <<- NA
	}else{
		seqESRwt = substr(seqPhysio,146,156)
		seqESRmut = substr(seqMutated,146,156)
		ESRscoreMut=NULL
		ESRscoreWT=NULL
		j=1
		if(SstypePhy=="acceptor"){
			rangeToCheck = (distSS-5):(distSS)
			if(min(rangeToCheck)<=0){
				rangeToCheck = 1:distSS
				seqESRwt = substr(seqESRwt,7-distSS,11)
				seqESRmut = substr(seqESRmut,7-distSS,11)
			}
			for (i in rangeToCheck){
				motifWT = substr(seqESRwt,j,j+5)
				motifMUT = substr(seqESRmut,j,j+5)
				j = j +1
				if(motifWT%in%ESRmotif){
					ESRscoreWT = c(ESRscoreWT,ESRlistScore[motifWT])
				}
				if(motifMUT%in%ESRmotif){
					ESRscoreMut = c(ESRscoreMut,ESRlistScore[motifMUT])
				}
			}
		}else{
			rangeToCheck = (distSS):(distSS+5)
			if(max(rangeToCheck)>=0){
				rangeToCheck = distSS:(-1)
				seqESRwt = substr(seqESRwt,1,11-(6+distSS))
				seqESRmut = substr(seqESRmut,1,11-(6+distSS))
			}
			for (i in rangeToCheck){
				motifWT = substr(seqESRwt,j,j+5)
				motifMUT = substr(seqESRmut,j,j+5)
				j = j +1
				if(motifWT%in%ESRmotif){
					ESRscoreWT = c(ESRscoreWT,ESRlistScore[motifWT])
				}
				if(motifMUT%in%ESRmotif){
					ESRscoreMut = c(ESRscoreMut,ESRlistScore[motifMUT])
				}
			}
		}
		ESRscore = sum(ESRscoreMut) - sum(ESRscoreWT)
	}
	return(ESRscore)
}

getSeqCons <- function(SstypePhy, distSS, seq, varType, varPos, ntChange){
	if (varType=="del"){
		nbDel <- abs(varPos[1]-varPos[2])+1
		if(SstypePhy=="acceptor"){
			if(distSS<0){
				if(nbDel>abs(distSS)){
					seqCons = substr(seq, 131, 153)
				}else{
					seqCons = substr(seq, 150-(19+distSS+nbDel), 150+(3-(distSS+nbDel)))
				}
			}else{
				seqCons = substr(seq, 151-(19+distSS), 151+(3-distSS))
			}
		}else{
			if(distSS<0){
				if(nbDel>abs(distSS)){
					seqCons = substr(seq, 148, 156)
				}else{
					seqCons = substr(seq, 150-(2+distSS+nbDel), 150+(6-(distSS+nbDel)))
				}
			}else{
				seqCons = substr(seq, 151-(2+distSS), 151+(6-distSS))
			}
		}
	}else if (varType=="dup" | varType=="ins"){
		if (varType=="ins"){
			nbIns <- nchar(gsub("ins", "", ntChange))
		}else{
			nbIns <- abs(varPos[1]-varPos[2])+1
		}
		if(SstypePhy=="acceptor"){
			if(distSS<0){
				seqCons = substr(seq, 150-(19+distSS) + nbIns, 150+(3-distSS) + nbIns)
			}else{
				seqCons = substr(seq, 151-(19+distSS), 151+(3-distSS))
			}
		}else{
			if(distSS<0){
				seqCons = substr(seq, 150-(2+distSS) + nbIns, 150+(6-distSS) + nbIns)
			}else{
				seqCons = substr(seq, 151-(2+distSS), 151+(6-distSS))
			}
		}
	}else if (varType=="delins"){
		nbDel <- abs(varPos[1]-varPos[2])+1
		nbIns <- nchar(gsub("delins", "", ntChange))
		if(SstypePhy=="acceptor"){
			if(distSS<0){
				if(nbDel>abs(distSS)){
					seqCons = substr(seq, 131, 153)
				}else{
					seqCons = substr(seq, 150-(19+distSS) + (nbIns-nbDel), 150+(3-distSS) + (nbIns-nbDel))
				}
			}else{
				seqCons = substr(seq, 151-(19+distSS), 151+(3-distSS))
			}
		}else{
			if(distSS<0){
				if(nbDel>abs(distSS)){
					seqCons = substr(seq, 148, 156)
				}else{
					seqCons = substr(seq, 150-(2+distSS) + (nbIns-nbDel), 150+(6-distSS) + (nbIns-nbDel))
				}
			}else{
				seqCons = substr(seq, 151-(2+distSS), 151+(6-distSS))
			}
		}
	}else{
		if(SstypePhy=="acceptor"){
			if(distSS<0){
				seqCons = substr(seq, 150-(19+distSS), 150+(3-distSS))
			}else{
				seqCons = substr(seq, 151-(19+distSS), 151+(3-distSS))
			}
		}else{
			if(distSS<0){
				seqCons = substr(seq, 150-(2+distSS), 150+(6-distSS))
			}else{
				seqCons = substr(seq, 151-(2+distSS), 151+(6-distSS))
			}
		}
	}
	return(seqCons)
}

getSPiCE <- function(SstypePhy, seqConsWT, seqConsMut){
	SSFwt = getSSF(paste(toupper(substr(SstypePhy,1,1)),substr(SstypePhy,2,3),sep=""),seqConsWT)
	SSFmut = getSSF(paste(toupper(substr(SstypePhy,1,1)),substr(SstypePhy,2,3),sep=""),seqConsMut)
	MESwt = getMES(paste(toupper(substr(SstypePhy,1,1)),substr(SstypePhy,2,3),sep=""),seqConsWT)
	MESmut = getMES(paste(toupper(substr(SstypePhy,1,1)),substr(SstypePhy,2,3),sep=""),seqConsMut)

	if (MESwt==0){
		deltaMES=0
	}else{
		deltaMES = (MESmut-MESwt)/MESwt
	}
	if (SSFwt==0){
		deltaSSF=0
	}else{
		deltaSSF = (SSFmut-SSFwt)/SSFwt
	}
    if(exp(-3.59-8.21*round(deltaMES,3)-32.30*round(deltaSSF,3))==Inf){
        SPiCEproba <<- 1
    }else{
        SPiCEproba <<- round(exp(-3.59-8.21*round(deltaMES,3)-32.30*round(deltaSSF,3))/(1+exp(-3.59-8.21*round(deltaMES,3)-32.30*round(deltaSSF,3))),5)
    }

	if(SPiCEproba<0.115){
		SPiCEinter_2thr <<- "low"
	}else if(SPiCEproba>0.749){
		SPiCEinter_2thr <<- "high"
	}else{
		SPiCEinter_2thr <<- "medium"
	}
}

getDeltaMES <- function(SstypePhy, seqConsWT, seqConsMut){
	MESwt = getMES(paste(toupper(substr(SstypePhy,1,1)),substr(SstypePhy,2,3),sep=""),seqConsWT)
	MESmut = getMES(paste(toupper(substr(SstypePhy,1,1)),substr(SstypePhy,2,3),sep=""),seqConsMut)
	if (MESwt==0){
		deltaMES=0
	}else{
		deltaMES = (MESmut-MESwt)/MESwt
	}
	return(deltaMES)
}

getBPParea <- function(varPos,transcript,genome){
	if(genome=="hg19"){
		dataBPannot = dataBPannot19
	}else{
		dataBPannot = dataBPannot38
	}

	tmpAnnotBP = dataBPannot[dataBPannot$transcrit==transcript,]
	if(nrow(tmpAnnotBP)==0){
		mutInPBareaBPP = "No BPP-predicted BP in intron"
	}else{
		if(length(varPos)==1){
			if(nrow(tmpAnnotBP[tmpAnnotBP$start<=varPos & tmpAnnotBP$end>=varPos,])==1){
				posBP = tmpAnnotBP[tmpAnnotBP$start<=varPos & tmpAnnotBP$end>=varPos,"posBP"]
				distToAcc = tmpAnnotBP[tmpAnnotBP$start<=varPos & tmpAnnotBP$end>=varPos,"distToAcc"]
				score = tmpAnnotBP[tmpAnnotBP$start<=varPos & tmpAnnotBP$end>=varPos,"score"]
				mutInPBareaBPP = paste("Yes g.",posBP," (",distToAcc,"): ",round(score,3), sep = "")
			}else{
				mutInPBareaBPP = "No"
			}
		}else if(length(varPos)==2) {
			if(nrow(tmpAnnotBP[(tmpAnnotBP$start<=varPos[1] & tmpAnnotBP$end>=varPos[1]) |
								(tmpAnnotBP$start<=varPos[2] & tmpAnnotBP$end>=varPos[2]) |
								(tmpAnnotBP$start>varPos[1] & tmpAnnotBP$end<varPos[2]),])==1){
				posBP = tmpAnnotBP[(tmpAnnotBP$start<=varPos[1] & tmpAnnotBP$end>=varPos[1]) |
								(tmpAnnotBP$start<=varPos[2] & tmpAnnotBP$end>=varPos[2]) |
								(tmpAnnotBP$start>varPos[1] & tmpAnnotBP$end<varPos[2]),"posBP"]
				distToAcc = tmpAnnotBP[(tmpAnnotBP$start<=varPos[1] & tmpAnnotBP$end>=varPos[1]) |
								(tmpAnnotBP$start<=varPos[2] & tmpAnnotBP$end>=varPos[2]) |
								(tmpAnnotBP$start>varPos[1] & tmpAnnotBP$end<varPos[2]),"distToAcc"]
				score = tmpAnnotBP[(tmpAnnotBP$start<=varPos[1] & tmpAnnotBP$end>=varPos[1]) |
								(tmpAnnotBP$start<=varPos[2] & tmpAnnotBP$end>=varPos[2]) |
								(tmpAnnotBP$start>varPos[1] & tmpAnnotBP$end<varPos[2]),"score"]
				mutInPBareaBPP = paste("Yes g.",posBP," (",distToAcc,"): ",round(score,3), sep = "")
			}else{
				mutInPBareaBPP = "No"
			}
		}
	}
	return(mutInPBareaBPP)
}

getMutInfo <- function(mutInput){
	if(length(grep("delins",mutInput ))>0){
		varPos = unlist(strsplit(mutInput,"delins"))[1]
		ntChange = paste("delins",unlist(strsplit(mutInput,"delins"))[2],sep="")
	}else if (length(grep("del",mutInput ))>0){
		varPos = unlist(strsplit(mutInput,"del"))[1]
		ntChange = "del"
	}else if(length(grep("ins",mutInput ))>0){
		varPos = unlist(strsplit(mutInput,"ins"))[1]
		ntChange = paste("ins",unlist(strsplit(mutInput,"ins"))[2],sep="")
	}else if(length(grep("dup",mutInput ))>0){
		varPos = unlist(strsplit(mutInput,"dup"))[1]
		ntChange = "dup"
	}else if(length(grep(">", mutInput, fixed = TRUE))>0){
		varPos = substr(mutInput,1,nchar(mutInput)-3)
		ntChange = substr(mutInput,nchar(mutInput)-2,nchar(mutInput))
	}else{
		print("Error of mut annotation")
	}
	varPos <<- varPos
	ntChange <<- ntChange
}

getVariantInfo <- function(varID){
    varID = as.character(varID)
	varDecomp=unlist(strsplit(varID,":"))
	if(length(varDecomp)!=3 & length(varDecomp)!=2){
		message("You must import variant as:Transcrit:position(:)nucleotidic change")
	}else{
		transcript = varDecomp[1]
		if(length(grep(".",transcript ))>0){
			transcript = unlist(strsplit(transcript,".",fixed = T))[1]
		}
		if(length(varDecomp)==2){
			getMutInfo(varDecomp[2])
		}else{
			varPos = as.character(varDecomp[2])
			ntChange = as.character(varDecomp[3])
		}
		if(length(grep("delins",ntChange ))>0){
			varType = "delins"
			if(length(grep("_",varPos ))>0){
				varPos = unlist(strsplit(varPos ,"_"))
				if(length(grep("c.",varPos[1]))>0){
					varPos[1] = substr(varPos[1],3,nchar(varPos[1]))
					convertcNomenIngNomen(transcript, varPos[1])
					varPos[1] = gVar
					convertcNomenIngNomen(transcript, varPos[2])
					varPos[2] = gVar
				}else{
					varPos[1]=substr(varPos[1],3,nchar(varPos[1]))
				}
			}else{
				if(length(grep("c.",varPos[1]))>0){
					varPos[1] = substr(varPos[1],3,nchar(varPos[1]))
					convertcNomenIngNomen(transcript, varPos[1])
					varPos = c(gVar,gVar)
				}else{
					varPos[1] = substr(varPos[1],3,nchar(varPos[1]))
					varPos = c(varPos,varPos)
				}
			}
		}else if (length(grep("del",ntChange ))>0){
			varType = "del"
			if(length(grep("_",varPos ))>0){
				varPos = unlist(strsplit(varPos ,"_"))
				if(length(grep("c.",varPos[1]))>0){
					varPos[1] = substr(varPos[1],3,nchar(varPos[1]))
					convertcNomenIngNomen(transcript, varPos[1])
					varPos[1] = gVar
					convertcNomenIngNomen(transcript, varPos[2])
					varPos[2] = gVar
				}else{
					varPos[1] = substr(varPos[1],3,nchar(varPos[1]))
				}
			}else{
				if(length(grep("c.",varPos[1]))>0){
					varPos[1] = substr(varPos[1],3,nchar(varPos[1]))
					convertcNomenIngNomen(transcript, varPos[1])
					varPos = c(gVar,gVar)
				}else{
					varPos[1] = substr(varPos[1],3,nchar(varPos[1]))
					varPos = c(varPos,varPos)
				}
			}
		}else if(length(grep("ins",ntChange ))>0){
			varType = "ins"
			if(length(grep("c.",varPos))>0){
				varPos = substr(varPos,3,nchar(varPos))
				convertcNomenIngNomen(transcript, varPos)
				varPos = gVar
			}else{
				varPos=substr(varPos,3,nchar(varPos))
			}
		}else if(length(grep("dup",ntChange ))>0){
			varType = "dup"
			if(length(grep("_",varPos ))>0){
				varPos = unlist(strsplit(varPos ,"_"))
				if(length(grep("c.",varPos[1]))>0){
					varPos[1] = substr(varPos[1],3,nchar(varPos[1]))
					convertcNomenIngNomen(transcript, varPos[1])
					varPos[1] = gVar
					convertcNomenIngNomen(transcript, varPos[2])
					varPos[2] = gVar
				}else{
					varPos[1] = substr(varPos[1],3,nchar(varPos[1]))
				}
			}else{
				if(length(grep("c.",varPos[1]))>0){
					varPos[1] = substr(varPos[1],3,nchar(varPos[1]))
					convertcNomenIngNomen(transcript, varPos[1])
					varPos = c(gVar,gVar)
				}else{
					varPos[1] = substr(varPos[1],3,nchar(varPos[1]))
					varPos = c(varPos,varPos)
				}
			}
		}else{
			varType = "substitution"
			if(length(grep("c.",varPos))>0){
				varPos=substr(varPos,3,nchar(varPos))
				convertcNomenIngNomen(transcript, varPos)
				varPos = gVar
			}else{
				varPos=substr(varPos,3,nchar(varPos))
			}
		}
	}
	varID <<- varID
	varPos <<- as.numeric(varPos)
	varType <<- varType
	transcript <<- transcript
	ntChange <<- ntChange
}

getGlobaInterpretation <- function(SPiCEinterpret, RegType, deltaMES, mutInPBareaBPP, classProbaMut, classProbaWT, distSS, deltaESR){
	interpretFinal = NULL
	if(SPiCEinterpret=="high" | SPiCEinterpret=="medium"){
		interpretFinal = c(interpretFinal,"Alter by SPiCE")
	}
	if (length(grep("PolyTC",RegType ))>0 & deltaMES<(-0.15)){
		interpretFinal = c(interpretFinal,"Alter by MES (Poly TC)")
	}
	if (length(grep("BP",RegType ))>0 & mutInPBareaBPP!="No"){
        if(mutInPBareaBPP=="No BPP-predicted BP in intron"){
            interpretFinal = c(interpretFinal,"Pos BP unknown")
        }else{
            interpretFinal = c(interpretFinal,"Alter BP")
        }
	}
    if (length(grep("Exon",RegType ))>0 & abs(distSS)<120 & deltaESR<(-1.10)){
		interpretFinal = c(interpretFinal,"Alter ESR")
	}
	if (classProbaMut=="Yes"){
        if(classProbaWT=="Yes"){
            if(abs(distSS)>=150 & length(grep("Intron",RegType ))>0) {
                interpretFinal =c(interpretFinal, "Alter by create Exon")
            }else{
                interpretFinal = c(interpretFinal,"Alter by create Cryptic")
            }
        }else{
            if(abs(distSS)>=150 & length(grep("Intron",RegType ))>0) {
    			interpretFinal =c(interpretFinal, "Alter by create de Novo Exon")
    		}else{
    			interpretFinal = c(interpretFinal,"Alter by create de Novo Cryptic")
    		}
        }
	}
	if(is.null(interpretFinal)){
		interpretFinal = "NTR"
	}else if (length(interpretFinal)>1){
		interpretFinal = paste(interpretFinal, collapse=" + ")
	}
	return(interpretFinal)
}

getPredConfident <- function(interpretFinal, RegType, distSS, SstypePhy){
	#proba from SNP + UV analysis (N = 66,000)
	probaInter = -1
	if (length(grep("Exon",RegType ))>0){
        if(SstypePhy=="acceptor"){
		    if(length(grep("Alter by SPiCE",interpretFinal ))>0){
                probaInter = "47.37 % [24.92 % ; 69.82 %]"
            }else if (length(grep("Alter ESR",interpretFinal ))>0){
                probaInter = "27.33 % [22.54 % ; 32.12 %]"
            }else if (length(grep("Alter by create Cryptic",interpretFinal ))>0){
                if(abs(distSS[1])<=120){
                    probaInter = "06.25 % [02.69 % ; 13.81 %]"
                }else{
                    probaInter = "04.17 % [02.13 % ; 08.01 %]"
                }
            }else if (length(grep("de Novo",interpretFinal ))>0){
                if(abs(distSS[1])<=120){
                    probaInter = "23.81 % [05.59 % ; 42.03 %]"
                }else{
                    probaInter = "03.39 % [00.94 % ; 11.54 %]"
                }
            }else{
                if (length(grep("Cons",RegType ))>0){
                    probaInter = "09.09 % [02.53 % ; 27.81 %]"
                }else if (abs(distSS[1])<=120){
                    probaInter = "06.33 % [04.59 % ; 08.67 %]"
                }else{
                    probaInter = "01.09 % [00.55 % ; 02.13 %]"
                }
            }
        }else{
            if(length(grep("Alter by SPiCE",interpretFinal ))>0){
                probaInter = "85.71 % [76.90 % ; 91.82 %]"
            }else if (length(grep("Alter ESR",interpretFinal ))>0){
                probaInter = "34.05 % [28.49 % ; 39.61 %]"
            }else if (length(grep("Alter by create Cryptic",interpretFinal ))>0){
                if(abs(distSS[1])<=120){
                    probaInter = "09.88 % [05.09 % ; 18.29 %]"
                }else{
                    probaInter = "04.17 % [02.13 % ; 08.01 %]"
                }
            }else if (length(grep("de Novo",interpretFinal ))>0){
                if(abs(distSS[1])<=120){
                    probaInter = "27.27 % [08.66% ; 45.88 %]"
                }else{
                    probaInter = "03.39 % [00.94 % ; 11.54 %]"
                }
            }else{
                if (length(grep("Cons",RegType ))>0){
                    probaInter = "07.14 % [01.98 % ; 22.64 %]"
                }else if (abs(distSS[1])<=120){
                    probaInter = "05.35 % [03.76 % ; 07.58 %]"
                }else{
                    probaInter = "01.09 % [00.55 % ; 02.13 %]"
                }
            }
        }
    }else{
        if(SstypePhy=="acceptor"){
            if(length(grep("Alter by SPiCE",interpretFinal ))>0){
                if(distSS[1]>=(-2)){
                    probaInter = "96.15 % [92.60 % ; 98.04 %]"
                }else{
                    probaInter = "69.57 % [61.89 % ; 77.25 %]"
                }
            }else if (length(grep("Alter by MES",interpretFinal ))>0){
                probaInter = "83.33 % [72.79 % ; 93.87 %]"
            }else if (length(grep("Alter BP",interpretFinal ))>0){
                probaInter = "48.57 % [36.86 % ; 60.28 %]"
            }else if(length(grep("Alter by create Cryptic",interpretFinal ))>0){
                if(length(grep("Cons",RegType ))>0 | length(grep("PolyTC",RegType))>0 | length(grep("BP",RegType))>0){
                    probaInter = "No available"
                }else{
                    probaInter = "02.84 % [01.11 % ; 07.07 %]"
                }
            }else if(length(grep("Alter by create Exon",interpretFinal ))>0){
                probaInter = "00.50 % [00.31 % ; 00.83 %]"
            }else if (length(grep("de Novo Cryptic",interpretFinal ))>0){
                if(length(grep("Cons",RegType ))>0){
                    probaInter = "No available"
                }else if (length(grep("PolyTC",RegType))>0){
                    probaInter = "63.64 % [35.21 % ; 92.07 %]"
                }else if (length(grep("BP",RegType))>0){
                    probaInter = "18.18 % [00.00 % ; 40.97 %]"
                }else{
                    probaInter = "09.62 % [04.17 % ; 20.61 %]"
                }
            }else if (length(grep("de Novo Exon",interpretFinal ))>0){
                probaInter = "02.26 % [01.45 % ; 03.51%]"
    		}else{
                if(length(grep("Cons",RegType ))>0){
                    probaInter = "00.93 % [00.17 % ; 05.11 %]"
                }else if (length(grep("PolyTC",RegType))>0){
                    probaInter = "06.90 % [03.20 % ; 14.24 %]"
                }else if (length(grep("BP",RegType))>0){
                    probaInter = "01.95 % [00.76 % ; 04.90 %]"
                }else if (abs(distSS[1])<=150){
                    probaInter = "00.13 % [00.02 % ; 00.74 %]"
                }else{
                    probaInter = "00.06 % [00.03 % ; 00.11 %]"
                }
            }
    	}else{
    		if(length(grep("Alter by SPiCE",interpretFinal ))>0){
                if(distSS[1]<=2){
                    probaInter = "98.42 % [96.00 % ; 99.38 %]"
                }else{
                    probaInter = "95.35 % [90.23 % ; 97.85 %]"
                }
    		}else if(length(grep("Alter by create Cryptic",interpretFinal ))>0){
                if(length(grep("Cons",RegType ))>0){
                    probaInter = "No available"
                }else{
                    probaInter = "03.63 % [01.76 % ; 07.29 %]"
                }
    		}else if(length(grep("Alter by create Exon",interpretFinal ))>0){
    			probaInter = "00.30 % [00.16 % ; 00.56 %]"
            }else if (length(grep("de Novo Cryptic",interpretFinal ))>0){
                if(length(grep("Cons",RegType ))>0){
                    probaInter = "No available"
                }else{
                    probaInter = "16.67 % [08.06 % ; 25.28 %]"
                }
            }else if (length(grep("de Novo Exon",interpretFinal ))>0){
    			probaInter = "03.31 % [02.29 % ; 04.77 %]"
    		}else{
        		if(length(grep("Cons",RegType ))>0){
                    probaInter = "00.00 % [00.00 % ; 25.88 %]"
                }else if(abs(distSS[1])<=150){
                    probaInter = "00.25 % [00.08 % ; 00.74 %]"
                }else{
                    probaInter = "00.04 % [00.02 % ; 00.08%]"
                }
            }
        }
	}
	return (probaInter)
}

getOutput <- function(){
	getPosSSphysio(transcript)
	getNearestPos(sens, varPos ,posDon, posAcc)

	start=varPos[1]-150
	end=varPos[1]+150
	seqPhysio <<- getSequencePhysio (genome,sens,chr,start,end)
	getSequenceMutated (varPos, sens, seqPhysio, ntChange, varType, genome, chr)
	seqPhysio <<- toupper(as.character(seqPhysio))
	seqMutated <<- toupper(as.character(seqMutated))
	tmpTableSeq <<- getSplitTableSeq(varID, chr, varPos[1], sens, seqPhysio, seqMutated, posDon, posAcc, nearestPosAll,varType)

	if(as.numeric(gregexpr("Exon",RegType))<0){
		ESRscore = NA
	}else{
		ESRscore = getDeltaESRseq(SstypePhy, distSS, seqPhysio, seqMutated)
	}

	if(as.numeric(gregexpr("Cons",RegType))<0){
		seqConsWT = ""
		seqConsMut = ""
		SPiCEproba = 0
		SPiCEinter_2thr = "Outside SPiCE Interpretation"
	}else{
		seqConsWT = getSeqCons(SstypePhy, distSS, seqPhysio, varType = "WT", varPos, ntChange)
		seqConsMut = getSeqCons(SstypePhy, distSS, seqMutated, varType, varPos, ntChange)
		getSPiCE(SstypePhy, seqConsWT, seqConsMut)
	}
	if(as.numeric(gregexpr("PolyTC",RegType))<0){
		seqConsWT = ""
		seqConsMut = ""
		deltaMES = 0
	}else{
		seqConsWT = getSeqCons(SstypePhy, distSS, seqPhysio, varType = "WT", varPos, ntChange)
		seqConsMut = getSeqCons(SstypePhy, distSS, seqMutated, varType, varPos, ntChange)
		deltaMES = getDeltaMES(SstypePhy, seqConsWT, seqConsMut)
	}
	if(as.numeric(gregexpr("BP",RegType))<0){
		mutInPBareaBPP = "No"
	}else{
		mutInPBareaBPP = getBPParea(varPos,transcript,genome)
	}

	tmpTableSeqNoPhyMut = tmpTableSeq[tmpTableSeq$Physio!="Yes" & tmpTableSeq$seqType == "Mut",]

	strand <- sens
	gNomen <- if(length(varPos)==1){varPos}else{paste(varPos,collapse=";")}
	NearestSS <- SstypePhy
	if(length(varPos)==1){
		DistSS <- distSS
	}else if(length(varPos)==2){
		if(distSS[1]==distSS2){
			DistSS <- distSS[1]
		}else{
			DistSS <- paste(distSS, distSS2,sep="; ")
		}
	}else{
		print("erreur varpos")
	}
	gene <- as.character(dataRefSeq$V13[dataRefSeq$V4==transcript])
	mutInPBarea <- mutInPBareaBPP
	deltaESRscore <- ESRscore

	if(nrow(tmpTableSeqNoPhyMut)>0){
		tmpTableSeqMAXMut = tmpTableSeqNoPhyMut[tmpTableSeqNoPhyMut$proba==max(tmpTableSeqNoPhyMut$proba),]
		sstypCrypt = tmpTableSeqMAXMut$sstype[1]
		tmpTableSeqNoPhyWT = tmpTableSeq[tmpTableSeq$Physio!="Yes" & tmpTableSeq$seqType == "WT"& tmpTableSeq$sstype==sstypCrypt,]
		tmpTableSeqMAXWT = tmpTableSeqNoPhyWT[tmpTableSeqNoPhyWT$proba==max(tmpTableSeqNoPhyWT$proba),]

		tmpTableSeqPhy = tmpTableSeq[tmpTableSeq$Physio=="Yes"& tmpTableSeq$seqType == "WT" & tmpTableSeq$sstype==sstypCrypt,]
		tmpTableSeqPhyMut = tmpTableSeq[tmpTableSeq$Physio=="Yes"& tmpTableSeq$seqType == "Mut" & tmpTableSeq$sstype==sstypCrypt,]

		getNearestPosCrypt(sens, as.numeric(tmpTableSeqMAXMut$pos) ,posDon, posAcc)
		posCryptMut <- tmpTableSeqMAXMut$pos
		sstypeCryptMut <- as.character(tmpTableSeqMAXMut$sstype)
		nearestSStoCrypt <- SstypePhyCrypt
		nearestPosSStoCrypt <- nearestPosPhyCrypt
		nearestDistSStoCrypt <- distSScrypt
		probaCryptMut <- tmpTableSeqMAXMut$proba
		classProbaCryptMut <- tmpTableSeqMAXMut$classProba
		if(nrow(tmpTableSeqMAXWT)>0){
			posCryptWT <- tmpTableSeqMAXWT$pos[1]
			probaCryptWT <- tmpTableSeqMAXWT$proba[1]
			classProbaCryptWT <- tmpTableSeqMAXWT$classProba[1]
		} else{
            posCryptWT <- 0
            probaCryptWT <- 0
            classProbaCryptWT <- "No"
        }
		if(nrow(tmpTableSeqPhy)>0){
			posSSPhysio <- tmpTableSeqPhy$pos[1]
			probaSSPhysio <- tmpTableSeqPhy$proba[1]
			classProbaSSPhysio <- tmpTableSeqPhy$classProba[1]
			probaSSPhysioMut <- tmpTableSeqPhy$proba[1]
			classProbaSSPhysioMut <- tmpTableSeqPhy$classProba[1]
		}else{
            posSSPhysio <- 0
            probaSSPhysio <- 0
            classProbaSSPhysio <- "No"
            probaSSPhysioMut <- 0
            classProbaSSPhysioMut <- "No"
        }
		if(nrow(tmpTableSeqPhyMut)>0){
			probaSSPhysioMut <- tmpTableSeqPhyMut$proba[1]
			classProbaSSPhysioMut <- tmpTableSeqPhyMut$classProba[1]
		}
	}else{
		posCryptMut <- 0
		sstypeCryptMut <- "No site"
		nearestSStoCrypt <- "No site"
		nearestPosSStoCrypt <- 0
		nearestDistSStoCrypt <- 0
		probaCryptMut <- 0
		classProbaCryptMut <- "No"
		posCryptWT <- 0
		probaCryptWT <- 0
		classProbaCryptWT <- "No"
		posSSPhysio <- 0
		probaSSPhysio <- 0
		classProbaSSPhysio <- "No"
		probaSSPhysioMut <- 0
		classProbaSSPhysioMut <- "No"
	}
	interpretation <- getGlobaInterpretation(SPiCEinter_2thr, RegType, deltaMES, mutInPBareaBPP, classProbaCryptMut, classProbaCryptWT, distSS[1], ESRscore)
	Interpretation <- interpretation
	InterConfident <- getPredConfident(interpretation, RegType, distSS, SstypePhy)
    ExonInfo <- getExonInfo(transcript,varPos[1])
    result <<- c(Interpretation, InterConfident, chr, strand, gNomen, varType, ntChange, ExonInfo, transcript, gene, NearestSS, DistSS, RegType, seqPhysio, seqMutated,
        SPiCEproba, SPiCEinter_2thr, deltaMES, mutInPBarea, deltaESRscore, posCryptMut, sstypeCryptMut, probaCryptMut,
        classProbaCryptMut, nearestSStoCrypt, nearestPosSStoCrypt, nearestDistSStoCrypt, posCryptWT, probaCryptWT,
        classProbaCryptWT, posSSPhysio, probaSSPhysio, classProbaSSPhysio, probaSSPhysioMut, classProbaSSPhysioMut)
}

splitRawToTable <- function(raw, sep = "\t", head = TRUE){
    if(head){
        columNames = unlist(strsplit(raw[1],sep,fixed=TRUE))
        nCol = length(columNames)
        splitRaw = unlist(strsplit(raw[2:length(raw)],sep,fixed=TRUE))
        data=as.data.frame(matrix(splitRaw, ncol = nCol,byrow = TRUE))
        colnames(data) <- columNames
    }else{
        if(as.numeric(unlist(gregexpr("\t", raw[1] ,fixed=TRUE)))<0){
            nCol=1
        }else{
            nCol = length(as.numeric(unlist(gregexpr("\t", raw[1] ,fixed=TRUE))))+1
        }
        splitRaw = unlist(strsplit(raw, sep, fixed=TRUE))
        data=as.data.frame(matrix(splitRaw, ncol = nCol, byrow = TRUE))
    }
    return(data)
}

convertLine2VCF <- function(varID,SPiPout=NULL,error=FALSE){
    ID = as.character(varID)
    QUAL = "."
    FILTER = "."
    if(!error){
        names(SPiPout) = c("Interpretation", "InterConfident", "chr", "strand", "gNomen", "varType", "ntChange", "ExonInfo", "transcript",
            "gene", "NearestSS", "DistSS", "RegType", "seqPhysio", "seqMutated",
            "SPiCEproba", "SPiCEinter_2thr", "deltaMES", "mutInPBarea", "deltaESRscore", "posCryptMut", "sstypeCryptMut", "probaCryptMut",
            "classProbaCryptMut", "nearestSStoCrypt", "nearestPosSStoCrypt", "nearestDistSStoCrypt", "posCryptWT", "probaCryptWT",
            "classProbaCryptWT", "posSSPhysio", "probaSSPhysio", "classProbaSSPhysio", "probaSSPhysioMut", "classProbaSSPhysioMut")
        seq = SPiPout["seqPhysio"]
        strand = SPiPout["strand"]
        varType = SPiPout["varType"]
        if(strand=="-"){seq = getRevSeq(seq)}

        CHROM = SPiPout["chr"]
        if(varType=="substitution"){
            POS = SPiPout["gNomen"]
            REF = if(strand=="+"){unlist(strsplit(ntChange,">",fixed=T))[1]}else if(strand=="-"){getRevSeq(unlist(strsplit(ntChange,">",fixed=T))[1])}
            ALT = if(strand=="+"){unlist(strsplit(ntChange,">",fixed=T))[2]}else if(strand=="-"){getRevSeq(unlist(strsplit(ntChange,">",fixed=T))[2])}
        }else if(varType=="ins"){
            POS = as.numeric(SPiPout["gNomen"])-1
            REF = substr(seq,151,151)
            ALT = paste0(REF,if(strand=="+"){gsub("ins","",ntChange)}else if(strand=="-"){getRevSeq(gsub("ins","",ntChange))})
        }else if(varType=="dup"){
            varPOS = as.numeric(unlist(strsplit(SPiPout["gNomen"],";",fixed=T)))
            if(length(varPOS)==1){varPOS=rep(varPOS,2)}
            POS = min(varPOS)-1
            mutSize = abs(varPOS[1]-varPOS[2])
            ntDup = substr(seq,151,151+mutSize)
            REF = paste0(substr(seq,150,150),ntDup)
            ALT = paste0(substr(seq,150,150),paste(rep(ntDup,2),collapse=""))
        }else if(varType=="del"){
            varPOS = as.numeric(unlist(strsplit(SPiPout["gNomen"],";",fixed=T)))
            if(length(varPOS)==1){varPOS=rep(varPOS,2)}
            POS = min(varPOS)-1
            mutSize = abs(varPOS[1]-varPOS[2])
            ntDel = substr(seq,151,151+mutSize)
            REF = paste0(substr(seq,150,150),ntDel)
            ALT = substr(seq,150,150)
        }else if(varType=="delins"){
            varPOS = as.numeric(unlist(strsplit(SPiPout["gNomen"],";",fixed=T)))
            if(length(varPOS)==1){varPOS=rep(varPOS,2)}
            POS = min(varPOS)-1
            mutSize = abs(varPOS[1]-varPOS[2])
            ntDel = substr(seq,151,151+mutSize)
            ntIns = if(strand=="+"){gsub("delins","",ntChange)}else if(strand=="-"){getRevSeq(gsub("delins","",ntChange))}
            REF = paste0(substr(seq,150,150),ntDel)
            ALT = paste0(substr(seq,150,150),ntIns)
        }

        INFO = paste(c("Interpretation", "InterConfident", "strand", "varType", "ntChange", "ExonInfo", "transcript", "gene", "NearestSS",
                    "DistSS", "RegType", "SPiCEproba", "SPiCEinter_2thr", "deltaMES", "mutInPBarea", "deltaESRscore", "posCryptMut", "sstypeCryptMut",
                    "probaCryptMut", "classProbaCryptMut", "nearestSStoCrypt", "nearestPosSStoCrypt", "nearestDistSStoCrypt", "posCryptWT",
                    "probaCryptWT", "classProbaCryptWT", "posSSPhysio", "probaSSPhysio", "classProbaSSPhysio", "probaSSPhysioMut", "classProbaSSPhysioMut"),
                    paste0("\"",SPiPout[c("Interpretation", "InterConfident", "strand", "varType", "ntChange", "ExonInfo", "transcript", "gene", "NearestSS",
                                "DistSS", "RegType", "SPiCEproba", "SPiCEinter_2thr", "deltaMES", "mutInPBarea", "deltaESRscore", "posCryptMut", "sstypeCryptMut",
                                "probaCryptMut", "classProbaCryptMut", "nearestSStoCrypt", "nearestPosSStoCrypt", "nearestDistSStoCrypt", "posCryptWT",
                                "probaCryptWT", "classProbaCryptWT", "posSSPhysio", "probaSSPhysio", "classProbaSSPhysio", "probaSSPhysioMut", "classProbaSSPhysioMut")],
                            "\""),
                    sep="=",collapse="|")
    }else{
        CHROM <- POS <- REF <- ALT <- INFO  <- "."
    }

    VCFline = c(CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO)
    return(VCFline)
}

SPiP <- function(varID,i){
    setTxtProgressBar(pb2, i)
    if(as.numeric(regexpr('no transcript',varID))>0|
        as.numeric(regexpr('mutUnknown',varID))>0)
    {
        if(printVCF){
            return(paste(convertLine2VCF(varID,error=T),collapse="\t"))
        }else{
            return(paste(rep("NA",35),collapse="\t"))
        }
    }else{
        tryCatch({
            getVariantInfo(as.character(varID))
            tmp <- getOutput()
            if(printVCF){
                return(paste(convertLine2VCF(varID,tmp,error=F),collapse="\t"))
            }else{
                return(paste(tmp,collapse="\t"))
            }
        },
        error=function(cond) {
            message(paste("Variant caused a error:", varID))
            if(printVCF){
                return(paste(convertLine2VCF(varID,error=T),collapse="\t"))
            }else{
                return(paste(rep("NA",35),collapse="\t"))
            }
        })
    }
}
#launch analysis

T1 <- as.numeric(format(Sys.time(), "%s"))

#import data
readVCF <- function(dataLine,i){
    setTxtProgressBar(pb, i)
    dataLine = unlist(strsplit(dataLine,split='\t')) #c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
    chrom <- dataLine[1]
    pos <- as.numeric(dataLine[2])
    ref <- dataLine[4]
    alt <- dataLine[5]
    if(ref=="."|alt=="."){
        variant = "mutUnknown"
    }else{
        if(substr(as.character(chrom[1]),1,3)!="chr"){
            if(substr(as.character(chrom[1]),1,3)=="NC_"){
                chrom = contigToChr(as.character(chrom))
            }else{
                chrom = paste("chr",chrom,sep="")
            }
        }
        transcript = as.character(dataRefSeq[dataRefSeq$V1==chrom & dataRefSeq$V2<=pos & dataRefSeq$V3>=pos,'V4'])
        strandTrans = as.character(dataRefSeq[dataRefSeq$V1==chrom & dataRefSeq$V2<=pos & dataRefSeq$V3>=pos,'V6'])
        if(length(transcript)==0){
            variant = paste("no transcript",pos,sep=":")
        }else{
            altSplit = unlist(strsplit(alt,',',fixed = TRUE))
            variant = NULL
            for(i in 1:length(transcript)){
                if(strandTrans[i]=="+"){
                    for(j in 1:length(altSplit)){
                        if(nchar(ref)==1 & nchar(altSplit[j])==1){
                            variant = c(variant,paste(transcript[i],':g.',pos,':',ref,'>',altSplit[j],sep=""))
                        }else{
                            variant = c(variant,paste(transcript[i],':g.',pos,"_",pos+nchar(ref)-1,':delins',altSplit[j],sep=""))
                        }
                    }
                }else if(strandTrans[i]=="-"){
                    refRev = getRevSeq(ref)
                    for(j in 1:length(altSplit)){
                        altRev = getRevSeq(altSplit[j])
                        if(nchar(refRev)==1 & nchar(altRev)==1){
                            variant = c(variant,paste(transcript[i],':g.',pos,':',refRev,'>',altRev,sep=""))
                        }else{
                            variant = c(variant,paste(transcript[i],':g.',pos+nchar(refRev)-1,"_",pos,':delins',altRev,sep=""))
                        }
                    }
                }
            }
        }
    }
    result <<- list(variant)
}

s<-0 # first iteration
input<-file(inputFile,"r")

message(paste(s+1,s+maxLines,sep=" to "))
data=NULL
rawInput<-readLines(input, n=maxLines)
if(length(rawInput)==0) stop()
rawInput<-rawInput[!nchar(rawInput)<3]

message(paste(sub("CET",":",Sys.time(),fixed=T),"Read lines..."))
if(fileFormat=="txt"){
    data = splitRawToTable(rawInput,head=TRUE)
    columNames <- names(data)
    rawInput = rawInput[-1]
    if(is.null(data$varID)){
        message("###########################")
        message("#Your data doesn't have the varID column")
        message("###########################")
        message(helpMessage)
        stop()
    }
}else if(fileFormat=="vcf"){
    if(as.numeric(regexpr("#",rawInput[length(rawInput)]))<0){
        mHeader = rawInput[grep("#",rawInput,fixed=TRUE)]
        columNames <- c(mHeader[length(mHeader)],"varID")
        mHeader = mHeader[-length(mHeader)]
        message(paste(sub("CET",":",Sys.time(),fixed=T),"Align VCF mutation on transcripts..."))
        rawInput = rawInput[-grep("#",rawInput,fixed=TRUE)]
        total <- length(rawInput)
        pb <- txtProgressBar(min = 0, max = total, initial = 1, char = "=", style = 3)
        tmpVCF = unlist(mcmapply(FUN = readVCF,rawInput,1:length(rawInput), mc.cores = threads, mc.preschedule = TRUE))
        data = data.frame(varID = as.character(tmpVCF))
        rawInput = paste(names(tmpVCF),data[,'varID'],sep="\t")
    }else{mHeader = rawInput}
}

if(!is.null(data)){
    message(paste("\n",gsub("CET",":",Sys.time(),fixed=T),"Score Calculation..."))
    total <- nrow(data)
    pb2 <- txtProgressBar(min = 0, max = total, initial = 1, char = "=", style = 3)
    rawResult<-foreach (i=1:nrow(data),.errorhandling='pass') %dopar% {
        SPiP(data[i,"varID"],i)
    }
#    print(rawResult[nchar(rawResult)<600])
#    rawResult[nchar(rawResult)<600] = paste(rep("NA",35),collapse="\t")
    message(paste("\n",sub("CET",":",Sys.time(),fixed=T),"Write results..."))

    colNames <- paste(c(columNames, "Interpretation", "InterConfident", "chr", "strand", "gNomen", "varType", "ntChange",
        "ExonInfo", "transcript", "gene", "NearestSS", "DistSS", "RegType", "seqPhysio", "seqMutated",
        "SPiCEproba", "SPiCEinter_2thr", "deltaMES", "mutInPBarea", "deltaESRscore", "posCryptMut", "sstypeCryptMut", "probaCryptMut",
        "classProbaCryptMut", "nearestSStoCrypt", "nearestPosSStoCrypt", "nearestDistSStoCrypt", "posCryptWT", "probaCryptWT",
        "classProbaCryptWT", "posSSPhysio", "probaSSPhysio", "classProbaSSPhysio", "probaSSPhysioMut", "classProbaSSPhysioMut"),collapse="\t")

    if(!printHead & !printVCF){output<-file(outputFile,"w")}else{output<-file(outputFile,"a")}
    if(printVCF){
        writeLines(c("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", unlist(rawResult)), con = output,sep="\n")
    }else{
        writeLines(c(if(fileFormat=="vcf"){mHeader},
            colNames,
            paste(rawInput,rawResult,sep="\t")),
            con = output,sep="\n")
    }
    flush(output)
    close(output)
}

s<-s+maxLines

# Remaining iteration
output<-file(outputFile,"a")
while(T){
    data=NULL
    printmHeader = FALSE
    message(paste(s+1,s+maxLines,sep=" to "))
    rawInput<-readLines(input, n=maxLines)
    if(length(rawInput)==0) break
    rawInput<-rawInput[!nchar(rawInput)<3]

    message(paste(sub("CET",":",Sys.time(),fixed=T),"Read lines..."))
    if(fileFormat=="txt"){
        data = splitRawToTable(rawInput,head=FALSE)
        colnames(data) <- columNames
    }else if(fileFormat=="vcf"){
        if(as.numeric(regexpr("#",rawInput[length(rawInput)]))<0){
            if(as.numeric(regexpr("#",rawInput[1]))>0){
                printmHeader = TRUE
                mHeader = c(mHeader,rawInput[grep("#",rawInput,fixed=TRUE)])
                columNames <- c(mHeader[length(mHeader)],"varID")
                mHeader = mHeader[-length(mHeader)]
                rawInput = rawInput[-grep("#",rawInput,fixed=TRUE)]
            }
            message(paste(sub("CET",":",Sys.time(),fixed=T),"Align VCF mutation on transcripts..."))
            total <- length(rawInput)
            pb <- txtProgressBar(min = 0, max = total, initial = 1, char = "=", style = 3)
            tmpVCF = unlist(mcmapply(FUN = readVCF,rawInput,1:length(rawInput), mc.cores = threads, mc.preschedule = TRUE))
            data = data.frame(varID = as.character(tmpVCF))
            rawInput = paste(names(tmpVCF),data[,'varID'],sep="\t")
        }else{mHeader = c(mHeader,rawInput);printmHeader = TRUE}
    }
    if(!is.null(data)){
        message(paste("\n",gsub("CET",":",Sys.time(),fixed=T),"Score Calculation..."))
        total <- nrow(data)
        pb2 <- txtProgressBar(min = 0, max = total, initial = 1, char = "=", style = 3)
        rawResult<-foreach (i=1:nrow(data),.errorhandling='pass') %dopar% {
            SPiP(data[i,"varID"],i)
        }
        print(rawResult[nchar(rawResult)<600])
        rawResult[nchar(rawResult)<600] = paste(rep("NA",35),collapse="\t")

        message(paste("\n",sub("CET",":",Sys.time(),fixed=T),"Write results..."))
        if(printVCF){
            writeLines(unlist(rawResult), con = output,sep="\n")
        }else{
            writeLines(c(if(printmHeader){mHeader},paste(rawInput,rawResult,sep="\t")), con = output, sep = "\n")
        }
        flush(output)
    }
    s<-s+maxLines
}

close(input)
close(output)

T2 <- as.numeric(format(Sys.time(), "%s"))

print(paste("Runtime:",round((T2 - T1),3),"sec"))
