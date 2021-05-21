#!/usr/bin/env Rscript

#######################
# SPiP software
#######################
# author Raphael Leman r.leman@baclesse.unicancer.fr, Center François Baclesse and Normandie University, Unicaen, Inserm U1245
# Copyright 2019 Center François Baclesse and Normandie University, Unicaen, Inserm U1245

# This software was developed from the work:
# SPiP, a comprehensive Splicing Prediction Pipeline for massive detection of exonic and intronic variant effect on mRNA splicing.
# Raphaël Leman, Béatrice Parfait, Dominique Vidaud, Emmanuelle Girodon, Laurence Pacot, Gérald Le Gac, Chandran Ka, Claude Ferec, Yann Fichou, Céline Quesnelle,
# Etienne Muller, Dominique Vaur, Laurent Castera, Agathe Ricou, Hélène Tubeuf, Omar Soukarieh, Pascaline Gaildrat, Florence Riant, Marine Guillaud-Bataille,
# Sandrine M. Caputo, Virginie Caux-Moncoutier, Nadia Boutry-Kryza, Françoise Bonnet-Dorion, Ines Schultz, Maria Rossing, Louis Goldenberg, Olivier Quenez,
# Valentin Harter, Michael T. Parsons, Amanda B. Spurdle, Thierry Frébourg, Alexandra Martins, Claude Houdayer, Sophie Krieger

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
# FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#launch analysis
T1 <- as.numeric(format(Sys.time(), "%s"))

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

tryCatch({
library(randomForest)
},
	error=function(cond) {
		message("Here's the original error message:")
		message(cond)
		message("*****You need to install \'randomForest\' library\nInstall it by: install.pakages(\'randomForest\')")
})

cat("
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
threads = 1
genome="hg19"
maxLines = 1000
printVCF = FALSE
printProcess = FALSE
pathToGene = NULL
pathToTranscript = NULL
pathToTranscriptome = NULL
version = "2.0"

#SPiP arguments
helpMessage=paste0("Usage: SPiPv",version,".r\n
    Mandatory \n
        -I, --input /path/to/inputFile\t\tlist of variants file (.txt or .vcf)
        -O, --output /path/to/outputFile\t\tName of ouput file (.txt or .vcf)\n
    Genome options \n
        -g, --GenomeAssenbly hg19\t\tGenome assembly version (hg19 or hg38) [default= ",genome,"] \n
    Parallel options \n
        -t, --threads N\t\tNumber of threads used for the calculation [default= ",threads,"]
        -l, --maxLines N\t\tNumber of lines read in each time [default= ",maxLines,"]\n
    Other options\n
        --geneList /path/to/geneList.txt\t\tlist of gene to study
        --transcriptList /path/to/transcriptList.txt\t\tlist of transcript to study
        --transcriptome /path/to/transcriptome_hgXX.RData\t\tTranscriptome path if file is not in /path/to/SPiP/RefFiles/
        --VCF\t\tPrint output in vcf format
        --header\t\tPrint meta-header info
        --verbose\t\tShow run process
    -h, --help\t\tPrint this help message and exit\n
   You could : Rscript SPiPv",version,".r -I ./testCrypt.txt -O ./outTestCrypt.txt")

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
        inputFile=normalizePath(args[i+1]);i = i+2
    }else if(args[i]=="-O"|args[i]=="--output"){
        outputFile=args[i+1];i = i+2
    }else if(args[i]=="--geneList"){
        pathToGene=normalizePath(args[i+1]);i = i+2
    }else if(args[i]=="--transcriptList"){
        pathToTranscript=normalizePath(args[i+1]);i = i+2
    }else if(args[i]=="--transcriptome"){
        pathToTranscriptome=normalizePath(args[i+1]);i = i+2
    }else if(args[i]=="-g"|args[i]=="--GenomeAssenbly"){
        genome=args[i+1];i = i+2
    }else if(args[i]=="-t"|args[i]=="--threads"){
        threads= as.numeric(args[i+1]);i = i+2
    }else if(args[i]=="-l"|args[i]=="--maxLines"){
        maxLines=as.numeric(args[i+1]);i = i+2
    }else if(args[i]=="--VCF"){
        printVCF=TRUE;i = i+1
    }else if(args[i]=="--verbose"){
        printProcess=TRUE;i = i+1
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

registerDoParallel(threads)
CMD = paste0(normalizePath(sub("--file=","",argsFull[substr(argsFull,1,7)=="--file="])),
        " --input ", inputFile,
        " --output ", outputFile,
        " --GenomeAssenbly ", genome,
        " --threads ", threads,
        " --maxLines ", maxLines,
        if(!is.null(pathToGene)){paste0(" --geneList ",pathToGene)},
        if(!is.null(pathToTranscript)){paste0(" --transcriptList ",pathToTranscript)},
        if(!is.null(pathToTranscriptome)){paste0(" --transcriptome ",pathToTranscriptome)},
        if(printVCF){" --VCF "})

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

cat("##################\n")
cat("#Your command:\n")
cat("##################\n")
cat(paste0(CMD,"\n"))

#Get Ref files
inputref = paste(scriptPath, "/RefFiles",sep="")

path2scripts = paste0(inputref,"/SPiP_libs/")
tryCatch({
    source(paste0(path2scripts,"SPiP_functions.r"))
},
	error=function(cond) {
		message("Here's the original error message:")
		message(cond)
		message("*****Error while importing the source file SPiP_functions.r")
	}
)

# Reads headers from corresponding files, stores them as vectors of strings
path2headers = paste0(inputref,"/headers/") # path from script to folder that contains the headers

headerHelp = readLines(paste(path2headers,"headerHelp.txt",sep=""))
headerHelp[1] = paste0("##SPiP output v",version) # modifies the dynamic line containing the version
headerHelp[2] = paste0("##SPiPCommand=",CMD) # modifies the dynamic line containing the CMD

cat("Check transcriptome sequences...\n")
if(!file.exists(paste(inputref,"/transcriptome_hg19.RData",sep="")) | !file.exists(paste(inputref,"/transcriptome_hg38.RData",sep=""))){
    if(!is.null(pathToTranscriptome)){
        cat(paste("Your transcriptome file:",pathToTranscriptome))
        load(pathToTranscriptome)
    }else{
        cat("You have to install the transcriptome file in /path/to/SPiP/RefFiles/\n")
        cat("transcriptome_hg19.RData available at : https://sourceforge.net/projects/splicing-prediction-pipeline/files/transcriptome/transcriptome_hg19.RData/download\n")
        cat("transcriptome_hg38.RData available at : https://sourceforge.net/projects/splicing-prediction-pipeline/files/transcriptome/transcriptome_hg38.RData/download\n")
        q(save="no")
    }
}
cat("Load transcriptome sequences...\n")
load(paste0(inputref, "/transcriptome_",genome,".RData"))

cat("Load SPiP model...\n")
load(paste0(inputref, "/model.RData"))

cat("Load VPP table...\n")
VPPtable = read.table(paste0(inputref, "/VPP_table.txt"),sep="\t",dec=",",header=TRUE)

cat("Load VPN table...\n")
VPNtable = read.table(paste0(inputref, "/VPN_table.txt"),sep="\t",dec=",",header=TRUE)

cat("Check RefSeq database...\n")
if(!file.exists(paste(inputref,"/dataRefSeqhg19.RData",sep="")) & !file.exists(paste(inputref,"/dataRefSeqhg38.RData",sep=""))){
    currentWD = getwd()
    setwd(scriptPath)
	cat("Create RefSeq database...\n")
	source(paste(inputref,"/getRefSeqDatabase.r",sep=""),local =TRUE)
    setwd(currentWD)
}
cat("Load RefSeq database...\n")
load(paste0(inputref, "/dataRefSeq",genome,".RData"))
load(paste0(inputref, "/RefFiles.RData"))

if(!is.null(pathToGene)){
    geneList = readLines(pathToGene)
    dataRefSeq = dataRefSeq[which(as.character(dataRefSeq$V13)%in%geneList),]
}
if(!is.null(pathToTranscript)){
    transcriptList = readLines(pathToTranscript)
    dataRefSeq = dataRefSeq[which(as.character(dataRefSeq$V4)%in%transcriptList),]
}

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
RegTypeToNumber <- data.frame(V1 = c(1:11),
						row.names = c("DeepIntron", "Exon",
						"ExonESR", "ExonESRCons", "Intron",
						"IntronBP", "IntronCons", "IntronConsPolyTC",
						"IntronConsPolyTCBP", "IntronPolyTC",
						"IntronPolyTCBP"))

thToSPiPexon = 0.18
thToSPiPintron = 0.035

####################
# READING INPUT FILE
####################

# io
input<-file(inputFile,"r")
output<-file(outputFile,"w")

if(fileFormat=="vcf"){
    cat("Reading header...\n")
    mHeader <- NULL
    rawInput <- readLines(input, n=1) # line by line # first line
    firstLine <- 1 # index of current line
    while (as.numeric(regexpr("^#CHROM",rawInput)) != 1 ) {
        # while one has not reached the last line of the header
        mHeader <- c(mHeader,rawInput)
        rawInput <- readLines(input, n=1)
        firstLine <- firstLine + 1 # index of current line
    }
    # last line of the header
    mHeader <- c(mHeader, headerHelp) # add lines for SPiP
    mHeader <- c(mHeader, rawInput) # add line #CHROM ...
    columnsNamesLine <- rawInput # line CHROM POS ID REF ALT QUAL FILTER INFO ...
    columnsNames <- unlist(strsplit(columnsNamesLine,split='\t')) #c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO",SUITE?)
    columnsNames[1] <- substr(columnsNames[1],2,nchar(columnsNames[1])) # delete the # of #CHROM
    lengthDataLine <- length(columnsNames)

    if(printVCF){
        writeLines(mHeader, con = output,sep="\n")
    }else{
        writeLines(paste(paste(columnsNames,collapse="\t"), "varID", "Interpretation", "InterConfident", "SPiPscore", "strand",
            "gNomen", "varType", "ntChange", "ExonInfo", "exonSize", "transcript", "gene", "NearestSS", "DistSS", "RegType", "SPiCEproba",
            "SPiCEinter_2thr", "deltaMES", "BP", "mutInPBarea", "deltaESRscore", "posCryptMut", "sstypeCryptMut", "probaCryptMut", "classProbaCryptMut",
            "nearestSStoCrypt", "nearestPosSStoCrypt", "nearestDistSStoCrypt", "posCryptWT", "probaCryptWT", "classProbaCryptWT", "posSSPhysio",
            "probaSSPhysio", "classProbaSSPhysio", "probaSSPhysioMut", "classProbaSSPhysioMut",sep="\t"),con = output,sep="\n")
    }

    # Start to read the variants

    rawInput <- readLines(input, n=maxLines) # read lines by block, starting right after the header
    # lengthDataLine <- length(unlist(strsplit(rawInput[1],split="\t")))
    local_dataRefSeq <- NULL

    while(length(rawInput) != 0) { # read all file # only supports VCF format
        # EOF in last iteration -> empty input -> stop
        cat(paste("Reading lines",firstLine,"to",firstLine + length(rawInput),"\n"))
        cat(paste("\n",sub("CET",":",Sys.time(),fixed=T),"Extracting sequences and variant informations","\n"))
        total <- length(rawInput)
        if(printProcess){pb <- txtProgressBar(min = 0, max = total, initial = 1, char = "=", style = 3)}

        tmpVCF <-foreach (i=1:length(rawInput),.errorhandling='pass') %dopar% {
            readVCF(rawInput[i],i)
        }

        # merging the results of mcmapply into a dataframe
        unlisted <- as.data.frame(matrix(unlist(tmpVCF), ncol = 5+lengthDataLine, byrow = TRUE))
        data <- unlisted[,c(1,2,3,5)]
        names(data) <- c("varID", "seqPhysio", "seqMutated", "altUsed")
        transcriptsDF <- data.frame(trID = unlisted[,c(4)])
        VCFinfo_DF <- unlisted[,c(6:(5+lengthDataLine))]
        VCFinfo_text <- apply(VCFinfo_DF,1,paste,collapse="\t")

        ###################
        # SCORE CALCULATION
        ###################

        if(!is.null(data)){
            # STEP 1 : deleting repetitions in the vcf storage
            colnames(VCFinfo_DF) <- columnsNames
            VCFinfo_DF = unique(VCFinfo_DF)
            VCFinfo_toPrint <- unique(VCFinfo_text)

            # STEP 2 : update of the local_dataRefSeq
            if (!is.null(transcriptsDF)) {
                setTranscript(transcriptsDF) # add transcripts information to local_dataRefSeq
            }

            # STEP 3 : Score computation
            cat(paste("\n",gsub("CET",":",Sys.time(),fixed=T),"Score Calculation...","\n"))
            total <- nrow(data)

            if(printProcess){pb2 <- txtProgressBar(min = 0, max = total, initial = 1, char = "=", style = 3)}

            rawResult<-foreach (i=1:nrow(data),.errorhandling='pass') %dopar% {
                getOutputToSPiPmodel(data[i,"varID"],i,data[i,"seqPhysio"], data[i,"seqMutated"])
            }
            rawAnnotation = as.data.frame(matrix(unlist(rawResult),ncol=35, byrow = TRUE))
            names(rawAnnotation) <- c("chr", "strand", "gNomen", "varType", "ntChange", "ExonInfo", "exonSize", "transcript",
                "gene", "NearestSS", "DistSS", "RegType", "seqPhysio", "seqMutated", "SPiCEproba", "SPiCEinter_2thr", "deltaMES", "BP", "mutInPBarea",
                "deltaESRscore", "posCryptMut", "sstypeCryptMut", "probaCryptMut", "classProbaCryptMut", "nearestSStoCrypt", "nearestPosSStoCrypt",
                "nearestDistSStoCrypt", "posCryptWT", "probaCryptWT", "classProbaCryptWT", "posSSPhysio", "probaSSPhysio", "classProbaSSPhysio",
                "probaSSPhysioMut", "classProbaSSPhysioMut")
            data = cbind(data,rawAnnotation)
            data = SPiP(data)

            #################
            # WRITING RESULTS
            #################

            cat(paste("\n",sub("CET",":",Sys.time(),fixed=T),"Write results...","\n"))
            INFO <- apply(data[,c("altUsed","varID","Interpretation", "InterConfident", "SPiPscore", "strand", "gNomen", "varType", "ntChange", "ExonInfo", "exonSize", "transcript",
                "gene", "NearestSS", "DistSS", "RegType", "SPiCEproba", "SPiCEinter_2thr", "deltaMES", "BP", "mutInPBarea",
                "deltaESRscore", "posCryptMut", "sstypeCryptMut", "probaCryptMut", "classProbaCryptMut", "nearestSStoCrypt", "nearestPosSStoCrypt",
                "nearestDistSStoCrypt", "posCryptWT", "probaCryptWT", "classProbaCryptWT", "posSSPhysio", "probaSSPhysio", "classProbaSSPhysio",
                "probaSSPhysioMut", "classProbaSSPhysioMut")],1,paste,collapse="|")
            if(length(grep(".",data$RegType,fixed=T))>0){
                INFO[data$RegType=="."] = paste0(data$alt_used[data$RegType=="."],"|","Variant caused an error in SPiP execution")
            }
            if(printVCF){
                spipResult = mapply(mergeSPiPresult,VCFinfo_toPrint)
                if(length(grep(".",VCFinfo_DF$INFO))>0){
                    VCFinfo_DF$INFO[grep(".",VCFinfo_DF$INFO)] = spipResult[grep(".",VCFinfo_DF$INFO)]
                    VCFinfo_DF$INFO[-grep(".",VCFinfo_DF$INFO)] = paste(VCFinfo_DF$INFO[-grep(".",VCFinfo_DF$INFO)],spipResult[-grep(".",VCFinfo_DF$INFO)],sep=";")
                }else{
                    VCFinfo_DF$INFO = paste(VCFinfo_DF$INFO,spipResult,sep=";")
                }
                rawToprint = apply(VCFinfo_DF,1,paste0,collapse="\t")
                writeLines(rawToprint, con = output,sep="\n")
            }else{
                rawResult = apply(data[,c("varID", "Interpretation", "InterConfident", "SPiPscore", "strand",
                    "gNomen", "varType", "ntChange", "ExonInfo", "exonSize", "transcript", "gene", "NearestSS", "DistSS", "RegType", "SPiCEproba",
                    "SPiCEinter_2thr", "deltaMES", "BP", "mutInPBarea", "deltaESRscore", "posCryptMut", "sstypeCryptMut", "probaCryptMut", "classProbaCryptMut",
                    "nearestSStoCrypt", "nearestPosSStoCrypt", "nearestDistSStoCrypt", "posCryptWT", "probaCryptWT", "classProbaCryptWT", "posSSPhysio",
                    "probaSSPhysio", "classProbaSSPhysio", "probaSSPhysioMut", "classProbaSSPhysioMut")],1,paste,collapse="\t")
                writeLines(paste(VCFinfo_text,rawResult,sep="\t"),con = output,sep="\n")
            }
            flush(output) # necessary ?
        }

        # read the following lines = while loop increment
        rawInput <- readLines(input, n=maxLines)
        firstLine <- firstLine + maxLines
    }

}else{
    firstLine <- 1 # index of current line
    rawInput <- readLines(input, n=firstLine) # line by line # first line
    columnsNamesLine <- rawInput # table names
    columnsNames <- c(unlist(strsplit(columnsNamesLine,split='\t')), "Interpretation", "InterConfident", "SPiPscore", "strand",
        "gNomen", "varType", "ntChange", "ExonInfo", "exonSize", "transcript", "gene", "NearestSS", "DistSS", "RegType", "SPiCEproba",
        "SPiCEinter_2thr", "deltaMES", "BP", "mutInPBarea", "deltaESRscore", "posCryptMut", "sstypeCryptMut", "probaCryptMut", "classProbaCryptMut",
        "nearestSStoCrypt", "nearestPosSStoCrypt", "nearestDistSStoCrypt", "posCryptWT", "probaCryptWT", "classProbaCryptWT", "posSSPhysio",
        "probaSSPhysio", "classProbaSSPhysio", "probaSSPhysioMut", "classProbaSSPhysioMut")
    lengthDataLine <- length(columnsNames)
    if(length(grep("varID",columnsNames))==0){
        message("###########################")
        message("#Your data doesn't have the varID column")
        message("###########################")
        message(helpMessage)
        stop()
    }
    if(printVCF){
         writeLines(c(headerHelp,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"), con = output,sep="\n")
    }else{
         writeLines(c(paste(columnsNames,collapse="\t")), con = output,sep="\n")
    }

    rawInput <- readLines(input, n=maxLines) # read lines by block, starting right after the header
    while(length(rawInput) != 0) { # read all file
        cat(paste("Reading lines",firstLine,"to",firstLine + length(rawInput),"\n"))
        data = splitRawToTable(rawInput,head=FALSE)
        names(data)=unlist(strsplit(columnsNamesLine,split='\t'))
        transcriptsDF <- data.frame(trID = unlist(mapply(getTranscript2varID,data$varID)))
        setTranscript(transcriptsDF) # add transcripts information to local_dataRefSeq

        message(paste("\n",gsub("CET",":",Sys.time(),fixed=T),"Score Calculation..."))
        total <- nrow(data)
        if(printProcess){pb2 <- txtProgressBar(min = 0, max = total, initial = 1, char = "=", style = 3)}
        rawAnnotation <- foreach (i=1:nrow(data),.errorhandling='pass') %dopar% {
            getOutputToSPiPmodel(data[i,"varID"],i)
        }
        rawAnnotation = as.data.frame(matrix(unlist(rawAnnotation),ncol=35, byrow = TRUE))
        names(rawAnnotation) <- c("chr", "strand", "gNomen", "varType", "ntChange", "ExonInfo", "exonSize", "transcript",
            "gene", "NearestSS", "DistSS", "RegType", "seqPhysio", "seqMutated", "SPiCEproba", "SPiCEinter_2thr", "deltaMES", "BP", "mutInPBarea",
            "deltaESRscore", "posCryptMut", "sstypeCryptMut", "probaCryptMut", "classProbaCryptMut", "nearestSStoCrypt", "nearestPosSStoCrypt",
            "nearestDistSStoCrypt", "posCryptWT", "probaCryptWT", "classProbaCryptWT", "posSSPhysio", "probaSSPhysio", "classProbaSSPhysio",
            "probaSSPhysioMut", "classProbaSSPhysioMut")
        data = cbind(data,rawAnnotation)
        data = SPiP(data)

        message(paste("\n",sub("CET",":",Sys.time(),fixed=T),"Write results..."))
        if(printVCF){
            rawOutputVCF <- foreach (i=1:nrow(data),.errorhandling='pass') %dopar% {
                paste(convertLine2VCF(data[i,]),collapse="\t")
            }
            writeLines(unlist(rawOutputVCF), con = output,sep="\n")
        }else{
            writeLines(c(apply(data[,columnsNames],1,paste,collapse="\t")), con = output, sep = "\n")
        }
        flush(output)
        # read the following lines = while loop increment
        rawInput <- readLines(input, n=maxLines)
        firstLine <- firstLine + maxLines
    }
}
