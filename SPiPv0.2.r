#import librairy
tryCatch({
library(RCurl)
},
	error=function(cond) {
		message("Here's the original error message:")
		message(cond)
		message("*****You need to install \'RCurl\' library")
})

myOpts <- curlOptions(connecttimeout = 10)
useOption = rep("Yes",5)
names(useOption) = c("SPiCE","MES","BP","Cryptic","ESR")
samPath=NULL
fastaFile=NULL
genome="hg19"
#SPiP arguments

helpMessage="Usage: SPiPv0.2.r\n
    [Mandatory] \n
        \t[-I|--input /path/to/inputFile] \n
        \t\tlist of variants file (.txt or .vcf)\n
        \t[-O|--output /path/to/outputFile] \n
        \t\tName of ouput file (.txt)\n
    [Options] \n
        \t[-g|--GenomeAssenbly hg19] \n
        \t\tGenome assembly version (hg19 or hg38) [default= hg19]\n
        \t[-s|--SamPath /path/to/samtools] \n
        \t\tPath to samtools, if you want to use Ensembl api keep this argument empty\n
        \t[-f|--fastaGenome /path/to/fastaGenome] \n
        \t\tfasta file of genome used by samtools\n
        \t[-c|--SPiCE]\n
        \t\tUse of consensus prediction (SPiCE) (Yes/No) [default= Yes]\n
        \t[-m|--MES]\n
        \t\tUse of MES prediction in {-20; -13} (Yes/No) [default= Yes]\n
        \t[-b|--BPP]\n
        \t\tUse of Branch point prediction (Yes/No) [default= Yes]\n
        \t[-d|--Cryptic]\n
        \t\tUse of de Novo/cryptic prediction (Yes/No) [default= Yes]\n
        \t[-e|--ESR]\n
        \t\tUse of ESR prediction (Yes/No) [default= Yes]\n
        \t[-h|--help]\n
        \t\tprint this help message and exit\n
   You could : Rscript SPiPv0.2.r -I ./testCrypt.txt -O ./outTestCrypt.txt"

#get script argument
args <- commandArgs(trailingOnly = TRUE)

if (length(args)<2){message(helpMessage);stop()}

i=1
while (i < length(args)){
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
    }else if(args[i]=="-c"|args[i]=="--SPiCE"){
        useOption["SPiCE"]=args[i+1];i = i+2
    }else if(args[i]=="-m"|args[i]=="--MES"){
        useOption["MES"]=args[i+1];i = i+2
    }else if(args[i]=="-b"|args[i]=="--BPP"){
        useOption["BP"]=args[i+1];i = i+2
    }else if(args[i]=="-d"|args[i]=="--Cryptic"){
        useOption["Cryptic"]=args[i+1];i = i+2
    }else if(args[i]=="-e"|args[i]=="--ESR"){
        useOption["ESR"]=args[i+1];i = i+2
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

if(is.null(samPath)){
	useEnsemblAPI="YES"
}else{
	useEnsemblAPI="NO"
	if(is.null(fastaFile)){
		message("###########################")
		message("#No fasta file for samtools")
		message("###########################")
        message(helpMessage)
		stop()
	}
}

message('Your options are:')
message(paste(names(useOption),collapse="\t"))
message(paste(useOption,collapse="\t"))

#Get Ref files
wd_R=getwd()
inputref = paste(wd_R, "/RefFiles",sep="")

load(paste(wd_R, "/RefFiles/RefFiles.RData",sep=""))
load(paste(wd_R, "/RefFiles/dataRefSeq",genome,".RData",sep=""))

dataRefSeq = dataRefSeq[-grep("_",dataRefSeq$V1),]

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

#import data

getTranscriptFromVCF <- function(chrom, pos){
	if(substr(as.character(chrom),1,3)!="chr"){
		chrom = paste("chr",chrom,sep="")
	}
	transcript = as.character(dataRefSeq[dataRefSeq$V1==chrom & dataRefSeq$V2<=pos & dataRefSeq$V3>=pos,'V4'])
	if(length(transcript)==0){
		transcript = "no transcript"
		return(transcript)
	}else{
		return(transcript)
	}
}

getMutationFromVCF <- function (REF, ALT){
	if(as.numeric(regexpr(',',ALT,fixed =TRUE))<0){
		if(REF=='.'){
			mut = paste('ins',ALT,sep="")
		}else if(ALT=='.'){
			mut = 'del'
		}else{
			if(nchar(REF)>1 | nchar(ALT)>1){
				mut = paste('delins',ALT,sep="")
			}else{
				mut = paste(REF,ALT,sep=">")
			}
		}
	}else{
		ALT = unlist(strsplit(ALT,',',fixed = TRUE))
		mut=NULL
		for (i in 1:length(ALT)){
			if(REF=='.'){
				mut[i] = paste('ins',ALT[i],sep="")
			}else if(ALT[i]=='.'){
				mut[i] = 'del'
			}else{
				if(nchar(REF)>1 | nchar(ALT[i])>1){
					mut[i] = paste('delins',ALT[i],sep="")
				}else{
					mut[i] = paste(REF,ALT[i],sep=">")
				}
			}
		}
	}
	return(mut)
}

getPositionFromVCF <- function(POS, REF){
	if(nchar(REF)>1){
		start= POS
		end = POS+(nchar(REF)-1)
		pos = paste('g.',paste(start,end,sep='_'),sep="")
	}else{
		pos = paste('g.',POS,sep="")
	}
	return(pos)
}

getVariantFromVCF <- function(transcript,position,mutation){
	transPos = paste(unlist(transcript),unlist(position),sep=":")
	mutation = unlist(mutation)
	transPosAdj = rep(transPos,length(mutation))
	variant = paste(transPosAdj,rep(mutation,each = length(transPos)),sep=":")
	return(variant)
}

readVCF <- function(dataLine){
	dataLine = unlist(strsplit(dataLine,split='\t')) #c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
	transcript = getTranscriptFromVCF(dataLine[1],as.numeric(dataLine[2]))
	position = getPositionFromVCF(as.numeric(dataLine[2]),dataLine[4])
	mutation =getMutationFromVCF(dataLine[4],dataLine[5])
	variant = getVariantFromVCF(transcript,position,mutation)
	result <<- list(variant)
}

fileFormat = tolower(substr(basename(inputFile),nchar(basename(inputFile))-2,nchar(basename(inputFile))))

if(fileFormat=="txt"){
	data = read.table(inputFile,sep="\t",header=T)
	if(is.null(data$varID)){
		message("###########################")
		message("#Your data doesn't have the varID column")
		message("###########################")
		print_help(opt_parser)
		stop()
	}
}else if(fileFormat=="vcf"){
	data=NULL
	VCF =readLines(inputFile)
	VCF = VCF[-grep("#",VCF,fixed=TRUE)]
	tmp = mapply(readVCF,VCF)
	data$varID = as.character(unlist(tmp))
	data$annot = gsub('\t',' ',names(unlist(tmp)),fixed=TRUE)
	data= as.data.frame(data)
	if (length(data[grep('no transcript',data[,'varID']),'varID'])>0){
		message(paste( "I find no transcrit for the mutation:",paste(data[grep('no transcript',data[,'varID']),'varID'],collapse=", ")))
		data = data[-grep('no transcript',data[,'varID']),]
	}
}else{
	message("###########################")
	message("#Incorrect format of input, please try again with a txt or vcf file")
	message("###########################")
	print_help(opt_parser)
	stop()
}

data$Interpretation <- NA
data$InterConfident <- NA
data$chr <- NA
data$strand <- NA
data$gNomen <- NA
data$seqPhysio <- NA
data$seqMutated <- NA
data$NearestSS <- NA
data$distSS <- NA
data$RegType <- NA
data$SPiCEproba <- NA
data$SPiCEinter_2thr <- NA
data$deltaMES <- NA
data$mutInPBarea <- NA
data$deltaESRscore <- NA
data$posCryptMut <- NA
data$sstypeCryptMut <- NA
data$nearestSStoCrypt <- NA
data$nearestPosSStoCrypt <- NA
data$nearestDistSStoCrypt <- NA
data$probaCryptMut <- NA
data$classProbaCryptMut <- NA
data$posCryptWT <- NA
data$probaCryptWT <- NA
data$classProbaCryptWT <- NA
data$posSSPhysio <- NA
data$probaSSPhysio <- NA
data$classProbaSSPhysio <- NA
data$probaSSPhysioMut <- NA
data$classProbaSSPhysioMut <- NA

#functions used in this pipeline

getPosSSphysio <- function(transcrit){

	if(dim(dataRefSeq[dataRefSeq$V4==transcrit,])[1]==0){
		tkmessageBox(message = paste("I don't find the transcript:",transcrit,"in the Refseq database"),
		icon = "error",type="ok", title="SPiCE alert")
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
				RegType <<- "Exon"
				distSS1 <<- distSS1 -1
				if(length(varPos)==2){
					distSS2 = distSS2 -1
				}
			}else{
				nearestPosAcc = posAcc[which(posDon==nearestPosDon)]
				RegType <<- "Intron"
			}
		}else{
			distSS1 <<- nearestPosDon - varPos1+1
			if(length(varPos)==2){
				distSS2 = nearestPosDon - varPos2+1
			}
			if(varPos1 <= nearestPosDon){
				nearestPosAcc = posAcc[which(posDon==nearestPosDon)]
				RegType <<- "Intron"
			}else{
				nearestPosAcc = posAcc[which(posDon==nearestPosDon)+1]
				RegType <<- "Exon"
				distSS1 <<- distSS1 -1
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
				RegType <<- "Intron"
			}else{
				nearestPosDon = posDon[which(posAcc==nearestPosAcc)+1]
				RegType <<- "Exon"
				distSS1 <<- distSS1 + 1
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
				RegType <<- "Exon"
				distSS1 <<- distSS1 + 1
				if(length(varPos)==2){
					distSS2 = distSS2 + 1
				}
			}else{
				nearestPosDon = posDon[which(posAcc==nearestPosAcc)]
				RegType <<- "Intron"
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
					abs(abs(posVar)-abs(dataConvert$cStart[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar])-1)
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
				gVar <<- dataConvert$gStart[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar]+
			(posVar-dataConvert$cStart[dataConvert$cStart<=posVar & dataConvert$cEnd>=posVar])
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
	lenSeq = nchar(sequence)
	seqRev = paste(rep("N",lenSeq),collapse="")

	for(i in 1:lenSeq){
		substr(seqRev,i,i) = inverseDic[substr(sequence,lenSeq-i+1,lenSeq-i+1),1]
	}
	return(seqRev)
}

getSeqFromSamtools <- function(chromosome,start,end,strand,command,fastaRef){
	posToAsk <- paste(chromosome,':',start,'-',end,sep="")
	path2script <- paste('faidx',fastaRef,posToAsk)
	samtoolsSeq = system2(command, args=path2script, stdout = TRUE)
	samtoolsSeq = samtoolsSeq[-1]

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
	return(seqDNA)
}

getSequencePhysio <- function(genome,sens,chr,start,end){
	if(useEnsemblAPI=="YES"){
		if(genome=="hg38"){
			server = "http://rest.ensembl.org/"
		}else if(genome=="hg19"){
		server = "http://grch37.rest.ensembl.org/"
		}

		if(sens=="+"){
			urls = paste(server ,"sequence/region/human/",chr,":",start,"..",end,":1.fasta",sep="")

		}else if(sens=="-"){
			urls = paste(server ,"sequence/region/human/",chr,":",start,"..",end,":-1.fasta",sep="")
		}
		seqFasta=""
		i = 0
		while(seqFasta==""){
			try(eval(parse(text="seqFasta = getURL(urls, .opts = myOpts)")))
			i = i+1
			if(seqFasta==""){
				Sys.sleep(3)
				print(i)
			}else{
				seqDNA <<- unlist(strsplit(as.character(seqFasta),"\n"))[-1]
			}
		}
		if(length(seqDNA)>1){
			seqDNApool=""
			for(i in 1:length(seqDNA)){
				seqDNApool=paste(seqDNApool,seqDNA[i],sep="")
			}
			seqDNA = seqDNApool
		}else{
			seqDNA = seqDNA
		}
	}else{
		command = samPath
		fastaRef = fastaFile
		gPos=c(start,end)
		gPos = gPos[order(gPos)]
		seqDNA = getSeqFromSamtools(chr, gPos[1], gPos[2], sens, command, fastaRef)
	}
	return(seqDNA)
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
	if(sstype=="Acc"){
		scoreMES = MESacceptor(seqCon)
	}else{
		if(substr(as.vector(seqCon),4,5)=="GT"){
			scoreMES = MESdonor(seqCon)
		}else{
			scoreMES = 0
		}
	}
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
	SPiCEproba <<- round(exp(-3.59-8.21*round(deltaMES,3)-32.30*round(deltaSSF,3))/(1+exp(-3.59-8.21*round(deltaMES,3)-32.30*round(deltaSSF,3))),5)
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
		mutInPBareaBPP = NA
	}else{
		if(length(varPos)==1){
			if(nrow(tmpAnnotBP[tmpAnnotBP$start<=varPos & tmpAnnotBP$end>=varPos,])==1){
				posBP = tmpAnnotBP[tmpAnnotBP$start<=varPos & tmpAnnotBP$end>=varPos,"posBP"]
				distToAcc = tmpAnnotBP[tmpAnnotBP$start<=varPos & tmpAnnotBP$end>=varPos,"distToAcc"]
				score = tmpAnnotBP[tmpAnnotBP$start<=varPos & tmpAnnotBP$end>=varPos,"score"]
				mutInPBareaBPP = paste("Yes, g.",posBP," (",distToAcc,"): ",round(score,3), sep = "")
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
				mutInPBareaBPP = paste("Yes, g.",posBP," (",distToAcc,"): ",round(score,3), sep = "")
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
	print(varDecomp)

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

getGlobaInterpretation <- function(SPiCEinterpret, RegType, deltaMES, mutInPBareaBPP, classProbaMut, distSS, deltaESR){
	interpretFinal = NULL
	if(SPiCEinterpret=="high" | SPiCEinterpret=="medium"){
		interpretFinal = c(interpretFinal,"Alter by SPiCE")
	}
	if (length(grep("PolyTC",RegType ))>0 & deltaMES<(-0.15)){
		interpretFinal = c(interpretFinal,"Alter by MES (Poly TC)")
	}
	if (length(grep("BP",RegType ))>0 & mutInPBareaBPP!="No"){
		interpretFinal = c(interpretFinal,"Alter BP")
	}
	if (classProbaMut=="Yes"){
		if(abs(distSS)>=150 & length(grep("Intron",RegType ))>0) {
			interpretFinal =c(interpretFinal, "Alter by create Exon")
		}else{
			interpretFinal = c(interpretFinal,"Alter by create Cryptic")
		}
	}
	if (length(grep("Exon",RegType ))>0 & abs(distSS)<120){
		if(deltaESR<(-1.25)) {
			interpretFinal = c(interpretFinal,"Alter ESR")
		}else{
			interpretFinal = c(interpretFinal,"NTR ESR")
		}
	}
	if(is.null(interpretFinal)){
		interpretFinal = "NTR"
	}else if (length(interpretFinal)>1){
		interpretFinal = paste(interpretFinal[interpretFinal!="NTR ESR"], collapse=" + ")
	}
	return(interpretFinal)
}

getPredConfident <- function(interpretFinal, RegType){
	#proba from SNP + UV analysis
	probaInter = -1
	if (length(grep("Exon",RegType ))>0){
		if(length(grep("Alter by SPiCE",interpretFinal ))>0){
			probaInter = "80.1 % +/- 6.8 %"
		}else if (length(grep("Alter ESR",interpretFinal ))>0){
			probaInter = "24.1 % +/- 4.3 %"
		}else if (length(grep("Alter by create Cryptic",interpretFinal ))>0){
			probaInter = "15.9 % +/- 3.0 %"
		}else if(length(grep("NTR ESR",interpretFinal ))>0){
			probaInter = "05.5 % +/- 1.4 %"
		}else if(length(grep("NTR",interpretFinal ))>0){
			probaInter = "01.1 % +/- 0.8 %"
		}
	}else{
		if(length(grep("Alter by SPiCE",interpretFinal ))>0){
			probaInter = "92.9 % +/- 2.1 %"
		}else if (length(grep("Alter by MES",interpretFinal ))>0){
			probaInter = "75.0 % +/- 12.3 %"
		}else if (length(grep("Alter BP",interpretFinal ))>0){
			probaInter = "52.4 % +/- 12.5 %"
		}else if(length(grep("Alter by create Cryptic",interpretFinal ))>0){
			probaInter = "08.4 % +/- 2.6 %"
		}else if(length(grep("Alter by create Exon",interpretFinal ))>0){
			probaInter = "02.3 % +/- 0.5 %"
		}else if(length(grep("NTR",interpretFinal ))>0){
			probaInter = "00.2 % +/- 0.1 %"
		}
	}
	return (probaInter)
}

getOutput <- function(data,i){
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

	varType <<- varType
	tmpTableSeqNoPhyMut = tmpTableSeq[tmpTableSeq$Physio!="Yes" & tmpTableSeq$seqType == "Mut",]

	data$chr[i] <<- chr
	data$strand[i] <<- sens
	data$gNomen[i] <<- varPos[1]
	data$NearestSS[i] <<- SstypePhy
	if(length(varPos)==1){
		data$distSS[i] <<- distSS
	}else if(length(varPos)==2){
		if(distSS[1]==distSS2){
			data$distSS[i] <<- distSS[1]
		}else{
			data$distSS[i] <<- paste(distSS, distSS2,sep="; ")
		}
	}else{
		print("erreur varpos")
	}
	data$RegType[i] <<- RegType
	data$seqPhysio[i] <<- seqPhysio
	data$seqMutated[i] <<- seqMutated
	data$SPiCEproba[i] <<- SPiCEproba
	data$SPiCEinter_2thr[i] <<- SPiCEinter_2thr
	data$deltaMES[i] <<- deltaMES
	data$mutInPBarea[i] <<- mutInPBareaBPP
	data$deltaESRscore[i] <<- ESRscore

	if(nrow(tmpTableSeqNoPhyMut)>0){
		tmpTableSeqMAXMut = tmpTableSeqNoPhyMut[tmpTableSeqNoPhyMut$proba==max(tmpTableSeqNoPhyMut$proba),]
		sstypCrypt = tmpTableSeqMAXMut$sstype[1]
		tmpTableSeqNoPhyWT = tmpTableSeq[tmpTableSeq$Physio!="Yes" & tmpTableSeq$seqType == "WT"& tmpTableSeq$sstype==sstypCrypt,]
		tmpTableSeqMAXWT = tmpTableSeqNoPhyWT[tmpTableSeqNoPhyWT$proba==max(tmpTableSeqNoPhyWT$proba),]

		tmpTableSeqPhy = tmpTableSeq[tmpTableSeq$Physio=="Yes"& tmpTableSeq$seqType == "WT" & tmpTableSeq$sstype==sstypCrypt,]
		tmpTableSeqPhyMut = tmpTableSeq[tmpTableSeq$Physio=="Yes"& tmpTableSeq$seqType == "Mut" & tmpTableSeq$sstype==sstypCrypt,]

		getNearestPosCrypt(sens, as.numeric(tmpTableSeqMAXMut$pos) ,posDon, posAcc)
		data$posCryptMut[i] <<- tmpTableSeqMAXMut$pos
		data$sstypeCryptMut[i] <<- as.character(tmpTableSeqMAXMut$sstype)
		data$nearestSStoCrypt[i] <<- SstypePhyCrypt
		data$nearestPosSStoCrypt[i] <<- nearestPosPhyCrypt
		data$nearestDistSStoCrypt[i] <<- distSScrypt
		data$probaCryptMut[i] <<- tmpTableSeqMAXMut$proba
		data$classProbaCryptMut[i] <<- tmpTableSeqMAXMut$classProba
		if(nrow(tmpTableSeqMAXWT)>0){
			data$posCryptWT[i] <<- tmpTableSeqMAXWT$pos[1]
			data$probaCryptWT[i] <<- tmpTableSeqMAXWT$proba[1]
			data$classProbaCryptWT[i] <<- tmpTableSeqMAXWT$classProba[1]
		}
		if(nrow(tmpTableSeqPhy)>0){
			data$posSSPhysio[i] <<- tmpTableSeqPhy$pos[1]
			data$probaSSPhysio[i] <<- tmpTableSeqPhy$proba[1]
			data$classProbaSSPhysio[i] <<- tmpTableSeqPhy$classProba[1]
			data$probaSSPhysioMut[i] <<- tmpTableSeqPhy$proba[1]
			data$classProbaSSPhysioMut[i] <<- tmpTableSeqPhy$classProba[1]
		}
		if(nrow(tmpTableSeqPhyMut)>0){
			data$probaSSPhysioMut[i] <<- tmpTableSeqPhyMut$proba[1]
			data$classProbaSSPhysioMut[i] <<- tmpTableSeqPhyMut$classProba[1]
		}
	}else{
		data$posCryptMut[i] <<- 0
		data$sstypeCryptMut[i] <<- "No site"
		data$nearestSStoCrypt[i] <<- "No site"
		data$nearestPosSStoCrypt[i] <<- 0
		data$nearestDistSStoCrypt[i] <<- 0
		data$probaCryptMut[i] <<- 0
		data$classProbaCryptMut[i] <<- "No"
		data$posCryptWT[i] <<- 0
		data$probaCryptWT[i] <<- 0
		data$classProbaCryptWT[i] <<- "No"
		data$posSSPhysio[i] <<- 0
		data$probaSSPhysio[i] <<- 0
		data$classProbaSSPhysio[i] <<- "No"
		data$probaSSPhysioMut[i] <<- 0
		data$classProbaSSPhysioMut[i] <<- "No"
	}
	interpretation <- getGlobaInterpretation(SPiCEinter_2thr, RegType, deltaMES, mutInPBareaBPP, data$classProbaCryptMut[i], distSS[1], ESRscore)
	data$Interpretation[i] <<- interpretation
	data$InterConfident[i] <<- getPredConfident(interpretation, RegType)
}

#launch analysis

T1 <- as.numeric(format(Sys.time(), "%s"))

for (i in 1:nrow(data)){
	tryCatch({
		print(i)
		getVariantInfo(as.character(data$varID[i]))
		},
	error=function(cond) {
		message(paste("Variant caused a error:", data[i,"varID"]))
		message("Here's the original error message of getVariantInfo function:")
		message(cond)
	})
	tryCatch({
		getOutput(data,i)
		},
	error=function(cond) {
		message(paste("Variant caused a error:", data[i,"varID"]))
		message("Here's the original error message of getOutput function:")
		message(cond)
	})
}
T2 <- as.numeric(format(Sys.time(), "%s"))

print(paste("Runtime:",round((T2 - T1),3),"sec"))

write.table(data,outputFile,sep="\t",dec=".",row.names=F)
