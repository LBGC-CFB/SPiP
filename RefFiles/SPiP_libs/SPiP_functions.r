tryCatch({
source(paste0(path2scripts,"SPiP_sequences.r"))
},
    error=function(cond) {
        message("Here's the original error message:")
		message(cond)
		message("*****Error while importing source file SPiP_sequences.r")
    }
)

tryCatch({
source(paste0(path2scripts,"SPiP_scores.r"))
},
    error=function(cond) {
        message("Here's the original error message:")
		message(cond)
		message("*****Error while importing source file SPiP_scores.r")
    }
)

setTranscript <- function (transcriptsDF) {
	transcriptsDF <- unique(transcriptsDF)
    transcript <- as.character(transcriptsDF[,"trID"])
    transcript <- transcript[transcript!=""]
	if (length(transcript)>0) {
		local_dataRefSeq <<- dataRefSeq[which(dataRefSeq$V4 %in% transcript),]
	}
}

getTranscript2varID <- function(varID){
    tmp = unlist(strsplit(as.character(varID),":",fixed = T))
    tr = tmp[1]
    if(length(grep(".",tr ))>0){
        tr = unlist(strsplit(tr,".",fixed = T))[1]
    }
    return(tr)
}

contigToChr <- function(text){
    text1=unlist(strsplit(text, ".", fixed = TRUE))[1]
    chr=text
    if(substr(text1,1,2)=="NC"){
        chr=paste('chr',as.numeric(substr(text1,4,nchar(text1))),sep="")
        if(chr=="chr23"){chr="chrX"}else if(chr=="chr24"){chr="chrY"}
    }
    return(chr)
}

getPosSSphysio <- function(transcrit){

	if(dim(local_dataRefSeq[local_dataRefSeq$V4==transcrit,])[1]==0){
		paste("I don't find the transcript:",transcrit,"in the Refseq database")
	}else{

		chr <<- as.character(local_dataRefSeq[local_dataRefSeq$V4==transcrit,1])
		sens <<- as.character(local_dataRefSeq[local_dataRefSeq$V4==transcrit,6])

	if(sens=="+"){
		posStart=local_dataRefSeq[local_dataRefSeq$V4==transcrit,2]
		tailleCum=local_dataRefSeq[local_dataRefSeq$V4==transcrit,12]
		tailleCum=strsplit(as.character(tailleCum),split=",")
		tailleCum=as.numeric(unlist(tailleCum))
		posAcc <<- posStart+tailleCum
		taille=local_dataRefSeq[local_dataRefSeq$V4==transcrit,11]
		taille=strsplit(as.character(taille),split=",")
		taille=as.numeric(unlist(taille))
		posDon <<- posAcc+taille
		if(length(posDon)>1 & length(posAcc)>1){
			posDon <<- posDon[-length(posDon)]
			posAcc <<- posAcc[-1]
		}
	}else if(sens=="-"){
		posEnd=local_dataRefSeq[local_dataRefSeq$V4==transcrit,2]
		tailleCum=local_dataRefSeq[local_dataRefSeq$V4==transcrit,12]
		tailleCum=strsplit(as.character(tailleCum),split=",")
		tailleCum=as.numeric(unlist(tailleCum))
		posDon <<- posEnd+tailleCum
		taille=local_dataRefSeq[local_dataRefSeq$V4==transcrit,11]
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

	sens = as.character(local_dataRefSeq[local_dataRefSeq$V4==transcrit,6])
	posStart=local_dataRefSeq[local_dataRefSeq$V4==transcrit,2]
	tailleExon = as.numeric(unlist(strsplit(as.character(local_dataRefSeq[local_dataRefSeq$V4==transcrit,11]),",")))
	tailleCum=local_dataRefSeq[local_dataRefSeq$V4==transcrit,12]
	tailleCum=strsplit(as.character(tailleCum),split=",")
	tailleCum=as.numeric(unlist(tailleCum))

	if(sens=="+"){

		gCDSstart = local_dataRefSeq[local_dataRefSeq$V4==transcrit,7]
		gCDSend = local_dataRefSeq[local_dataRefSeq$V4==transcrit,8]
		posAcc = posStart+tailleCum + 1
		posDon = posAcc+tailleExon - 1

		dataConvert=data.frame(idEx=c(1:length(tailleExon )),lenEx=tailleExon,gStart=posAcc,gEnd=posDon)

		if(nrow(dataConvert[dataConvert$gStart<=posVar & dataConvert$gEnd>=posVar,])>=1){
			ExonInfo = c(paste0("Exon ", dataConvert$idEx[dataConvert$gStart<=posVar & dataConvert$gEnd>=posVar]),
					dataConvert$lenEx[dataConvert$gStart<=posVar & dataConvert$gEnd>=posVar])
		}else{
			ExonInfo = c(paste0("Intron ", dataConvert$idEx[which(dataConvert$gEnd>=posVar)[1]-1]),
            				abs(dataConvert$gEnd[which(dataConvert$gEnd>=posVar)[1]-1]-
								dataConvert$gStart[which(dataConvert$gEnd>=posVar)[1]])-1)
		}
	}else if(sens=="-"){

		gCDSstart = local_dataRefSeq[local_dataRefSeq$V4==transcrit,8]
		gCDSend = local_dataRefSeq[local_dataRefSeq$V4==transcrit,7]
		posDon = posStart+tailleCum + 1
		posAcc = posDon+tailleExon - 1

		dataConvert=data.frame(idEx=c(length(tailleExon):1),lenEx=tailleExon,gStart=posAcc,gEnd=posDon)

		if(nrow(dataConvert[dataConvert$gStart>=posVar & dataConvert$gEnd<=posVar,])>=1){
			ExonInfo = c(paste0("Exon ", dataConvert$idEx[dataConvert$gStart>=posVar & dataConvert$gEnd<=posVar]),
					dataConvert$lenEx[dataConvert$gStart>=posVar & dataConvert$gEnd<=posVar])
		}else{
			ExonInfo = c(paste0("Intron ", dataConvert$idEx[which(dataConvert$gEnd>=posVar)[1]]),
							abs(dataConvert$gEnd[which(dataConvert$gEnd>=posVar)[1]]-
								dataConvert$gStart[which(dataConvert$gEnd>=posVar)[1]-1])-1)
		}
	}
	return(ExonInfo)
}

convertcNomenIngNomen <- function(transcrit,posVar){

	sens = as.character(local_dataRefSeq[local_dataRefSeq$V4==transcrit,6])
	posStart=local_dataRefSeq[local_dataRefSeq$V4==transcrit,2]
	tailleExon = as.numeric(unlist(strsplit(as.character(local_dataRefSeq[local_dataRefSeq$V4==transcrit,11]),",")))
	tailleCum=local_dataRefSeq[local_dataRefSeq$V4==transcrit,12]
	tailleCum=strsplit(as.character(tailleCum),split=",")
	tailleCum=as.numeric(unlist(tailleCum))

	if(sens=="+"){

		gCDSstart = local_dataRefSeq[local_dataRefSeq$V4==transcrit,7]
		gCDSend = local_dataRefSeq[local_dataRefSeq$V4==transcrit,8]

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

		gCDSstart = local_dataRefSeq[local_dataRefSeq$V4==transcrit,8]
		gCDSend = local_dataRefSeq[local_dataRefSeq$V4==transcrit,7]
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

getSplitTableSeq <- function(varName, chr, varPos, sens, seqPhysio, seqMutated, posDon, posAcc, nearestPosAll, varType){

	relPosAccPhy = as.numeric(gregexpr("AG",seqPhysio)[[1]])
	relPosDonPhy = as.numeric(gregexpr("GT",seqPhysio)[[1]])
	relPosAccMut = as.numeric(gregexpr("AG",seqMutated)[[1]])
	relPosDonMut = as.numeric(gregexpr("GT",seqMutated)[[1]])

	relPosAccPhyFilt = relPosAccPhy[relPosAccPhy >= 134 & relPosAccPhy <= 201]
	relPosDonPhyFilt = relPosDonPhy[relPosDonPhy >= 99 & relPosDonPhy <= 166]
	relPosAccMutFilt = relPosAccMut[relPosAccMut >= 146 & relPosAccMut <= 162]
	relPosDonMutFilt = relPosDonMut[relPosDonMut >= 144 & relPosDonMut <= 157]

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
    	tmpTableSeq$classProba[tmpTableSeq$sstype=="Don" & tmpTableSeq$proba >= 0.0262] = "Yes"
        tmpTableSeq$classProba[tmpTableSeq$sstype=="Acc" & tmpTableSeq$proba >= 0.0405] = "Yes"
        tmpPhysio = tmpTableSeq[tmpTableSeq$Physio=="Yes",]

    	tmp = merge(tmpTableSeq[tmpTableSeq$seqType=="WT",c("sstype","pos","proba","Physio")],tmpTableSeq[tmpTableSeq$seqType=="Mut",c("sstype","pos","proba","Physio")],
    					by.x = c("sstype","pos"),by.y = c("sstype","pos"),all.x=F,all.y=T)
    	tmp$proba.x[is.na(tmp$proba.x)]=0
    	tmp$Physio.x[is.na(tmp$Physio.x)] = "No"
    	tmp = tmp[order(tmp$proba.y,decreasing=T),]
    	tmp = tmp[!duplicated(tmp$pos),]

        tmpTableSeq = tmpTableSeq[which(tmpTableSeq$pos%in%tmp$pos),]
        tmpTableSeq = rbind(tmpTableSeq,tmpPhysio)
    	if(length(which(tmpTableSeq$seqType=="Mut" & tmpTableSeq$pos%in%tmp$pos[tmp$proba.x >= tmp$proba.y & tmp$Physio.x!="Yes"]))>0){
            tmpTableSeq = tmpTableSeq[-which(tmpTableSeq$seqType=="Mut" & tmpTableSeq$pos%in%tmp$pos[tmp$proba.x>=tmp$proba.y & tmp$Physio.x!="Yes"]),]
    	}
    }
	return(tmpTableSeq)
}

getDeltaESRseq <- function(SstypePhy, distSS, seqPhysio, seqMutated){
	if(abs(distSS)>120){
		ESRscore <<- 10
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
        if(length(varDecomp)==2 & varDecomp[1] != 'no transcript' & varDecomp[1] != 'no sequence'){
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
                varPos <- as.numeric(varPos)
				nbDel <<- abs(varPos[1]-varPos[2])+1
                ntIns <<- gsub("delins", "", ntChange)
                nbIns <<- nchar(ntIns)

			}else{
				if(length(grep("c.",varPos[1]))>0){
					varPos[1] = substr(varPos[1],3,nchar(varPos[1]))
					convertcNomenIngNomen(transcript, varPos[1])
					varPos = c(gVar,gVar)
				}else{
					varPos[1] = substr(varPos[1],3,nchar(varPos[1]))
					varPos = c(varPos,varPos)
				}
				nbDel <<- 1
                ntIns <<- gsub("delins", "", ntChange)
                nbIns <<- nchar(ntIns)
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
                varPos = as.numeric(varPos)
                nbDel <<- abs(varPos[1]-varPos[2])+1
			}else{
				if(length(grep("c.",varPos[1]))>0){
					varPos[1] = substr(varPos[1],3,nchar(varPos[1]))
					convertcNomenIngNomen(transcript, varPos[1])
					varPos = c(gVar,gVar)
				}else{
					varPos[1] = substr(varPos[1],3,nchar(varPos[1]))
					varPos = c(varPos,varPos)
				}
				nbDel <<- 1
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
			ntIns <<- gsub("ins", "", ntChange)
			nbIns <<- nchar(ntIns)
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
            # que faire ?
            varPos = as.numeric(varPos)
			nbIns <<- abs(varPos[1]-varPos[2])+1
		}else{
			varType = "substitution"
			if(length(grep("c.",varPos))>0){
				varPos=substr(varPos,3,nchar(varPos))
				convertcNomenIngNomen(transcript, varPos)
				varPos = gVar
			}else{
				varPos=substr(varPos,3,nchar(varPos))
			}
			ntMut <<- as.character(unlist(strsplit(ntChange,">")))[2]
		}
	}
	varID <<- varID
	varPos <<- as.numeric(varPos)
	varType <<- varType
	transcript <<- transcript
	ntChange <<- ntChange
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

getAnnotation <- function(seqPhysio, seqMutated){
	getPosSSphysio(transcript)
	getNearestPos(sens, varPos ,posDon, posAcc)

	start=varPos[1]-150
	end=varPos[1]+150
    if(fileFormat=="txt"){
        seqPhysio <- getSequencePhysio(genome,sens,chr,start,end)
        seqPhysio <<- seqPhysio
        seqMutated = getSequenceMutated(varPos, sens, seqPhysio, ntChange, varType, genome, chr)
    }
	seqPhysio <<- toupper(as.character(seqPhysio))
	seqMutated <<- toupper(as.character(seqMutated))
	tmpTableSeq <<- getSplitTableSeq(varID, chr, varPos[1], sens, seqPhysio, seqMutated, posDon, posAcc, nearestPosAll,varType)

	if(as.numeric(gregexpr("Exon",RegType))<0){
		ESRscore = 10
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
	gNomen <- if(length(varPos)==1){varPos}else{paste(varPos,collapse="_")}
	NearestSS <- SstypePhy
	if(length(varPos)==1){
		DistSS <- distSS
	}else if(length(varPos)==2){
		DistSS <- distSS[1]
	}else{
		print("erreur varpos")
	}
	gene <- as.character(local_dataRefSeq$V13[local_dataRefSeq$V4==transcript])
	mutInPBarea <- mutInPBareaBPP
    if(is.na(mutInPBareaBPP)){BP <- 0}else if(mutInPBareaBPP=="No"){BP <- 0}else{BP <- 1}
	deltaESRscore <- ESRscore

	if(nrow(tmpTableSeqNoPhyMut)>0){
		tmpTableSeqMAXMut = tmpTableSeqNoPhyMut[tmpTableSeqNoPhyMut$proba==max(tmpTableSeqNoPhyMut$proba),]
		sstypCrypt = tmpTableSeqMAXMut$sstype[1]
		tmpTableSeqNoPhyWT = tmpTableSeq[tmpTableSeq$Physio!="Yes" & tmpTableSeq$seqType == "WT"& tmpTableSeq$sstype==sstypCrypt,]

		tmpTableSeqPhy = tmpTableSeq[tmpTableSeq$Physio=="Yes"& tmpTableSeq$seqType == "WT" & tmpTableSeq$sstype==sstypCrypt,]
		tmpTableSeqPhyMut = tmpTableSeq[tmpTableSeq$Physio=="Yes"& tmpTableSeq$seqType == "Mut" & tmpTableSeq$sstype==sstypCrypt,]
        if(nrow(tmpTableSeqMAXMut)>1){tmpTableSeqMAXMut =tmpTableSeqMAXMut[1,]}

		getNearestPosCrypt(sens, as.numeric(tmpTableSeqMAXMut$pos) ,posDon, posAcc)
		posCryptMut <- tmpTableSeqMAXMut$pos
		sstypeCryptMut <- as.character(tmpTableSeqMAXMut$sstype)
		nearestSStoCrypt <- SstypePhyCrypt
		nearestPosSStoCrypt <- nearestPosPhyCrypt
		nearestDistSStoCrypt <- distSScrypt
		probaCryptMut <- tmpTableSeqMAXMut$proba
		classProbaCryptMut <- tmpTableSeqMAXMut$classProba
		if(nrow(tmpTableSeqNoPhyWT)>0){
            tmpTableSeqMAXWT = tmpTableSeqNoPhyWT[tmpTableSeqNoPhyWT$proba==max(tmpTableSeqNoPhyWT$proba),]
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
    tmp_ExonInfo <- getExonInfo(transcript,varPos[1])
    ExonInfo <- tmp_ExonInfo[1]
    exonSize <- tmp_ExonInfo[2]
    result <- c(chr, strand, gNomen, varType, ntChange, ExonInfo, exonSize, transcript,
        gene, NearestSS, DistSS, RegType, seqPhysio, seqMutated, SPiCEproba, SPiCEinter_2thr, deltaMES, BP, mutInPBarea,
        deltaESRscore, posCryptMut, sstypeCryptMut, probaCryptMut, classProbaCryptMut, nearestSStoCrypt, nearestPosSStoCrypt,
        nearestDistSStoCrypt, posCryptWT, probaCryptWT, classProbaCryptWT, posSSPhysio, probaSSPhysio, classProbaSSPhysio,
        probaSSPhysioMut, classProbaSSPhysioMut)
    return(result)
}

getOutputToSPiPmodel <- function(varID,i,seqPhysio = "", seqMutated = ""){
    if(printProcess){setTxtProgressBar(pb2, i)}
    if(as.numeric(regexpr('no transcript',varID))>0|
        as.numeric(regexpr('mutUnknown',varID))>0|
        as.numeric(regexpr('no sequence',varID))>0)    {
        tmp <- rep("NA",35)
        return(tmp)
    }else{
        tryCatch({
            getVariantInfo(as.character(varID))
            tmp <- getAnnotation(seqPhysio, seqMutated)
            return(tmp)
        },
        error=function(cond) {
            message(paste("Variant caused a error:", varID))
            tmp <- rep("NA",35)
            return(tmp)
        })
    }
}

getGlobaInterpretation <- function(SPiCEinterpret, RegType, deltaMES, mutInPBareaBPP, classProbaMut, classProbaCryptWT, distSS, deltaESR){
	if(SPiCEinterpret=="."){
		interpretFinal = NA
	}else{
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
		if (classProbaMut=="Yes"){
			if(abs(distSS)>=150 & length(grep("Intron",RegType ))>0) {
				interpretFinal =c(interpretFinal, "Alter by create New Exon")
				if(classProbaCryptWT=="Yes") {
					interpretFinal = NULL
				}
			}else{
				interpretFinal = c(interpretFinal,"Alter by create New splice site")
			}
		}
		if (length(grep("Exon",RegType ))>0 & abs(distSS)<120 & deltaESR<(-0.415)){
			interpretFinal = c(interpretFinal,"Alter ESR")
		}
		if(is.null(interpretFinal)){
			interpretFinal = "Alter by complex event"
		}else if (length(interpretFinal)>1){
			interpretFinal = paste(interpretFinal, collapse=" + ")
		}
	}
	return(interpretFinal)
}

getVPP <- function(score = 0){
    if(score==(-1)){
		probaInter = NA
	}else{
		VPP = round(VPPtable$propPos[score>=VPPtable$minScore & score<=VPPtable$maxScore][1]*100,2)
		ICmin = round(VPPtable$confint95a[score>=VPPtable$minScore & score<=VPPtable$maxScore][1]*100,2)
		ICmax = round(VPPtable$confint95b[score>=VPPtable$minScore & score<=VPPtable$maxScore][1]*100,2)
		if(VPP<10){VPP=paste0("0",VPP)}
		if(ICmin<10){ICmin=paste0("0",ICmin)}
		if(ICmax<10){ICmax=paste0("0",ICmax)}
		probaInter = paste0(VPP," % [",ICmin," % - ",ICmax," %]")
		return(probaInter)
	}
}

getVPN <- function(distSS,RegType){
    if(RegType=="."){
		probaInter = NA
	}else{
		if(length(grep("Intron",RegType))>0 & distSS>0){
			tmpVPNtable = VPNtable[VPNtable$region=="intron5p",]
		}else if(length(grep("Intron",RegType))>0 & distSS<0){
			tmpVPNtable = VPNtable[VPNtable$region=="intron3p",]
		}else if(length(grep("Exon",RegType))>0 & distSS>0){
			tmpVPNtable = VPNtable[VPNtable$region=="exon3p",]
		}else if(length(grep("Exon",RegType))>0 & distSS<0){
			tmpVPNtable = VPNtable[VPNtable$region=="exon5p",]
		}
		VPN = round(tmpVPNtable$VPN[distSS>tmpVPNtable$rangeInf & distSS<=tmpVPNtable$rangeSupp][1]*100,2)
		ICmin = round(tmpVPNtable$confint95a[distSS>tmpVPNtable$rangeInf & distSS<=tmpVPNtable$rangeSupp][1]*100,2)
		ICmax = round(tmpVPNtable$confint95b[distSS>tmpVPNtable$rangeInf & distSS<=tmpVPNtable$rangeSupp][1]*100,2)
		if(VPN<10){VPN=paste0("0",VPN)}
		if(ICmin<10){ICmin=paste0("0",ICmin)}
		if(ICmax<10){ICmax=paste0("0",ICmax)}
		probaInter = paste0(VPN," % [",ICmin," % - ",ICmax," %]")
	}
	return(probaInter)
}


getPredConfident <- function(RegType, distSS, SPiPscore, interpretation){
	#proba from SNP + UV analysis (N = 99,616)
	probaInter = -1
    if(SPiPscore>thToSPiPintron & length(grep("Intron",RegType))>0){
		probaInter = getVPP(SPiPscore)
    }else if(SPiPscore>thToSPiPexon & length(grep("Exon",RegType))>0){
		probaInter = getVPP(SPiPscore)
    }else{
		tryCatch({
			probaInter = getVPN(distSS[1],RegType)
		},
		error=function(cond) {
			probaInter = NA
		})
    }
    if(!is.na(interpretation)){
        if(SPiPscore<thToComplexEvent & interpretation=="Alter by complex event"){
            tryCatch({
                probaInter = getVPN(distSS[1],RegType)
            },
            error=function(cond) {
                probaInter = NA
            })
        }
    }else{probaInter = NA}

	return (probaInter)
}

SPiP <- function(data){

    data$deltaMES = as.numeric(as.character(data$deltaMES))
    data$BP = as.numeric(as.character(data$BP))
    data$probaSSPhysio = as.numeric(as.character(data$probaSSPhysio))
    data$probaCryptMut = as.numeric(as.character(data$probaCryptMut))
    data$DistSS = as.numeric(as.character(data$DistSS))
    data$exonSize = as.numeric(as.character(data$exonSize))
    data$probaCryptWT = as.numeric(as.character(data$probaCryptWT))
    data$deltaESRscore = as.numeric(as.character(data$deltaESRscore))
    data$SPiCEproba = as.numeric(as.character(data$SPiCEproba))
    data$probaSSPhysioMut = as.numeric(data$probaSSPhysioMut)

    data$RegTypeNum = unlist(lapply(list(data$RegType),function(x) RegTypeToNumber[x,1]))
	data$probaSSPhysioMut[is.na(data$probaSSPhysioMut)] = 0
	data$probaSSPhysio[is.na(data$probaSSPhysio)] = 0
	data$exonSize[is.na(data$exonSize)] = 0
	prediction = predict(fit.rf,newdata= data,type="prob" )
	data$SPiPscore = prediction[,2]

	# remove NAs
    data$SPiPscore[is.na(data$SPiPscore)] <- (-1)
	data$SPiCEinter_2thr[is.na(data$SPiCEinter_2thr)] <- "."
	data$RegType[is.na(data$RegType)] <- "."
	data$deltaMES[is.na(data$deltaMES)] <- (-1)
	data$mutInPBarea[is.na(data$mutInPBarea)] <- "."
	data$classProbaCryptMut[is.na(data$classProbaCryptMut)] <- "."
	data$classProbaCryptWT[is.na(data$classProbaCryptWT)] <- "."
	data$DistSS[is.na(data$DistSS)] <- (-1)
	data$deltaESRscore[is.na(data$deltaESRscore)] <- (-1)

	# get interpretation
    data$Interpretation[data$SPiPscore<=thToSPiPintron & as.numeric(regexpr("Intron",data$RegType))>0] = "NTR"
	data$Interpretation[data$SPiPscore<=thToSPiPexon & as.numeric(regexpr("Exon",data$RegType))>0] = "NTR"

	data$Interpretation[data$SPiPscore>thToSPiPintron & as.numeric(regexpr("Intron",data$RegType))>0] = unlist(mapply(getGlobaInterpretation,
							data$SPiCEinter_2thr[data$SPiPscore>thToSPiPintron & as.numeric(regexpr("Intron",data$RegType))>0],
							data$RegType[data$SPiPscore>thToSPiPintron & as.numeric(regexpr("Intron",data$RegType))>0],
							data$deltaMES[data$SPiPscore>thToSPiPintron & as.numeric(regexpr("Intron",data$RegType))>0],
							data$mutInPBarea[data$SPiPscore>thToSPiPintron & as.numeric(regexpr("Intron",data$RegType))>0],
							data$classProbaCryptMut[data$SPiPscore>thToSPiPintron & as.numeric(regexpr("Intron",data$RegType))>0],
							data$classProbaCryptWT[data$SPiPscore>thToSPiPintron & as.numeric(regexpr("Intron",data$RegType))>0],
							data$DistSS[data$SPiPscore>thToSPiPintron & as.numeric(regexpr("Intron",data$RegType))>0],
							data$deltaESRscore[data$SPiPscore>thToSPiPintron & as.numeric(regexpr("Intron",data$RegType))>0]))

	data$Interpretation[data$SPiPscore>thToSPiPexon & as.numeric(regexpr("Exon",data$RegType))>0] = unlist(mapply(getGlobaInterpretation,
							data$SPiCEinter_2thr[data$SPiPscore>thToSPiPexon & as.numeric(regexpr("Exon",data$RegType))>0],
							data$RegType[data$SPiPscore>thToSPiPexon & as.numeric(regexpr("Exon",data$RegType))>0],
							data$deltaMES[data$SPiPscore>thToSPiPexon & as.numeric(regexpr("Exon",data$RegType))>0],
							data$mutInPBarea[data$SPiPscore>thToSPiPexon & as.numeric(regexpr("Exon",data$RegType))>0],
							data$classProbaCryptMut[data$SPiPscore>thToSPiPexon & as.numeric(regexpr("Exon",data$RegType))>0],
							data$classProbaCryptWT[data$SPiPscore>thToSPiPexon & as.numeric(regexpr("Exon",data$RegType))>0],
							data$DistSS[data$SPiPscore>thToSPiPexon & as.numeric(regexpr("Exon",data$RegType))>0],
							data$deltaESRscore[data$SPiPscore>thToSPiPexon & as.numeric(regexpr("Exon",data$RegType))>0]))
    data$InterConfident = unlist(mapply(getPredConfident, data$RegType, data$DistSS, data$SPiPscore, data$Interpretation))
    data$Interpretation[data$SPiPscore<thToComplexEvent & data$Interpretation=="Alter by complex event"] = "NTR"

    oldNames = names(data)[-which(names(data)%in%c("varID", "Interpretation", "InterConfident", "SPiPscore", "chr", "strand", "gNomen", "varType", "ntChange", "ExonInfo", "exonSize", "transcript",
    "gene", "NearestSS", "DistSS", "RegType", "seqPhysio", "seqMutated", "SPiCEproba", "SPiCEinter_2thr", "deltaMES", "BP", "mutInPBarea", "deltaESRscore",
    "posCryptMut", "sstypeCryptMut", "probaCryptMut", "classProbaCryptMut", "nearestSStoCrypt", "nearestPosSStoCrypt", "nearestDistSStoCrypt", "posCryptWT",
    "probaCryptWT", "classProbaCryptWT", "posSSPhysio", "probaSSPhysio", "classProbaSSPhysio", "probaSSPhysioMut", "classProbaSSPhysioMut","RegTypeNum"))]

    data = data[,c(oldNames,"varID", "Interpretation", "InterConfident", "SPiPscore", "chr", "strand", "gNomen", "varType", "ntChange", "ExonInfo", "exonSize", "transcript",
    "gene", "NearestSS", "DistSS", "RegType", "seqPhysio", "seqMutated", "SPiCEproba", "SPiCEinter_2thr", "deltaMES", "BP", "mutInPBarea", "deltaESRscore",
    "posCryptMut", "sstypeCryptMut", "probaCryptMut", "classProbaCryptMut", "nearestSStoCrypt", "nearestPosSStoCrypt", "nearestDistSStoCrypt", "posCryptWT",
    "probaCryptWT", "classProbaCryptWT", "posSSPhysio", "probaSSPhysio", "classProbaSSPhysio", "probaSSPhysioMut", "classProbaSSPhysioMut")]

    return(data)
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


mergeSPiPresult<-function(x){
    c = paste0("SPiP=",paste0(INFO[which(VCFinfo_text%in%x)],collapse=";"))
    return(c)
}

convertLine2VCF <- function(line){
    ID = as.character(line$varID)
    QUAL = "."
    FILTER = "."
    if(line$RegType!="."){
        seq = line$seqPhysio
        strand = line$strand
        varType = line$varType
        ntChange = line$ntChange
        if(strand=="-"){seq = getRevSeq(seq)}

        CHROM = line$chr
        if(varType=="substitution"){
            POS = line$gNomen
            REF = if(strand=="+"){unlist(strsplit(ntChange,">",fixed=T))[1]}else if(strand=="-"){getRevSeq(unlist(strsplit(ntChange,">",fixed=T))[1])}
            ALT = if(strand=="+"){unlist(strsplit(ntChange,">",fixed=T))[2]}else if(strand=="-"){getRevSeq(unlist(strsplit(ntChange,">",fixed=T))[2])}
        }else if(varType=="ins"){
            POS = as.numeric(line$gNomen)-1
            REF = substr(seq,151,151)
            ALT = paste0(REF,if(strand=="+"){gsub("ins","",ntChange)}else if(strand=="-"){getRevSeq(gsub("ins","",ntChange))})
        }else if(varType=="dup"){
            varPOS = as.numeric(unlist(strsplit(line$gNomen,"_",fixed=T)))
            if(length(varPOS)==1){varPOS=rep(varPOS,2)}
            POS = min(varPOS)-1
            mutSize = abs(varPOS[1]-varPOS[2])
            ntDup = substr(seq,151,151+mutSize)
            REF = paste0(substr(seq,150,150),ntDup)
            ALT = paste0(substr(seq,150,150),paste(rep(ntDup,2),collapse=""))
        }else if(varType=="del"){
            varPOS = as.numeric(unlist(strsplit(line$gNomen,"_",fixed=T)))
            if(length(varPOS)==1){varPOS=rep(varPOS,2)}
            POS = min(varPOS)-1
            mutSize = abs(varPOS[1]-varPOS[2])
            ntDel = substr(seq,151,151+mutSize)
            REF = paste0(substr(seq,150,150),ntDel)
            ALT = substr(seq,150,150)
        }else if(varType=="delins"){
            varPOS = as.numeric(unlist(strsplit(line$gNomen,"_",fixed=T)))
            if(length(varPOS)==1){varPOS=rep(varPOS,2)}
            POS = min(varPOS)-1
            mutSize = abs(varPOS[1]-varPOS[2])
            ntDel = substr(seq,151,151+mutSize)
            ntIns = if(strand=="+"){gsub("delins","",ntChange)}else if(strand=="-"){getRevSeq(gsub("delins","",ntChange))}
            REF = paste0(substr(seq,150,150),ntDel)
            ALT = paste0(substr(seq,150,150),ntIns)
        }

        INFO = paste(line[c("Interpretation", "InterConfident", "SPiPscore", "strand", "gNomen", "varType", "ntChange", "ExonInfo", "exonSize", "transcript",
                "gene", "NearestSS", "DistSS", "RegType", "SPiCEproba", "SPiCEinter_2thr", "deltaMES", "BP", "mutInPBarea",
                "deltaESRscore", "posCryptMut", "sstypeCryptMut", "probaCryptMut", "classProbaCryptMut", "nearestSStoCrypt", "nearestPosSStoCrypt",
                "nearestDistSStoCrypt", "posCryptWT", "probaCryptWT", "classProbaCryptWT", "posSSPhysio", "probaSSPhysio", "classProbaSSPhysio",
                "probaSSPhysioMut", "classProbaSSPhysioMut")],collapse="|")
    }else{
        CHROM <- POS <- REF <- ALT <- INFO  <- "."
    }

    VCFline = c(CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO)
    return(VCFline)
}

#import data
readVCF <- function(txtLine,i){
    if(printProcess){setTxtProgressBar(pb, i)}
    dataLine = unlist(strsplit(txtLine,split='\t')) #c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO",SUITE?)
    chrom <- dataLine[1]
    pos <- as.numeric(dataLine[2])
    ref <- dataLine[4]
    alt <- dataLine[5]
	lengthDataLine <- length(dataLine)

    if(ref=="."|alt=="."){
		resultMatrix <- c("mutUnknown","","","","")
		for (k in 1:lengthDataLine) {
			resultMatrix <- c(resultMatrix,dataLine[k])
		}
    }else{
        if(substr(as.character(chrom[1]),1,3)!="chr"){
            if(substr(as.character(chrom[1]),1,3)=="NC_"){
                chrom = contigToChr(as.character(chrom))
            }else{
                chrom = paste("chr",chrom,sep="")
            }
        }

		bonus <- nchar(ref) # in case ref is deleted, we need some additional bases to compute the deleted seq

        transcriptLines <- dataRefSeq[dataRefSeq$V1==chrom & dataRefSeq$V2<=pos & dataRefSeq$V3>=pos,]
		transcript <- transcriptLines[,'V4']
        strandTrans <- transcriptLines[,'V6']

        if(length(transcript)==0){
            resultMatrix <- matrix(c(paste("no transcript",pos,sep=":"),"","","","",dataLine),ncol = 5+lengthDataLine)
        }else{
            totalSequence <- getSequencePhysio(genome, NULL, chrom, max(1,pos - 150 - bonus), pos+150+bonus)
            if( length(unlist(strsplit(totalSequence,"|")))>0){
                seqStart <- max(1,pos-150-bonus)
                seqLength <- nchar(totalSequence)
                if (ref == "-") {
                    posIndex <- c(pos-seqStart+1,pos-seqStart) # special case for comment below : ref == - => last = first-1
                } else {
                    posIndex <- c(pos-seqStart+1,pos-seqStart+nchar(ref)) # posIndex[1] : first index of ref / posIndex[2] : last index of ref
                }

                # computed once for all '+' transcripts
                seqPhysio <- substr(totalSequence, posIndex[1]-150, posIndex[1]+150)
                seqMutated1 <- substr(totalSequence, posIndex[1]-150, posIndex[1]-1)
                seqMutated2 <- substr(totalSequence, posIndex[2]+1, posIndex[2]+150)

                # computed once for all '-' transcripts
                revSeqPhysio <- getRevSeq(substr(totalSequence, posIndex[2]-150, posIndex[2]+150))
                revSeqMutated1 <- getRevSeq(substr(totalSequence, posIndex[2]+1, posIndex[2]+150))
                revSeqMutated2 <- getRevSeq(substr(totalSequence, posIndex[1]-150, posIndex[1]-1))
                revRef <- getRevSeq(ref)

                altSplit = unlist(strsplit(alt,',',fixed = TRUE))
                nbAlt <- length(altSplit)

                resultMatrix <- matrix(rep("",length(transcript)*nbAlt*(5+lengthDataLine)), ncol = 5+lengthDataLine)
                # matrix dedicated to the storage of [varID, seqPhysio, seqMutated, transcript, altUsed, +infosOfInput] for each variant

                for (i in 1:length(transcript)) {
                    if(strandTrans[i] =="+"){
                        for(j in 1:nbAlt){
                            resultMatrix[(i-1)*nbAlt + j,2] <- seqPhysio
                            resultMatrix[(i-1)*nbAlt + j,4] <- as.character(transcript[i])
                            resultMatrix[(i-1)*nbAlt + j,5] <- altSplit[j]
                            for (k in 1:lengthDataLine) {
                                resultMatrix[(i-1)*nbAlt + j,5+k] <- dataLine[k]
                            }
                            if(nchar(ref)==1 & nchar(altSplit[j])==1){
                                resultMatrix[(i-1)*nbAlt + j,1] <- paste(transcript[i],':g.',pos,':',ref,'>',altSplit[j],sep="")
                                resultMatrix[(i-1)*nbAlt + j,3] <- paste(seqMutated1,altSplit[j],seqMutated2,sep="")
                            }else{ # all other variants than subsitution are considered as delins
                                resultMatrix[(i-1)*nbAlt + j,1] <- paste(transcript[i],':g.',pos,"_",pos+nchar(ref)-1,':delins',altSplit[j],sep="")
                                resultMatrix[(i-1)*nbAlt + j,3] <- substr(paste(seqMutated1,altSplit[j],seqMutated2,sep=""),1,301)
                            }
                        }
                    } else if (strandTrans[i]=="-"){
                        for(j in 1:nbAlt){
                            altRev = getRevSeq(altSplit[j])
                            resultMatrix[(i-1)*nbAlt + j,2] <- revSeqPhysio
                            resultMatrix[(i-1)*nbAlt + j,4] <- as.character(transcript[i])
                            resultMatrix[(i-1)*nbAlt + j,5] <- altSplit[j] # not reverse-complemented
                            for (k in 1:lengthDataLine) {
                                resultMatrix[(i-1)*nbAlt + j,5+k] <- dataLine[k]
                            }
                            if(nchar(revRef)==1 & nchar(altRev)==1){
                                resultMatrix[(i-1)*nbAlt + j,1] <- paste(transcript[i],':g.',pos,':',revRef,'>',altRev,sep="")
                                resultMatrix[(i-1)*nbAlt + j,3] <- paste(revSeqMutated1,altRev,revSeqMutated2,sep="")
                            }else{
                                resultMatrix[(i-1)*nbAlt + j,1] <- paste(transcript[i],':g.',pos+nchar(revRef)-1,"_",pos,':delins',altRev,sep="")
                                resultMatrix[(i-1)*nbAlt + j,3] <- substr(paste(revSeqMutated1,altRev ,revSeqMutated2,sep=""),1,301)
                            }
                        }
                    }
                }
            }else{
                resultMatrix <- matrix(c(paste("no sequence",pos,sep=":"),"","","","",dataLine),ncol = 5+lengthDataLine)
            }

        }
    }
    result <<- unlist(as.list(t(resultMatrix)))
}
