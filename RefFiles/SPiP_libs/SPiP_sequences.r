getRevSeq <- function(sequence){
    splitRevSeq = rev(unlist(strsplit(sequence,'|')))
    seqRev = paste(unlist(lapply(list(splitRevSeq),function(x) inverseDic[x,1])),sep="",collapse="")
    return(seqRev)
}

getSequencePhysio <- function(genome,sens=NULL,chr,start,end){
	if(is.null(sens)){
        idx = which(transcriptome_idx$chr==chr & transcriptome_idx$seqStart<=start & transcriptome_idx$seqEnd>=end)
    } else{
        idx = which(transcriptome_idx$chr==chr & transcriptome_idx$strand==sens & transcriptome_idx$seqStart<=start & transcriptome_idx$seqEnd>=end)
    }
	if(length(idx)>1){idx = idx[1]}
    if(length(idx)>0){
        start2 = start;end2 = end
        seqStart = transcriptome_idx$seqStart[idx]-1
        if(transcriptome_idx$strand[idx]=="-"){
            start2=end;end2=start
            seqStart = transcriptome_idx$seqEnd[idx]+1
        }
        transcript_sequence <- transcriptome_seq[idx,2]
        seqDNA = substr(transcript_sequence,abs(start2-seqStart),abs(end2-seqStart))
        if(is.null(sens) & transcriptome_idx$strand[idx]=="-"){seqDNA = getRevSeq(seqDNA)} # by default all sequence should be forward strand
    }else{
        seqDNA=""
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
			seqMutated2=getSequencePhysio(genome,sens,chr,start,end)
			seqMutated <<- paste(seqMutated1 , seqMutated2,sep="")
		}else if(sens=="-"){
			seqMutated1 = substr(seqPhysio,1,150)
			start=varPos[2]-1
			end=varPos[2]-151
			seqMutated2=getSequencePhysio(genome,sens,chr,end,start)
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
			seqMutated2 = getSequencePhysio(genome,sens,chr,start,end)
			seqMutated <<- substr(paste(seqMutated1,ntIns ,seqMutated2,sep=""),1,301)

		}else if(sens=="-"){

			seqMutated1 = substr(seqPhysio,1,150)
			start = varPos[2]-1
			end = varPos[2]-151
			seqMutated2 = getSequencePhysio(genome,sens,chr,end,start)
			seqMutated <<- substr(paste(seqMutated1,ntIns ,seqMutated2,sep=""),1,301)

		}
	}
    return(seqMutated)
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
