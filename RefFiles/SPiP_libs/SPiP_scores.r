checkCalculableScore <- function(seq){
    scoreCalculable<<-TRUE
    if(as.numeric(regexpr("N",seq))>0){
        print(paste("Unknown sequence:",seq))
        scoreCalculable<<-FALSE
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
    			SSF = SSFdonGC(seqCon)
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
