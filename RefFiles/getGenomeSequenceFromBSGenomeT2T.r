library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hs1")

genome <- getBSgenome("BSgenome.Hsapiens.UCSC.hs1")

chr1 = as.character(genome$chr1)
chr2 = as.character(genome$chr2)
chr3 = as.character(genome$chr3)
chr4 = as.character(genome$chr4)
chr5 = as.character(genome$chr5)
chr6 = as.character(genome$chr6)
chr7 = as.character(genome$chr7)
chr8 = as.character(genome$chr8)
chr9 = as.character(genome$chr9)
chr10 = as.character(genome$chr10)
chr11 = as.character(genome$chr11)
chr12 = as.character(genome$chr12)
chr13 = as.character(genome$chr13)
chr14 = as.character(genome$chr14)
chr15 = as.character(genome$chr15)
chr16 = as.character(genome$chr16)
chr17 = as.character(genome$chr17)
chr18 = as.character(genome$chr18)
chr19 = as.character(genome$chr19)
chr20 = as.character(genome$chr20)
chr21 = as.character(genome$chr21)
chr22 = as.character(genome$chr22)
chrX = as.character(genome$chrX)
chrY = as.character(genome$chrY)

## get transcriptome
load("dataRefSeqT2T.RData")
cat("loaded refseq")

seqStart = NULL
seqEnd = NULL
chr = NULL
Gene = NULL
strand = NULL

for(gene in unique(dataRefSeq$V13)){
	seqStart = c(seqStart,min(dataRefSeq$V2[dataRefSeq$V13==gene])-150)
	seqEnd = c(seqEnd,max(dataRefSeq$V3[dataRefSeq$V13==gene])+150)
	chr = c(chr,unique(as.character(dataRefSeq$V1[dataRefSeq$V13==gene])))
	Gene = c(Gene,gene)
	strand = c(strand,unique(as.character(dataRefSeq$V6[dataRefSeq$V13==gene])))
}

transcriptome = data.frame(chr,seqStart,seqEnd,strand,Gene)
transcriptome = transcriptome[transcriptome$chr!="chrM",]

inverseDic <- data.frame(V1 = c('T','G','C','A','N'),row.names = c('A','C','G','T','N'))
inverseDic$V1 <- as.character(inverseDic$V1)

getRevSeq <- function(sequence){
    splitRevSeq = rev(unlist(strsplit(sequence,'|')))
    seqRev = paste(unlist(lapply(list(splitRevSeq),function(x) inverseDic[x,1])),sep="",collapse="")
    return(seqRev)
}

getSeq <- function(chr,start,end){
	eval(parse(text = paste("sequence = substr(",chr,",",start,",",end,")",sep="")))
	return(sequence)
}
tmpSeq = unlist(mapply(getSeq,transcriptome$chr,transcriptome$seqStart,transcriptome$seqEnd))

transcriptome$seq = tmpSeq

revSeq = unlist(lapply(transcriptome$seq,getRevSeq))

transcriptome$seq[transcriptome$strand=="-"] = revSeq[transcriptome$strand=="-"]

transcriptome_idx = transcriptome[,c(1:4)]
transcriptome_seq = transcriptome[,c(5:6)]
cat("finished filling transcriptome")

save(transcriptome_idx, transcriptome_seq,file = "transcriptome_T2T.RData")
cat("saved transcriptome")
