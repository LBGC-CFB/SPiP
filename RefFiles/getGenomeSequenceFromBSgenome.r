# use Rv4.0.3

library("BSgenome")
install.packages("BiocManager")

BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library("BSgenome.Hsapiens.UCSC.hg38")

genome <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38")

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

save(chr1,file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr1.RData")
save(chr2,file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr2.RData")
save(chr3,file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr3.RData")
save(chr4,file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr4.RData")
save(chr5,file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr5.RData")
save(chr6,file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr6.RData")
save(chr7,file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr7.RData")
save(chr8,file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr8.RData")
save(chr9,file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr9.RData")
save(chr10,file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr10.RData")
save(chr11,file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr11.RData")
save(chr12,file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr12.RData")
save(chr13,file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr13.RData")
save(chr14,file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr14.RData")
save(chr15,file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr15.RData")
save(chr16,file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr16.RData")
save(chr17,file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr17.RData")
save(chr18,file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr18.RData")
save(chr19,file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr19.RData")
save(chr20,file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr20.RData")
save(chr21,file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr21.RData")
save(chr22,file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr22.RData")
save(chrX,file = "C:/Users/LEMRAP/Desktop/SPP_tool/chrX.RData")
save(chrY,file = "C:/Users/LEMRAP/Desktop/SPP_tool/chrY.RData")

q(save="no")

load(file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr1.RData")
load(file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr2.RData")
load(file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr3.RData")
load(file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr4.RData")
load(file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr5.RData")
load(file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr6.RData")
load(file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr7.RData")
load(file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr8.RData")
load(file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr9.RData")
load(file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr10.RData")
load(file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr11.RData")
load(file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr12.RData")
load(file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr13.RData")
load(file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr14.RData")
load(file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr15.RData")
load(file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr16.RData")
load(file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr17.RData")
load(file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr18.RData")
load(file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr19.RData")
load(file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr20.RData")
load(file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr21.RData")
load(file = "C:/Users/LEMRAP/Desktop/SPP_tool/chr22.RData")
load(file = "C:/Users/LEMRAP/Desktop/SPP_tool/chrX.RData")
load(file = "C:/Users/LEMRAP/Desktop/SPP_tool/chrY.RData")

save.image(file = "C:/Users/LEMRAP/Desktop/SPP_tool/hg38.RData")

# get transcriptome

load("C:/Users/LEMRAP/Desktop/SPP_tool/hg38.RData")
load("C:/Users/LEMRAP/Desktop/SPP_tool/R-4.0.2/RefFiles/dataRefSeqhg38.RData")

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

save(transcriptome_idx, transcriptome_seq,file = "C:/Users/LEMRAP/Desktop/SPP_tool/transcriptome_hg38.RData")
