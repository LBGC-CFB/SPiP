########
#get RefSeq dataBase from UCSC
########
options(stringsAsFactors=FALSE)
message("Use the URL: http://hgdownload.cse.ucsc.edu/goldenPath/")
wd_R=getwd()
inputref = paste(wd_R, "/RefFiles",sep="")
#format of refseq: bin, transcrit, chr, strand, gstart, gend, CDSstart, CDSend, nEx, posExStart, posCum, Score (0), geneName, Annot1, Annot2, frame

#open connection
message("open connexion...")
hg38con <- gzcon(url("http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/ncbiRefSeqCurated.txt.gz"))
hg19con <- gzcon(url("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/ncbiRefSeqCurated.txt.gz"))

#uncompress
message("parse files...")
lineshg38 <- readLines(hg38con)
message("parse files...")
lineshg19 <- readLines(hg19con)

close(hg38con)
close(hg19con)

#read file
splitRawToTable <- function(raw, sep = "\t", head = FALSE){
    if(head){
        columNames = unlist(strsplit(raw[1],sep,fixed=TRUE))
        nCol = length(columNames)
        splitRaw = unlist(strsplit(raw[2:length(raw)],sep,fixed=TRUE))
        data=as.data.frame(matrix(splitRaw, ncol = nCol,byrow = TRUE))
        colnames(data) <- columNames
    }else{
        nCol = length(as.numeric(unlist(gregexpr("\t", raw[1] ,fixed=TRUE))))+1
        splitRaw = unlist(strsplit(raw, sep, fixed=TRUE))
        data=as.data.frame(matrix(splitRaw, ncol = nCol, byrow = TRUE))
    }
    return(data)
}

message("read files...")
dathg38 <- splitRawToTable(lineshg38)
dathg38$V5 = as.numeric(dathg38$V5)

dathg19 <- splitRawToTable(lineshg19)
dathg19$V5 = as.numeric(dathg19$V5)

#remove irrelevant chr
message("Data formatting...")
dathg38 = dathg38[-grep("_",dathg38$V3),]
dathg19 = dathg19[-grep("_",dathg19$V3),]

#remove duplicate transcrit between chrX and chrY
dathg38 = dathg38[order(dathg38$V3),]
dathg19 = dathg19[order(dathg19$V3),]

dathg38 = dathg38[!duplicated(dathg38$V2),]
dathg19 = dathg19[!duplicated(dathg19$V2),]

#remove transcrit version
splitTranshg38 = unlist(strsplit(dathg38$V2,".",fixed=TRUE))
splitTranshg19 = unlist(strsplit(dathg19$V2,".",fixed=TRUE))
dathg38$V2 = splitTranshg38[seq(from=1,to=length(splitTranshg38),by=2)]
dathg19$V2 = splitTranshg19[seq(from=1,to=length(splitTranshg19),by=2)]

getExonInfo <- function(start, listStart, listEnd){
	splitStart = unlist(strsplit(listStart,","))
	splitEnd = unlist(strsplit(listEnd,","))
	splitStart = as.numeric(splitStart)
	splitEnd = as.numeric(splitEnd)
	lenEx = splitEnd - splitStart
	relPos = splitStart - start
	result <<- list(paste(paste(lenEx,collapse=","),",",sep=""), paste(paste(relPos,collapse=","),",",sep=""))
}

tmpResulthg38 = mapply(getExonInfo,dathg38$V5,dathg38$V10,dathg38$V11)
tmpResulthg19 = mapply(getExonInfo,dathg19$V5,dathg19$V10,dathg19$V11)

dathg38$V17 = unlist(tmpResulthg38[1,])
dathg38$V18 = unlist(tmpResulthg38[2,])
dathg19$V17 = unlist(tmpResulthg19[1,])
dathg19$V18 = unlist(tmpResulthg19[2,])

#save data
message("Save data...")
dataRefSeq = dathg38[,c("V3","V5","V6","V2","V12","V4","V7","V8","V12","V9","V17","V18")]
names(dataRefSeq) = paste("V",c(1:12),sep="")
dataRefSeq$V1 = as.factor(dataRefSeq$V1)
dataRefSeq$V3 = as.numeric(dataRefSeq$V3)
dataRefSeq$V4 = as.factor(dataRefSeq$V4)
dataRefSeq$V5 = as.numeric(dataRefSeq$V5)
dataRefSeq$V6 = as.factor(dataRefSeq$V6)
dataRefSeq$V7 = as.numeric(dataRefSeq$V7)
dataRefSeq$V8 = as.numeric(dataRefSeq$V8)
dataRefSeq$V9 = as.numeric(dataRefSeq$V9)
dataRefSeq$V10 = as.numeric(dataRefSeq$V10)
dataRefSeq$V11 = as.factor(dataRefSeq$V11)
dataRefSeq$V12 = as.factor(dataRefSeq$V12)
#dataRefSeq$V13 = as.factor(dataRefSeq$V13)
save(dataRefSeq,file = paste(inputref,"/dataRefSeqhg38.RData",sep=""))

dataRefSeq = dathg19[,c("V3","V5","V6","V2","V12","V4","V7","V8","V12","V9","V17","V18")]
names(dataRefSeq) = paste("V",c(1:12),sep="")
dataRefSeq$V1 = as.factor(dataRefSeq$V1)
dataRefSeq$V3 = as.numeric(dataRefSeq$V3)
dataRefSeq$V4 = as.factor(dataRefSeq$V4)
dataRefSeq$V5 = as.numeric(dataRefSeq$V5)
dataRefSeq$V6 = as.factor(dataRefSeq$V6)
dataRefSeq$V7 = as.numeric(dataRefSeq$V7)
dataRefSeq$V8 = as.numeric(dataRefSeq$V8)
dataRefSeq$V9 = as.numeric(dataRefSeq$V9)
dataRefSeq$V10 = as.numeric(dataRefSeq$V10)
dataRefSeq$V11 = as.factor(dataRefSeq$V11)
dataRefSeq$V12 = as.factor(dataRefSeq$V12)
#dataRefSeq$V13 = as.factor(dataRefSeq$V13)
save(dataRefSeq,file = paste(inputref,"/dataRefSeqhg19.RData",sep=""))
