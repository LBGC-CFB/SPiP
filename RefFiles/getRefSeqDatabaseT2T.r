########
#get RefSeq dataBase from UCSC
########
options(stringsAsFactors=FALSE)
message("Use the URL: http://hgdownload.cse.ucsc.edu/goldenPath/")
argsFull <- commandArgs()
Rscript <- argsFull[1]

inputref=dirname(normalizePath(sub("--file=","",argsFull[substr(argsFull,1,7)=="--file="])))

#old format vs new format
# |  1 | bin          |  1 | chrom        |
# |  2 | name         |  2 | chromStart   |
# |  3 | chrom        |  3 | chromEnd     |
# |  4 | strand       |  4 | name         |
# |  5 | txStart      |  5 | score        |
# |  6 | txEnd        |  6 | strand       |
# |  7 | cdsStart     |  7 | thickStart   |
# |  8 | cdsEnd       |  8 | thickEnd     |
# |  9 | exonCount    |  9 | reserved     |
# | 10 | exonStarts   | 10 | blockCount   |
# | 11 | exonEnds     | 11 | blockSizes   |
# | 12 | score        | 12 | chromStarts  |
# | 13 | name2        | 13 | name2        |
# | 14 | cdsStartStat | 14 | cdsStartStat |
# | 15 | cdsEndStat   | 15 | cdsEndStat   |
# | 16 | exonFrames   | 16 | exonFrames   |
# |    |              | 17 | type         |
# |    |              | 18 | geneName     |
# |    |              | 19 | geneName2    |
# |    |              | 20 | geneType     |

#open connection
## message("open connexion...")
## system("curl -O https://hgdownload.soe.ucsc.edu/gbdb/hs1/ncbiRefSeq/ncbiRefSeqCurated.bb")
## Requires http://hgdownload.soe.ucsc.edu/admin/exe/biobig
## Don't write header as it misses a tab somewhere and break the file reading
system("./bigBedToBed ncbiRefSeqCurated.bb ncbiRefSeqCurated.bed")

message("read files...")
dat <- read.csv("ncbiRefSeqCurated.bed", sep="\t", header=FALSE)
dat$V7 = as.numeric(dat$V7)

# No need to rename chromosome in T2T

#remove duplicate transcrit between chrX and chrY
dat = dat[order(dat$V1),]
dat = dat[!duplicated(dat$V4),]

#remove transcrit version
splitTrans = unlist(strsplit(dat$V4,".",fixed=TRUE))
dat$V4 = splitTrans[seq(from=1,to=length(splitTrans),by=2)]

# Exon relative position and length is already stored
## dathg19$V17 = unlist(tmpResulthg19[1,])
## dathg19$V18 = unlist(tmpResulthg19[2,])

#save data
## message("Save data...")

dataRefSeq = dat[,c(
  "V1", # chrom
  "V2", # txStart
  "V3", # txEnd
  "V4", # name
  "V5", # score
  "V6", # strand
  "V7", # cdsstart ?
  "V8", # cdsend ?
  "V5", # score again
  "V10", # exoncount
  "V11", # exon length
  "V12", # exon relative pos
  "V13" # name3
  )]
## names(dataRefSeq) = paste("V",c(1:13),sep="")

save(dataRefSeq,file = paste(inputref,"dataRefSeqT2T.RData",sep=""))
