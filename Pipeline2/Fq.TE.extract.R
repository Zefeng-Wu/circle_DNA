suppressMessages(require(ShortRead))
suppressMessages(require(Biostrings))
suppressMessages(require(stringr))
suppressMessages(library(GenomicAlignments))

args<-commandArgs(T)

all_reads_file <-args[1]
all_reads<-readFastq(all_reads_file) #"/home/wuzefeng/Myresearch/ZZ/2020.3.28/zz.fq"
all_fq_name<-as.character(subseq(all_reads@id,1,36))

TEs_align_reads_file <- args[2]
TEs_align_reads <- readGAlignments(TEs_align_reads_file,use.names = TRUE)

TE_align_reads_id <-unique(names(TEs_align_reads))
TE_aligned_reads_fq <- all_reads[match(TE_align_reads_id,all_fq_name)]
non_TE_align_reads_fq <- all_reads[-match(TE_align_reads_id,all_fq_name)]

outname_prefix1 <- paste(str_sub(basename(all_reads_file),start = 1,end = -3),"TE_reads.fa",sep = "")
outname_prefix2 <- paste(str_sub(basename(all_reads_file),start = 1,end = -3),"non_TE_reads.fa",sep = "") 

writeFasta(object = TE_aligned_reads_fq, file = outname_prefix1)
writeFasta(object = non_TE_align_reads_fq, file = paste("../2non_TE_reads/",outname_prefix2,sep=""))

