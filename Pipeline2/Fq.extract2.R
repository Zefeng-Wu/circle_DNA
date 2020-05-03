suppressMessages(require(ShortRead))
suppressMessages(require(Biostrings))
suppressMessages(require(stringr))
suppressMessages(library(GenomicAlignments))

args<-commandArgs(T)

all_reads<-readFastq(args[1]) #"/home/wuzefeng/Myresearch/ZZ/2020.3.28/zz.fq"
all_fq_name<-as.character(subseq(all_reads@id,1,36))

repeat_reads_cons.fa_file <- args[2]
repeat_reads_cons.fa <- readDNAStringSet(repeat_reads_cons.fa_file,format = "fasta") #"/home/wuzefeng/Myresearch/ZZ/2020.3.28/1_zz.cons.fa"
names(repeat_reads_cons.fa) <- str_split_fixed(names(repeat_reads_cons.fa),pattern = " ",n = 2)[,1]
repeat_reads_ids<- unique(str_split_fixed(names(repeat_reads_cons.fa),pattern = "_",n = 2)[,1])

align_reads_file <- args[3]
message(align_reads_file)
align_reads <- readGAlignments(align_reads_file,use.names = TRUE)
align_reads_id <- unique(str_split_fixed(names(align_reads),pattern = "_",n = 2)[,1])

aligned_reads_fq <- all_reads[match(align_reads_id,all_fq_name)]
non_align_reads_fq <- all_reads[match(repeat_reads_ids[!repeat_reads_ids%in%align_reads_id],all_fq_name)]


outname_prefix1 <- paste("1_",str_sub(basename(repeat_reads_cons.fa_file),start = 3,end = -8),"mapped.fa",sep = "")
outname_prefix2 <- paste("1_",str_sub(basename(repeat_reads_cons.fa_file),start = 3,end = -8),"unmapped.fa",sep = "") 


writeFasta(object = aligned_reads_fq, file = paste("mapping/",outname_prefix1,sep=""))
writeFasta(object = non_align_reads_fq, file = paste("mapping/",outname_prefix2,sep=""))

