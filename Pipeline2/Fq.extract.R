suppressMessages(require(ShortRead))
suppressMessages(require(Biostrings))
suppressMessages(require(stringr))

args<-commandArgs(T)

all_reads<-readFastq(args[1]) #"/home/wuzefeng/Myresearch/ZZ/2020.3.28/zz.fq"
all_fq_name<-as.character(subseq(all_reads@id,1,36))

repeat_reads_cons.fa_file<-args[2]
repeat_reads_cons.fa <- readDNAStringSet(repeat_reads_cons.fa_file,format = "fasta") #"/home/wuzefeng/Myresearch/ZZ/2020.3.28/1_zz.cons.fa"
repeat_reads_ids<- unique(str_split_fixed(names(repeat_reads_cons.fa),pattern = "_",n = 2)[,1])

repeat_reads_fq <- all_reads[match(repeat_reads_ids,all_fq_name)]
non_repeat_reads_fq <- all_reads[-match(repeat_reads_ids,all_fq_name)]

outname_prefix <- paste("1_",str_sub(basename(repeat_reads_cons.fa_file),start = 3,end = -8),"original_reads.fa",sep = "") 
writeFasta(object = repeat_reads_fq, file = paste("1repeated_reads/",outname_prefix,sep=""))        # 1_zz.fa
writeFasta(object = non_repeat_reads_fq,file = paste("1non-repeated_reads/",outname_prefix,sep="")) # 1_zz.fa

