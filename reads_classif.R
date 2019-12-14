library(ShortRead)
library(Biostrings)
library(GenomicAlignments)
library(stringr)
all_fq<-readFastq("/home/wuzefeng/opt/zz/191018_filter.fq") #285536
with_repeat_seqs<-names(readDNAStringSet("/home/wuzefeng/opt/zz/191018.cons.fa"))
with_repeat_seqs<-unique(str_split_fixed(with_repeat_seqs,pattern = "_",n = 2)[,1]) # 41499

all_fq_name<-as.character(subseq(all_fq@id,1,36))
fq_without_repeat<-all_fq[-match(with_repeat_seqs,all_fq_name)] # 244037
writeFasta(fq_without_repeat,file = "191018_without_repeat.fa")

boxplot(width(all_fq[match(with_repeat_seqs,all_fq_name)]),width(all_fq[-match(with_repeat_seqs,all_fq_name)]),
        names=c("Orignal reads with repeat","Orignal reads without repeat"),ylab="Reads length",col="steelblue")

## inspect mapped reads length 
mapped_reads<-readGAlignments("/home/wuzefeng/opt/zz/191018.attB.aligned.sorted.bam",use.names=TRUE)
plot(log(qwidth(mapped_reads),2),log(width(mapped_reads),2),xlab="Read length (log2)",ylab="Mapping length (log2)")
hist(qwidth(mapped_reads),xlab = "Mapped reads length",col="steelblue",main="Mapped reads length")

mappd_reads_name<-unique(str_split_fixed(names(mapped_reads),pattern = "_",n=2)[,1]) #28560
unmapped_read_name <-with_repeat_seqs[!with_repeat_seqs%in%mappd_reads_name] #12939
unmapped_reads <-all_fq[match(unmapped_read_name,all_fq_name)] #12939
boxplot(qwidth(mapped_reads),width(unmapped_reads),names=c("Mapped reads with repeat","Unmapped reads with repeat"),ylab="Reads length",col="steelblue")
