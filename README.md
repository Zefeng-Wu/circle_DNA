## This pipeline is used to identify potential circle DNA from nanopore sequencing technology 

## method one: TRF
mkdir pass_merge
unzip pass.zip
for m in $(ls *.gz); do gunzip $m ; done
cat *.fastq > ../pass_merge/all.fastq

### reads quality statistics
NanoPlot -t 40 --fastq  all.fastq --plots hex dot

## fastq convert to fasta
sed -n '1~4s/^@/>/p;2~4p' all.fastq > all.fasta

## perform tandem repeat finder to seach TRs in reads
~/soft/trf409.linux64 all.fasta 2 7 7 80 10 50 2000 -l 6 -h > out.txt
# match score: 2
# mismatch score: 7
# indel score: 7
# acth probality: 80 
# min alignemt score to report: 50
# max period: 2000
# -l longest repeat size(million) 

### format convert
python test.py in.dat out.txt



### method 2: TideHunter
TideHunter all.fastq -t 40  -o 4TideHunter_result/cons.fa

## minimap2
minimap2 -ax map-ont ../reference_genome/Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa  4TideHunter_result/cons.fa -t 30 --sam-hit-only --secondary=no -o 5minimap2_result/minimap2.sam

## sam to bam
samtools view -S -b minimap2.sam > minimap2.bam

## obtain aligned reads
samtools view -bF 4 minimap2.bam > 3minimap2_aligned.bam

## sort sam
java -jar ~/soft/picard.jar  SortSam I= 3minimap2_aligned.bam O= 4minimap2_aligned_sorted.bam SORT_ORDER=coordinate

## load into R 
library(genomicFeatures)
library(bamsignals)

tr<-makeTxDbFromGFF("../../reference_genome/Drosophila_melanogaster.BDGP6.22.97.gtf")
genes<-genes(tr)

covSigs <- bamCount("4minimap2_aligned_sorted.bam", genes)
genes$reads_count <-covSigs
write.table(as.data.frame(genes),"out",sep="\t",col.names=TRUE,row.names=TRUE,quot=FALSE)

