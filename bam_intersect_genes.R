library(GenomicAlignments)
library(GenomicRanges)
library(GenomicFeatures)
aln1 <- readGAlignments("~/Myresearch/ZZ/2020.2.7/200117/4minimap2.sorted.bam",use.names = TRUE)
reads_names<-str_split_fixed(names(aln1),pattern  = "_",n =3)[,1] # 758818 hits
reads_num <-length(unique(reads_names)) #249188 reads

### load genes
tr<-makeTxDbFromGFF("~/Myresearch/Genomes/GTF/Drosophila_melanogaster.BDGP6.28.99.gtf")
genes <- genes(tr)
reads_op_genes_hits <- findOverlaps(aln1, genes, ignore.strand=TRUE)
reads_op_genes<-aln1[unique(queryHits(reads_op_genes_hits))]
reads_number<-length(unique(str_split_fixed(names(reads_op_genes),pattern  = "_",n =3)[,1])) # 245583 (normal chromsomes:5952)

### mitogenes
genes <-genes[seqnames(genes)=="mitochondrion_genome"]
reads_op_genes_hits <- findOverlaps(aln1, genes, ignore.strand=TRUE)
reads_op_genes<-aln1[unique(queryHits(reads_op_genes_hits))]
reads_number<-length(unique(str_split_fixed(names(reads_op_genes),pattern  = "_",n =3)[,1])) # 239783


### equal to in shell
export.bed(genes,"genes.bed")
system("samtools view select.bam  | cut -f1 | awk '{u = split($1,m,"_"); print m[1]}' | sort | uniq | wc -l")
