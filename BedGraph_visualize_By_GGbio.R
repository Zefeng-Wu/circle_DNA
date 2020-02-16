library(rtracklayer)
library(ggbio)
library(Biostrings)
library(stringr)

## obatin the chromsome length from genome fasta file
dna <- readDNAStringSet("~/Myresearch/Genomes/DNA/Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa")
names(dna) <- str_split_fixed(names(dna),pattern = " ",n = 2)[,1]
dna <-dna[names(dna)%in%c("2L","2R","3L","3R","4","X","Y","mitochondrion_genome")]
gr <- GRanges(seqnames = names(dna),IRanges(start = 0, end=width(dna)))
seqlengths(gr)<-width(dna)
seqlevels(gr)[8]<-"Mt"
## import bedgraph (treat)
bg <-import.bedGraph("Myresearch/ZZ/2020.2.7/200116/5minimap2.bedgraph")
bg<-bg[seqnames(bg)%in%c("2L","2R","3L","3R","4","X","Y","mitochondrion_genome")]
bg$score <-log(bg$score+1,2)
seqlevels(bg)[8]<-"Mt"
seqlevels(bg)<-seqlevels(gr)
seqlengths(bg)<-seqlengths(gr)

## import bedgraph (control)
bg2 <-import.bedGraph("Myresearch/ZZ/2020.2.7/200117/5minimap2.bedgraph")
bg2<-bg2[seqnames(bg2)%in%c("2L","2R","3L","3R","4","X","Y","mitochondrion_genome")]
bg2$score <-log(bg2$score+1,2)
seqlevels(bg2)[8]<-"Mt"
seqlevels(bg2)<-seqlevels(gr)
seqlengths(bg2)<-seqlengths(gr)


## make track name
gr1<-gr[1]
gr1$name <-"200116"

gr2<-gr[1]
gr2$name <-"200117"

#plot
p<-ggplot()+ layout_circle(gr, geom = "ideo", fill = aes(fill=seqnames), radius = 7, trackWidth = 2)+
  layout_circle(gr, geom = "scale",radius = 9, size = 4,trackWidth = 0.5)+
  layout_circle(gr, geom = "text", aes(label = seqnames), vjust = 0, size = 8,radius = 5)+
  layout_circle(bg, geom = "point", color = "red",alpha=0.8, radius = 5,trackWidth = 2, size=1,grid = TRUE, aes(y = score),ylim=c(0,20))+
  layout_circle(gr1, geom = "text", aes(label = name), vjust = 0, size = 5,radius = 3)+
  layout_circle(bg2, geom = "point", color = "red",alpha=0.8, radius = 3,trackWidth = 2, size=1,grid = TRUE, aes(y = score),ylim=c(0,20))+
  layout_circle(gr2, geom = "text", aes(label = name), vjust = 0, size = 5,radius = 1)+
  theme(legend.position = "none")

