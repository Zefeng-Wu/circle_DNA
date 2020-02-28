library(genoPlotR)
library(ggbio)
data<-read_comparison_from_blast("~/Desktop/constoTEs_beagle.blout")
data$s2<-ifelse(data$start2<data$end2,data$start2,data$end2)
data$e2<-ifelse(data$start2<data$end2,data$end2,data$start2)
data<-as.data.frame(data)
data$direction<-ifelse(data$direction==1,"+","-")
library(GenomicRanges)
gr<-makeGRangesFromDataFrame(df = data[data$name2=="HMS-Beagle",],seqnames.field = "name2",
                             strand.field = "direction",
                             start.field = "s2",
                             end.field = "e2",
                             keep.extra.columns = TRUE)
gr$direct <-as.character(strand(gr))
autoplot(gr,aes(fill=direct))

reduced_gr<-unlist(reduce(split(gr, elementMetadata(gr)$name1)))
reduced_gr$id<-names(reduced_gr)
reduced_gr$direction<-as.character(strand(reduced_gr))
autoplot(reduced_gr,aes(fill=direction)) # remove second alignment for each read in a similar position
autoplot(reduced_gr[1:100], layout = "circle",aes(col=direction))

## grangelist visuliation
xx<-reduce(split(gr, elementMetadata(gr)$name1))
ggplot()+geom_alignment(xx[1:100],type="id")
autoplot(grl, aes(type = id))
#008da576-7b24-4000-8f90-ae702802d2d4_cons0_3343_95_2947_958_3.1_0 # an interesting case


## another visuliaztion
autoplot(temp_ga, colour=as.factor(strand(temp_ga)))
