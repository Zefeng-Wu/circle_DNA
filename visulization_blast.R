library(genoPlotR)
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
