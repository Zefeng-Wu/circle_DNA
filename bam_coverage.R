library(ggbio)
library(GenomicAlignments)

ga_list<-c()
ga_name<-c()
file_list<-dir(path = "/home/wuzefeng/opt/zz",full.names = TRUE,pattern = ".bam")
for (file in file_list){
  if (!endsWith(file,".bai")&!endsWith(file,".input.bam")){
    message(file)
    temp_ga<-readGAlignments(file,use.names=TRUE)
    temp_plot<-autoplot(temp_ga,geom="line",stat="coverage")+ylab(paste(stringr::str_split_fixed(basename(file),pattern = "\\.",n = 3)[,c(1,2)],collapse = "_"))
    ga_list<-c(ga_list,temp_plot)
    ga_name<-c(ga_name,basename(file))
     }
}

tracks(ga_list)
tracks(Reads=p1, Coverage=p2,Transcripts=p4, heights = c(0.3, 0.2, 0.35),track.plot.color=c("red","blue","orange"),track.bg.color="lightyellow",label.bg.fill="red",theme=theme_minimal()) + ylab("")

