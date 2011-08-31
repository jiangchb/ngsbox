args <- commandArgs()
eco=args[5]
input=args[6]
pdffile=args[7]

pdf(file=pdffile, width=16, height=10)

data<-read.table(input)
blast<-read.table("/ebio/abt6_analysis/nobackup/data/Shore/Plants/ATH/DM2/significant_blasthits.out")

layoutmat<-matrix(data=c(1,2,3,4), nrow=4, ncol=1)
layout(layoutmat)

ymax=25

plot(data$V2[data$V1=="Bla-1"], data$V3[data$V1=="Bla-1"], type="l", col="red", xlab="Bla-1", ylab="Coverage", axes=F, main=eco, ylim=c(-5, min(max(max(data$V3[data$V1=="Bla-1"]), max(data$V4[data$V1=="Bla-1"])), ymax))) 
lines(data$V2[data$V1=="Bla-1"], data$V4[data$V1=="Bla-1"], type="l", col="blue")
axis(1)
axis(2)
segments(blast$V2[blast$V1=="Bla-1"], -4, blast$V3[blast$V1=="Bla-1"], -4)

plot(data$V2[data$V1=="Ler-1"], data$V3[data$V1=="Ler-1"], type="l", col="red", xlab="Ler-1", ylab="Coverage", axes=F, ylim=c(-5, min(max(max(data$V3[data$V1=="Ler-1"]), max(data$V4[data$V1=="Ler-1"])), ymax)))
lines(data$V2[data$V1=="Ler-1"], data$V4[data$V1=="Ler-1"], type="l", col="blue")
axis(1)
axis(2)
segments(blast$V2[blast$V1=="Ler-1"], -4, blast$V3[blast$V1=="Ler-1"], -4)

plot(data$V2[data$V1=="Uk1"], data$V3[data$V1=="Uk1"], type="l", col="red", xlab="Uk1", ylab="Coverage", axes=F, ylim=c(-5, min(max(max(data$V3[data$V1=="Uk1"]), max(data$V4[data$V1=="Uk1"])), ymax)))
lines(data$V2[data$V1=="Uk1"], data$V4[data$V1=="Uk1"], type="l", col="blue")
axis(1)
axis(2)
segments(blast$V2[blast$V1=="Uk1"], -4, blast$V3[blast$V1=="Uk1"], -4)

plot(data$V2[data$V1=="Col-0"], data$V3[data$V1=="Col-0"], type="l", col="red", xlab="Col-0", ylab="Coverage", axes=F, ylim=c(-5, min(max(max(data$V3[data$V1=="Col-0"]), max(data$V4[data$V1=="Col-0"])), ymax)))
lines(data$V2[data$V1=="Col-0"], data$V4[data$V1=="Col-0"], type="l", col="blue")
axis(1)
axis(2)
segments(blast$V2[blast$V1=="Col-0"], -4, blast$V3[blast$V1=="Col-0"], -4)

dev.off()

