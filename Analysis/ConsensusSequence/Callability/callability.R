args <- commandArgs()
chrsize<-read.table(args[5])

pdf(file=args[6], width=10, height=12)
layoutmat=matrix(data=c(1,2,3,4,5), ncol=1, nrow=5)
layout(layoutmat)

options(scipen=999999999)
ls=5000000

data<-read.table(args[7])
max=max(data$V4)
winstep=args[8]
winsize=args[9]


for (chr in 1:(length(chrsize$V1))) {

	chrname=chrsize$V1[chr]

	plot(data$V2[data$V1[]==chrname], data$V4[data$V1[]==chrname], ylim=c(0, max), xlim=c(0, max(chrsize$V2)), type="l", axes=F, xlab=paste("Chromosome ", chrname, sep=""), ylab="", main=paste("winstep:", winstep, " winsize:", winsize, sep=""))

	labels=c(1, seq(ls, chrsize$V2[chrsize$V1[]==chrname], by=ls), chrsize$V2[chrsize$V1[]==chrname])
	axis(1, label=labels, at=labels)
	axis(2, las=1, labels=c("0", "1"), at=c(0, max))
}

dev.off()

