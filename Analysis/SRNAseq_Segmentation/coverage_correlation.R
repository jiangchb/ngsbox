
args <- commandArgs()
filename<-(args[5])
sample1name<-(args[6])
sample2name<-(args[7])

cov<-read.table(paste(filename, ".txt", sep=""))

pdf(file=paste(filename, ".pdf", sep=""), height=12, width=12)

xlim<-c(0, max(log2(cov$V1), log2(cov$V2)))
ylim<-c(0, max(log2(cov$V1), log2(cov$V2)))

plot(log2(cov$V1), log2(cov$V2), type="p", main=paste("sRNA Count Sample ", sample1name, " vs. Sample ", sample2name) , xlim=xlim, ylim=ylim, col="blue", xlab=sample1name, ylab=sample2name, pch=20)

correlation<-cor.test(cov$V1, cov$V2, method="pearson")
text( x=1, y=max(log2(cov$V1), log2(cov$V2)) - 1, labels=paste("r square: ", correlation$estimate) )
text( x=1, y=max(log2(cov$V1), log2(cov$V2)) - 2, labels=paste("p value: ", correlation$p.value))

dev.off()
