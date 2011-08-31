
args <- commandArgs()
chr<-as.numeric(args[5])
beg<-as.numeric(args[6])
end<-as.numeric(args[7])

cov<-read.table("cov.txt")

pdf(file=paste("cov_", chr, "_", beg, "_", end, ".pdf", sep=""), width=15, height=6)

ylim<-c(0, max(cov$V2))
xlim<-c(min(cov$V1), max(cov$V1))
ratio<-cov$V2/cov$V3
ratio_col<-cov$V7/cov$V8

### Plot observed coverage and expected coverage
plot(cov$V1, cov$V2, type="l", main="Coverage vs. Expected Coverage in Duplicated Regions", xlim=xlim, ylim=ylim, col="blue", axes=F, xlab="Position", ylab="Coverage")
lines(cov$V1, cov$V3, type="l", xlim=xlim, ylim=ylim, col="firebrick")
lines(cov$V1, cov$V6, type="h", xlim=xlim, ylim=ylim, col="green", lwd=2)
#lines(cov$V1, cov$V7, type="l", xlim=xlim, ylim=ylim, col="grey")
#lines(cov$V1, cov$V8, type="l", xlim=xlim, ylim=ylim, col="green")
axis(1)
axis(2, at=c(0, max(cov$V2)), labels=c("0", max(cov$V2)))
legend(min(cov$V1) + 10, max(cov$V2), pch=1, col=c("blue", "firebrick", "green"), legend=c("Obs Cov", "Exp Cov", "CVP"))


### Plot coverage ratio
#ylim<-c(0,3)
#plot(cov$V1, ratio, type="l", main="Coverage Ration in Duplicated Regions", xlim=xlim, ylim=ylim, col="blue", axes=F, xlab="Position", ylab="Obs Cov / Exp Cov")
#lines(cov$V1, ratio_col, type="l", xlim=xlim, ylim=ylim, col="firebrick")
#lines(cov$V1, cov$V6, type="h", xlim=xlim, ylim=ylim, col="green", lwd=2)
#axis(1)
#axis(2, at=c(0, 3), labels=c("0", "3"))
#legend(min(cov$V1) + 10, 3, pch=1, col=c("blue"), legend=c("Obs/Exp"))


### Plot average hits and mismatches
#ylim<-c(0,5)

#plot(cov$V1, cov$V4, type="l", main="Average Hits and Mismatches in Duplicated Regions", xlim=xlim, ylim=ylim, col="blue", axes=F, xlab="Position", ylab="Count")
#lines(cov$V1, cov$V5, type="l", xlim=xlim, ylim=ylim, col="firebrick")
#lines(cov$V1, cov$V6, type="h", xlim=xlim, ylim=ylim, col="green", lwd=2)
#axis(1)
#axis(2, at=c(0, 5), labels=c("0", "5"))
#legend(min(cov$V1) + 10, 5, pch=1, col=c("blue", "firebrick", "green"), legend=c("Avg Hits", "Avg MM", "CVP"))


#axis(4, at=c(0, max(cov$V3)), labels=c("0", "1"))

dev.off()
