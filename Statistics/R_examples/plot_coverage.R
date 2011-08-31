cov<-read.table("cov.txt")

pdf(file="coverage.pdf", height=8, width=15)

ylim<-c(0, max(cov$V2))
xlim<-c(min(cov$V1), max(cov$V1))

### Plot observed coverage and expected coverage
plot(cov$V1, cov$V2, type="l", main="Coverage vs. Expected Coverage in Duplicated Regions", xlim=xlim, ylim=ylim, col="blue", axes=F, xlab="Position", ylab="Coverage")
lines(cov$V1, cov$V4, type="l", xlim=xlim, ylim=ylim, col="firebrick")
axis(1)
axis(2, at=c(0, max(cov$V2)), labels=c("0", max(cov$V2)))


###Plot observed nonrepetitive coverage and expected coverage
plot(cov$V1, cov$V3, type="l", xlim=xlim, ylim=ylim, col="blue", axes=F, xlab="Position", ylab="Nonrepetitive  Coverage")
lines(cov$V1, cov$V4, type="l", xlim=xlim, ylim=ylim, col="firebrick")
axis(1)
axis(2, at=c(0, max(cov$V2)), labels=c("0", max(cov$V2)))


#axis(4, at=c(0, max(cov$V3)), labels=c("0", "1"))

dev.off()
