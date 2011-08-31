pseudo<-read.table("quality.dist.pseudo_snps")
snps<-read.table("quality.hist.snps")
all<-read.table("scaffold_1.qual.sort")

pdf(file="lyrata.snps.by.quality.pdf", width=16, height=16)

layoutmat<-matrix(data=c(1,2,3,4,5,6), nrow=3, ncol=2, byrow=T)
layout(layoutmat)


hist(snps$V3, breaks=100, xlab="RefQuality at predicted European fixation",col="blue", main="Reference Qualities at changes in the European Arabidosis lyrata genome pool (scaffold 1 only)")
hist(snps$V3, breaks=100, xlab="RefQuality at predicted European fixation", ylim=c(0,100), col="blue", main="")

hist(pseudo$V3, breaks=100, xlab="RefQuality at majority vote against ref base", col="red", main="")
hist(pseudo$V3, breaks=100, xlab="RefQuality at majority vote against ref base", ylim=c(0,100), col="red", main="")

plot(all$V2, all$V1, type="h", ylab="Frequency", xlab="All quality scores on scaffold 1", col="limegreen", axes=F, lwd=2, xlim=c(0,60))
axis(1)
axis(2)

plot(all$V2, all$V1, type="h", ylab="Frequency", xlab="All quality scores on scaffold 1", col="limegreen", axes=F, ylim=c(0,10000), lwd=2, xlim=c(0,60))
axis(1)
axis(2)





dev.off()
