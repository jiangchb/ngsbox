pebur<-read.table("sv_indel.pe.bur-0.txt")
mpbur<-read.table("sv_indel.mp.bur-0.txt")
pecol<-read.table("sv_indel.pe.col-0.txt")

pdf(file="sv_length_dist.pdf", width=12, height=8)

layoutmat=matrix(data=c(1,3,2,4), nrow=2, ncol=2)
layout(layoutmat)

col="darkgrey"

hist(pebur$V1, breaks=35000, xlim=c(-500,1500), col=col, main="Paired-end analysis of Bur-0", xlab="Length (bp) of insertions (-) and deletions (+)")
hist(mpbur$V1, breaks=3000, xlim=c(-5000,15000), col=col, main="Mate pair analysis of Bur-0", xlab="Length (bp) of insertions (-) and deletions (+)")


hist(pecol$V1, breaks=35000, xlim=c(-500,1500), col=col, main="Paired-end analysis of reference re-sequencing", xlab="Length (bp) of insertions (-) and deletions (+)")



dev.off()

