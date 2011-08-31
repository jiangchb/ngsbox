pdf(file="cov_hist.pdf", width=8, height=8)

all<-read.table("cov_hist.all.out")
maxx<-16
maxy<-max(all$V2) + 8000000

lambda=sum(all$V1*all$V2) / sum(all$V2) 
xpois<-seq(0,50)
ypois<-c()
for (i in xpois) {
	ypois[i+1]<-(lambda^i * exp(-lambda) / factorial(i)) * sum(all$V2)
}

plot(all$V1[1:maxx], all$V2[1:maxx], type="n", xlim=c(0, maxx), ylim=c(0, maxy), axes=T, xlab="Coverage", ylab="Frequency", lwd=3)

polygon(c(0, xpois), c(0, ypois), col="darkgrey", border=NA )

lines(all$V1[1:maxx], all$V2[1:maxx], col="goldenrod1", lwd=3)

legend("topright", c(paste("Poisson distribution (lambda=", lambda, ")", sep=""), "Coverage distribution"), col=c("darkgrey", "goldenrod1"), bty="n", lwd=3)

dev.off()
