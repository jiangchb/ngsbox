pdf(file="cov_hist.pdf", width=8, height=8)

all<-read.table("cov_hist.all.out")
low<-read.table("cov_hist.low.out")
up<-read.table("cov_hist.up.out")
mid<-read.table("cov_hist.mid.out")
random<-read.table("cov_artificial_sampling.hist")

#randompoisson<-rpois(sum(all$V2), all$V1[all$V2==max(all$V2)])
#poishist<-hist(randompoisson, breaks=c(0:maxx), plot=F)

lambda=15.0130 # select avg(coverage) from consensus_col where seg_flag = "u" and ref_base != "N";
xpois<-seq(0,50)
ypois<-c()
for (i in xpois) {
	ypois[i+1]<-(lambda^i * exp(-lambda) / factorial(i)) * sum(all$V2)
}


maxy<-12000000
maxx<-51
plot(all$V1[1:maxx], all$V2[1:maxx], type="n", xlim=c(0, maxx), ylim=c(0, maxy), axes=F, xlab="Coverage", ylab="Frequency [10^6]", lwd=3)
axis(2, labels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), at=c(0, 1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 8000000, 9000000, 10000000, 11000000, 12000000), las=1)
axis(1, labels=c(0,10,20,30,40,50,60,70,80), at=c(0,10,20,30,40,50,60,70,80))


polygon(c(0, xpois), c(0, ypois), col="lightgrey", border=NA )
polygon(c(0, random$V2), c(0, random$V1), col="darkgrey", border=NA)
lines(xpois, ypois, col="lightgrey", lty=2, lwd=2)


lines(all$V1[1:maxx], all$V2[1:maxx], col="black", lwd=3)
lines(up$V1[1:maxx], up$V2[1:maxx], lwd=3, col="orange3")
lines(low$V1[1:maxx], low$V2[1:maxx], lwd=3, col="firebrick")
lines(mid$V1[1:maxx], mid$V2[1:maxx], lwd=3, col="navajowhite4")


legend("topright", c(paste("Poisson distribution (lambda=", lambda, ")", sep=""), "Random sampling of 36 mers", "Coverage distribution", "Coverage distribution (GC content <= 25)", "Coverage distribution (GC content >= 45)", "Coverage distribution (GC content between 32 and 40)"), col=c("lightgrey", "darkgrey", "black", "firebrick", "orange3", "navajowhite4"), bty="n", lwd=3)




dev.off()

