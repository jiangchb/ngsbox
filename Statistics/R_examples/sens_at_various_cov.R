sens<-read.table("sens_at_various_cov.txt")
spec<-read.table("spec_at_various_cov.txt")


postscript(file="sens_and_spec_at_various_min_covs.ps", width=16, height=8)

layoutmatrix<-matrix(data=c(1, 2), nrow=1, ncol=2)
layout(layoutmatrix)

plot(sens$V1, sens$V2, xlim=c(0,7), ylim=c(0,0.90), col="red", type="l", axes=F, xlab="Minimal coverage", ylab="Sensitivity", main="Sensitivity")
lines(sens$V1, sens$V3, col="limegreen")
lines(sens$V1, sens$V4, col="orange")
lines(sens$V1, sens$V5, col="yellow")
#lines(sens$V1, sens$V6, col="white")
axis(1)
axis(2, labels=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90), at=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), las=2)


legend("topright", c("24 fold", "18 fold", "12 fold", "6 fold"), col=c("red", "limegreen", "orange", "yellow"), bty="n", lwd=3)

abline(h=0.9, lty=2, col="grey")
abline(h=0.8, lty=2, col="grey")
abline(h=0.7, lty=2, col="grey")
abline(h=0.6, lty=2, col="grey")
abline(h=0.5, lty=2, col="grey")
abline(h=0.4, lty=2, col="grey")
abline(h=0.3, lty=2, col="grey")
abline(h=0.2, lty=2, col="grey")
abline(h=0.1, lty=2, col="grey")
abline(h=0, lty=2, col="grey")

plot(spec$V1, spec$V2, xlim=c(0,7), ylim=c(0.98,1.0), col="red", type="l", axes=F, xlab="Minimal coverage", ylab="Specificity", main="Specificity")
lines(spec$V1, spec$V3, col="limegreen")
lines(spec$V1, spec$V4, col="orange")
lines(spec$V1, spec$V5, col="yellow")
#lines(spec$V1, spec$V6, col="white")
axis(1)
axis(2, labels=c(98, 99, 100), at=c(0.98, 0.99, 1.0), las=2)

legend("topright", c("24 fold", "18 fold", "12 fold", "6 fold"), col=c("red", "limegreen", "orange", "yellow"), bty="n", lwd=3)

abline(h=1.0, lty=2, col="grey")
abline(h=0.99, lty=2, col="grey")
abline(h=0.98, lty=2, col="grey")



dev.off()





