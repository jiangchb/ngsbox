sens<-read.table("sens_vs_depth.txt")
spec<-read.table("spec_vs_depth.txt")

postscript(file="sens_spec_subsampling.ps", width=16, height=8)

layoutmatrix<-matrix(data=c(1,2), nrow=1, ncol=2)
layout(layoutmatrix)

plot(sens$V1, sens$V2, type="n", col="red", lwd=3, axes=F, xlab="Sequencing depth", ylab="Sensitivity", xlim=c(2,25), ylim=c(0.3, 1), main="Sensitivity")

for (i in 1:10) {
        abline(h=i/10, col="grey", lty=2, lwd=2)
}

lines(sens$V1, sens$V2, type="b", col="red", lwd=3)
lines(sens$V1, sens$V3, type="b", col="limegreen", lwd=3)
lines(sens$V1, sens$V4, type="b", col="orange", lwd=3)

axis(1, labels=c("3", "6", "12", "18", "24"), at=c("3", "6", "12", "18", "24"))
axis(2, labels=c("30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%"), at=c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), las=1)

legend(18, 0.5, c("Contigs", "SNPs", "1 bp indels"), col=c("orange", "red", "limegreen", "orange"), bty="n", lwd=3)



plot(spec$V1, spec$V2, type="n", col="red", lwd=3, axes=F, xlab="Sequencing depth", ylab="Specitivity", xlim=c(2,25), ylim=c(0.95, 1), main="Specitivity")

for (i in 1:100) {
        abline(h=i/100, col="grey", lty=2, lwd=2)
}

lines(spec$V1, spec$V2, type="b", col="red", lwd=3)
lines(spec$V1, spec$V3, type="b", col="limegreen", lwd=3)
lines(spec$V1, spec$V4, type="b", col="orange", lwd=3)

axis(1, labels=c("3", "6", "12", "18", "24"), at=c("3", "6", "12", "18", "24"))
axis(2, labels=c("95%", "96%", "97%", "98%", "99%", "100%"), at=c(0.95, 0.96, 0.97, 0.98, 0.99, 1.0), las=1)

legend(18, 0.5, c("Contigs", "SNPs", "1 bp indels"), col=c("orange", "red", "limegreen", "orange"), bty="n", lwd=3)

legend(18, 0.96, c("Contigs", "SNPs", "1 bp indels"), col=c("orange", "red", "limegreen", "orange"), bty="n", lwd=3)

dev.off()





