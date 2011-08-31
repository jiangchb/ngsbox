pdf(file="gc_by_cov_corr.pdf", width=12, height=12)

gcs<-read.table("gcs_by_cov.out")
cov<-as.double(gcs$V3)
gc41<-as.double(gcs$V4)
gc71<-as.double(gcs$V5)
gc101<-as.double(gcs$V6)
gc151<-as.double(gcs$V7)
gc201<-as.double(gcs$V8)
gc501<-as.double(gcs$V9)
gc1001<-as.double(gcs$V10)

corr41pear<-cor.test(cov, gc41, method="pearson")
corr71pear<-cor.test(cov, gc71, method="pearson")
corr101pear<-cor.test(cov, gc101, method="pearson")
corr151pear<-cor.test(cov, gc151, method="pearson")
corr201pear<-cor.test(cov, gc201, method="pearson")
corr501pear<-cor.test(cov, gc501, method="pearson")
corr1001pear<-cor.test(cov, gc1001, method="pearson")


corr41spearman<-cor.test(cov, gc41, method="spearman")
corr71spearman<-cor.test(cov, gc71, method="spearman")
corr101spearman<-cor.test(cov, gc101, method="spearman")
corr151spearman<-cor.test(cov, gc151, method="spearman")
corr201spearman<-cor.test(cov, gc201, method="spearman")
corr501spearman<-cor.test(cov, gc501, method="spearman")
corr1001spearman<-cor.test(cov, gc1001, method="spearman")

corr41ken<-cor.test(cov, gc41, method="kendall")
corr71ken<-cor.test(cov, gc71, method="kendall")
corr101ken<-cor.test(cov, gc101, method="kendall")
corr151ken<-cor.test(cov, gc151, method="kendall")
corr201ken<-cor.test(cov, gc201, method="kendall")
corr501ken<-cor.test(cov, gc501, method="kendall")
corr1001ken<-cor.test(cov, gc1001, method="kendall")






layoutmaxtrix<-matrix(data=c(1,2,3,4,5,6,7,8,9,10,11,12), nrow=4, ncol=3, byrow=T)
layout(layoutmaxtrix)

plot(cov, gc41, xlim=c(0, 100), ylim=c(0,100), axes=F, col="grey65", cex=0.3, pch=20, xlab="coverage", ylab="GC content (ws=41)")
axis(1)
axis(2)
plot(cov, gc71, xlim=c(0, 100), ylim=c(0,100), axes=F, col="grey65", cex=0.3, pch=20, xlab="coverage", ylab="GC content (ws=71)")
axis(1)
axis(2)
plot(cov, gc101, xlim=c(0, 100), ylim=c(0,100), axes=F, col="grey65", cex=0.3, pch=20, xlab="coverage", ylab="GC content (ws=101)")
axis(1)
axis(2)
plot(cov, gc151, xlim=c(0, 100), ylim=c(0,100), axes=F, col="grey65", cex=0.3, pch=20, xlab="coverage", ylab="GC content (ws=151)")
axis(1)
axis(2)
plot(cov, gc201, xlim=c(0, 100), ylim=c(0,100), axes=F, col="grey65", cex=0.3, pch=20, xlab="coverage", ylab="GC content (ws=201)")
axis(1)
axis(2)
plot(cov, gc501, xlim=c(0, 100), ylim=c(0,100), axes=F, col="grey65", cex=0.3, pch=20, xlab="coverage", ylab="GC content (ws=501)")
axis(1)
axis(2)
plot(cov, gc1001, xlim=c(0, 100), ylim=c(0,100), axes=F, col="grey65", cex=0.3, pch=20, xlab="coverage", ylab="GC content (ws=1001)")
axis(1)
axis(2)


x<-c(1,2,3,4,5,6,7)
speary<-c(corr41spearman$estimate, corr71spearman$estimate, corr101spearman$estimate, corr151spearman$estimate, corr201spearman$estimate, corr501spearman$estimate, corr1001spearman$estimate)
peary<-c(corr41pear$estimate, corr71pear$estimate, corr101pear$estimate, corr151pear$estimate, corr201pear$estimate, corr501pear$estimate, corr1001pear$estimate)
keny<-c(corr41ken$estimate, corr71ken$estimate, corr101ken$estimate, corr151ken$estimate, corr201ken$estimate, corr501ken$estimate, corr1001ken$estimate)

plot(c(1), c(1), type="n", axes=F, ylab="", xlab="")
plot(c(1), c(1), type="n", axes=F, ylab="", xlab="")

plot(x, speary, col="firebrick", lwd=1, pch=2, axes=F, xlab="Window size", ylab="Spearman's rank correlation rho")
axis(1, at=c(1,2,3,4,5,6,7), labels=c("41", "71", "101", "151", "201", "501", "1001"))
axis(2, las=1)

plot(x, peary, col="orange3", lwd=1, pch=2, axes=F, main="GC window size correlations", xlab="Window size", ylab="Pearson's product-moment correlation")
axis(1, at=c(1,2,3,4,5,6,7), labels=c("41", "71", "101", "151", "201", "501", "1001"))
axis(2, las=1)

plot(x, keny, col="black", lwd=1, pch=2, axes=F, xlab="Window size", ylab="Kendall's rank correlation tau")
axis(1, at=c(1,2,3,4,5,6,7), labels=c("41", "71", "101", "151", "201", "501", "1001"))
axis(2, las=1)



dev.off()

