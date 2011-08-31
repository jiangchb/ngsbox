col<-read.table("col.cov_vs_exp_cov.txt")

pdf(file="exp_vs_obs_col.pdf", width=8, height=8)

layoutmat<-matrix(data=c(1, 1), nrow=1, ncol=2, byrow=T)
layout(layoutmat)

##################################################################################
####### Col-0 with exp_cov

axeslim=1000

coldepth=16.0

x<-1:max(col$V2)

allcor<-cor.test(col$V2, col$V1, method="spearman")

xlabel<-paste("Expected coverage ex_p\nSpearman's rho: ", round(allcor$estimate, digits=6), " (p-value=", allcor$p.value,")", sep="")

plot(col$V2, col$V1, xlab=xlabel, ylab="Observed coverage", xlim=c(0,axeslim), ylim=c(0,axeslim), type="n", axes=F, main="Col-0 - Observed vs. expected coverage")
points(col$V2, col$V1, cex=0.2, pch=19, col="grey50")


axis(1)
axis(2, las=2)

#####################################################

#axeslim=1000

#coldepth=16.0

#x<-1:max(col$V2)

#allcor<-cor.test(col$V3, col$V1, method="spearman")

#xlabel<-paste("Expected coverage\nCorrelations: ", round(allcor$estimate, digits=2), " (p-value=", allcor$p.value,")", sep="")

#plot(col$V3, col$V1, xlab=xlabel, ylab="Observed coverage", xlim=c(0,axeslim), ylim=c(0,axeslim), type="n", axes=F, main="Col-0 - Observed vs. expected coverage")
#points(col$V3, col$V1, cex=0.2, pch=19, col="grey50")


#axis(1)
#axis(2)



dev.off()

