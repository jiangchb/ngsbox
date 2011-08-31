args <- commandArgs()
mat<-args[5]

pdf(file=paste(mat, ".pdf", sep="") , width=18, height=16)

layoutmat<-matrix(data=c(1,2), nrow=2, ncol=1)
layout(layoutmat)

########## Read in data ##############

data<-read.table(mat)

gc<-data$V1
cov<-data$V4
std<-data$V6
supp<-data$V2
avgcov<-sum(cov*supp)/(sum(supp))
maxgc<-max(cov+std)

allcov<-data$V9
allstd<-data$V11
allsupp<-data$V7
allavgcov<-sum(allcov*allsupp)/(sum(allsupp))
allmaxgc<-max(allcov+allstd)

########## Plot nonrep, non-oversampled ##########

plot(gc, cov, ylab="avg(coverage)", xlab="Local GC content [%]", ylim=c(0, max(cov+std)), col="orange3", type="h", lwd=7, axes=F, main="Coverage and GC content at uniquely mapped sites")

axis(2)
axis(1, at=c(0, 10,20,30,40,50,60,70,80,90,100), label=c("0", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"))
axis(4, at=c(0, max(cov+std)/2, max(cov+std)), labels=c("0", max(supp)/2, max(supp)))

lines(c(0, 100), c(avgcov+1, avgcov+1), lty=2)


lines(seq(0,100), ((supp/max(supp))*max(cov+std)),lty=2, col="firebrick", lwd=2)

for (i in 1:101) {
	if (cov[i] > 0) {
	        segments(gc[i],cov[i],gc[i],cov[i]+std[i], col="grey", lwd = 2)
	}
}

text(gc, cov, adj=c(-.5,.7), labels=supp, cex=0.65, srt=90, col="black")
text(0, avgcov+1+(max(cov+std)/40), labels=paste("Obs. avg. cov=", round(avgcov+1, digits=1), sep=""), pos=4)

legend("topleft", legend=c("Observed coverage and its", "  Stddev",  "  Support distribution", "  Average coverage"), bty="n", fill=c("orange3", "grey", "firebrick", "black"))

########## Plot all ###########

plot(gc, allcov, ylab="avg(coverage)", xlab="Local GC content [%]", ylim=c(0, max(allcov+allstd)), col="orange3", type="h", lwd=7, axes=F, main="Coverage and GC content at all sites")

axis(2)
axis(1, at=c(0, 10,20,30,40,50,60,70,80,90,100), label=c("0", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"))
axis(4, at=c(0, max(allcov+allstd)/2, max(allcov+allstd)), labels=c("0", max(allsupp)/2, max(allsupp)))

lines(c(0, 100), c(allavgcov+1, allavgcov+1), lty=2)

lines(seq(0,100), ((allsupp/max(allsupp))*max(allcov+allstd)),lty=2, col="firebrick", lwd=2)

for (i in 1:101) {
	if (allcov[i] > 0) {
	        segments(gc[i],allcov[i],gc[i],allcov[i]+allstd[i], col="grey", lwd = 2)
	}
}

text(gc, allcov, adj=c(-.5,.7), labels=allsupp, cex=0.65, srt=90, col="black")
text(0, allavgcov+1+(max(allcov+allstd)/40), labels=paste("Obs. avg. cov=", round(allavgcov+1, digits=1), sep=""), pos=4)

legend("topleft", legend=c("Observed coverage and its", "  Stddev",  "  Support distribution", "  Average coverage"), bty="n", fill=c("orange3", "grey", "firebrick", "black"))


#################################

dev.off()
