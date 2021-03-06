pdf(file="SV_sim_insertion_comp.pdf", width=10, height=10)

#szr="49-72"
szr="73-100"

red="#F90000"
darkred="#CD2626"
green="#86CC00"
darkgreen="#6E8B3D"
yellow="#FFC125"
darkyellow="#FF8B00"


#data<-read.table(paste("insertions.", szr,".txt", sep=""))
data<-read.table("tmp")

layoutmat<-matrix(data=c(1,3,2,4), nrow=2, ncol=2)
layout(layoutmat)
breaks=20
max=80
max2=75

#hist(data$V5-data$V2, breaks=breaks, xlim=c(-30,30), col="darkcyan", xlab="0 = left of insertion ", main="Distance to insertion start", ylim=c(0,max), axes=F)
hist(data$V5-data$V2, breaks=5000, xlim=c(-20,20), col="darkcyan", xlab="0 = left of insertion ", main="Distance to insertion start", axes=F, ylim=c(0,max))
#axis(1, labels=c(-20, -10, 0, 10, 20), at=c(-19.5, -9.5, 0.5, 10.5, 20.5))
axis(1)
axis(2)

#hist(((data$V2+1)-data$V6)*(-1), breaks=breaks, xlim=c(-30,30), col="darkcyan", xlab="0 = right of insertion", main="Distance to insertion end", ylim=c(0,max), axes=F)
hist(((data$V2+1)-data$V6)*(-1), breaks=5000, xlim=c(-20,20), col="darkcyan", xlab="0 = right of insertion", main="Distance to insertion end", axes=F, ylim=c(0,max))   
#axis(1, labels=c(-20, -10, 0, 10, 20), at=c(-19.5, -9.5, 0.5, 10.5, 20.5))
axis(1)
axis(2)

hist(data$V7-data$V4, breaks=19, xlim=c(-30,30), col="darkkhaki", xlab="0 = Matches insertion length", main="Distance to predicted length", ylim=c(0,max2))



conc<-read.table(paste("insertions.insertsizes.", szr,".txt", sep=""))
disc<-read.table(paste("insertions.insertsizes_discordant.", szr,".txt", sep=""))

layoutmat<-matrix(data=c(1,2), nrow=2, ncol=1)
layout(layoutmat)

conc_hist<-hist(conc$V1, breaks=200, plot=F)
abline(v=mean(conc$V1))
disc_hist<-hist(disc$V1, breaks=100, plot=F)
abline(v=mean(disc$V1))

plot(conc_hist$mids, conc_hist$density, type="l", xlim=c(50, 250), ylim=c(0, 0.06), axes=F, xlab="Insert size", ylab="Frequency", main="Insert distributions of deletions vs. background", col=darkyellow, lwd=2)
lines(disc_hist$mids, disc_hist$density, col=darkred, lwd=2)

axis(1)
axis(2)
abline(v=mean(disc), col=darkred)
abline(v=mean(conc), col=darkyellow)


dev.off()

