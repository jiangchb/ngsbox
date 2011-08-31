args <- commandArgs()

pdf(file="plot_errors.pdf", width=20, height=10)

########## Read in files ############

read_length_dist<-read.table(args[5])
mismatch_dist<-read.table(args[6])
hits_dist<-read.table(args[7])
prb<-read.table(args[8])
qCal<-read.table(args[9])
chas<-read.table(args[10])
position<-read.table(args[11])

read_length_dist_2<-read.table(args[12])
mismatch_dist_2<-read.table(args[13])
hits_dist_2<-read.table(args[14])
prb_2<-read.table(args[15])
qCal_2<-read.table(args[16])
chas_2<-read.table(args[17])
position_2<-read.table(args[18])



################ Layout for the first slides ##################################################################################################

layoutmat<-matrix(data=c(1,2,3,4), nrow=2, ncol=2, byrow=F)
layout(layoutmat)



############# Error by PRB -- read 1 ###############

supp<-(prb$V3+prb$V2)
rel<-prb$V2/supp
rel[is.nan(rel)]<-0

plot(prb$V1, rel, type="h", lwd=3, col="orange3", axes=F, ylab="Percent errors of observed positions [%]", xlab="prb (Qs)", main="Errors in mapped reads at uniquely mapped positions, excluding ambiguous calls", ylim=c(0, max(rel)+0.05))
axis(1, at=seq(-5, 40, by=5))
axis(2, at=c(0, max(rel)/2, max(rel)), labels=c(0, round(max(rel)/2*100, 2), round(max(rel)*100, 2)))
axis(4, at=c(0, max(rel)/2, max(rel)), labels=c("0", max(supp[1:45])/2, max(supp[1:45])))

text(prb$V1, rel, adj=c(-.5,.7), labels=prb$V2, cex=0.65, srt=90, col="black")

lines(seq(-5,40), supp/max(supp[1:45])*max(rel),lty=2, col="firebrick", lwd=2)

legend("topright", legend=c("Percent of observed errors", "Number of observed errors", "Support distribution"), bty="n", fill=c("orange3", "black", "firebrick"))


######### Error by PRB (cumulative) --read 1 ########

cumerror<-c()
sumerror<-0
for (i in 1:length(rel)) {
	sumerror<-sumerror+prb$V2[i]
	cumerror[i]<-sumerror
}
cumerror<-cumerror/sumerror

cumcorr<-c()
sumcorr<-0
for (i in 1:length(rel)) {
        sumcorr<-sumcorr+prb$V3[i]
        cumcorr[i]<-sumcorr
}
cumcorr<-cumcorr/sumcorr

plot(prb$V1, cumerror, type="h", lwd=3, col="orange3", axes=F, ylab="Cumulative percent errors in observed positions[%]", xlab="prb values", main="", ylim=c(0, max(cumerror)+0.1))
lines(prb$V1+0.2, cumcorr, type="h", col="darkolivegreen4", lwd=3)
axis(1, at=seq(-5, 40, by=5))
axis(2, at=c(0, max(cumerror)/20, max(cumerror)/10, max(cumerror)/2, max(cumerror)/10*8, max(cumerror)/10*9, max(cumerror)), labels=c(0, round(max(cumerror)/20*100, 2), round(max(cumerror)/10*100, 2), round(max(cumerror)/2*100, 2), round(max(cumerror)/10*8*100, 2), round(max(cumerror)/10*9*100, 2), round(max(cumerror)*100, 2)))

text(prb$V1, cumerror, adj=c(-.5,.7), labels=prb$V2, cex=0.65, srt=90)
text(prb$V1+0.2, cumcorr, adj=c(-.5,.7), labels=prb$V3, cex=0.65, srt=90)

abline(h=0.8, lty=3, col="gray")
abline(h=0.9, lty=3, col="gray")
abline(h=1, lty=3, col="gray")
abline(h=0.1, lty=2, col="gray")
abline(h=0.05, lty=2, col="gray")

legend("topleft", legend=c("Cumulative percent of observed errors", "Number of observed errors"), bty="n", fill=c("orange3", "black"))


############# Error by PRB -- read 2 ###############

supp<-(prb_2$V3+prb_2$V2)
rel<-prb_2$V2/supp
rel[is.nan(rel)]<-0

plot(prb_2$V1, rel, type="h", lwd=3, col="orange3", axes=F, ylab="Percent errors of observed positions [%]", xlab="prb (Qs)", main="Errors in mapped reads at uniquely mapped positions, excluding ambiguous calls", ylim=c(0, max(rel)+0.05))
axis(1, at=seq(-5, 40, by=5))
axis(2, at=c(0, max(rel)/2, max(rel)), labels=c(0, round(max(rel)/2*100, 2), round(max(rel)*100, 2)))
axis(4, at=c(0, max(rel)/2, max(rel)), labels=c("0", max(supp[1:45])/2, max(supp[1:45])))

text(prb_2$V1, rel, adj=c(-.5,.7), labels=prb_2$V2, cex=0.65, srt=90, col="black")

lines(seq(-5,40), supp/max(supp[1:45])*max(rel),lty=2, col="firebrick", lwd=2)

legend("topright", legend=c("Percent of observed errors", "Number of observed errors", "Support distribution"), bty="n", fill=c("orange3", "black", "firebrick"))


######### Error by PRB (cumulative) --read 2 ########

cumerror<-c()
sumerror<-0
for (i in 1:length(rel)) {
	sumerror<-sumerror+prb_2$V2[i]
	cumerror[i]<-sumerror
}
cumerror<-cumerror/sumerror

cumcorr<-c()
sumcorr<-0
for (i in 1:length(rel)) {
        sumcorr<-sumcorr+prb_2$V3[i]
        cumcorr[i]<-sumcorr
}
cumcorr<-cumcorr/sumcorr

plot(prb_2$V1, cumerror, type="h", lwd=3, col="orange3", axes=F, ylab="Cumulative percent errors in observed positions[%]", xlab="prb values", main="", ylim=c(0, max(cumerror)+0.1))
lines(prb_2$V1+0.2, cumcorr, type="h", col="darkolivegreen4", lwd=3)
axis(1, at=seq(-5, 40, by=5))
axis(2, at=c(0, max(cumerror)/20, max(cumerror)/10, max(cumerror)/2, max(cumerror)/10*8, max(cumerror)/10*9, max(cumerror)), labels=c(0, round(max(cumerror)/20*100, 2), round(max(cumerror)/10*100, 2), round(max(cumerror)/2*100, 2), round(max(cumerror)/10*8*100, 2), round(max(cumerror)/10*9*100, 2), round(max(cumerror)*100, 2)))

text(prb_2$V1, cumerror, adj=c(-.5,.7), labels=prb_2$V2, cex=0.65, srt=90)
text(prb_2$V1+0.2, cumcorr, adj=c(-.5,.7), labels=prb_2$V3, cex=0.65, srt=90)

abline(h=0.8, lty=3, col="gray")
abline(h=0.9, lty=3, col="gray")
abline(h=1, lty=3, col="gray")
abline(h=0.1, lty=2, col="gray")
abline(h=0.05, lty=2, col="gray")

legend("topleft", legend=c("Cumulative percent of observed errors", "Number of observed errors"), bty="n", fill=c("orange3", "black"))


################ Layout for the second slide ##################################################################################################



############# Error by qCal -- Read 1 ###############

supp<-(qCal$V3+qCal$V2)
rel<-qCal$V2/supp
rel[is.nan(rel)]<-0

plot(qCal$V1, rel, type="h", lwd=3, col="orange3", axes=F, ylab="Percent errors of observed positions [%]", xlab="qCal (Calibrated Qs)", main="Errors in mapped reads at uniquely mapped positions, excluding ambiguous calls", ylim=c(0, max(rel)+0.05))
axis(1, at=seq(-5, 40, by=5))
axis(2, at=c(0, max(rel)/2, max(rel)), labels=c(0, round(max(rel)/2*100, 2), round(max(rel)*100, 2)))
axis(4, at=c(0, max(rel)/2, max(rel)), labels=c("0", max(supp[1:46])/2, max(supp[1:46])))

text(qCal$V1, rel, adj=c(-.5,.7), labels=qCal$V2, cex=0.65, srt=90, col="black")

lines(seq(-5,40), supp/max(supp[1:46])*max(rel),lty=2, col="firebrick", lwd=2)


legend("topright", legend=c("Percent of observed errors", "Number of observed errors", "Support distribution"), bty="n", fill=c("orange3", "black", "firebrick"))


######### Error by qCal (cumulative) -- Read 1 ########

cumerror<-c()
sumerror<-0
for (i in 1:length(rel)) {
        sumerror<-sumerror+qCal$V2[i]
        cumerror[i]<-sumerror
}
cumerror<-cumerror/sumerror

cumcorr<-c()
sumcorr<-0
for (i in 1:length(rel)) {
        sumcorr<-sumcorr+qCal$V3[i]
        cumcorr[i]<-sumcorr
}
cumcorr<-cumcorr/sumcorr


plot(qCal$V1, cumerror, type="h", lwd=3, col="orange3", axes=F, ylab="Cumulative percent errors in observed positions[%]", xlab="qCal (Calibrated Qs)", main="", ylim=c(0, max(cumerror)+0.1))
lines(qCal$V1+0.2, cumcorr, type="h", col="darkolivegreen4", lwd=3)
axis(1, at=seq(-5, 40, by=5))
axis(2, at=c(0, max(cumerror)/20, max(cumerror)/10, max(cumerror)/2, max(cumerror)/10*8, max(cumerror)/10*9, max(cumerror)), labels=c(0, round(max(cumerror)/20*100, 2), round(max(cumerror)/10*100, 2), round(max(cumerror)/2*100, 2), round(max(cumerror)/10*8*100, 2), round(max(cumerror)/10*9*100, 2), round(max(cumerror)*100, 2)))

text(qCal$V1, cumerror, adj=c(-.5,.7), labels=qCal$V2, cex=0.65, srt=90)
text(qCal$V1+0.2, cumcorr, adj=c(-.5,.7), labels=qCal$V3, cex=0.65, srt=90)

abline(h=0.8, lty=3, col="gray")
abline(h=0.9, lty=3, col="gray")
abline(h=1, lty=3, col="gray")
abline(h=0.1, lty=2, col="gray")
abline(h=0.05, lty=2, col="gray")

legend("topleft", legend=c("Cumulative percent of observed errors", "Number of observed errors"), bty="n", fill=c("orange3", "black"))

############# Error by qCal -- Read 2 ###############

supp<-(qCal_2$V3+qCal_2$V2)
rel<-qCal_2$V2/supp
rel[is.nan(rel)]<-0

plot(qCal_2$V1, rel, type="h", lwd=3, col="orange3", axes=F, ylab="Percent errors of observed positions [%]", xlab="qCal (Calibrated Qs)", main="Errors in mapped reads at uniquely mapped positions, excluding ambiguous calls", ylim=c(0, max(rel)+0.05))
axis(1, at=seq(-5, 40, by=5))
axis(2, at=c(0, max(rel)/2, max(rel)), labels=c(0, round(max(rel)/2*100, 2), round(max(rel)*100, 2)))
axis(4, at=c(0, max(rel)/2, max(rel)), labels=c("0", max(supp[1:46])/2, max(supp[1:46])))

text(qCal_2$V1, rel, adj=c(-.5,.7), labels=qCal_2$V2, cex=0.65, srt=90, col="black")

lines(seq(-5,40), supp/max(supp[1:46])*max(rel),lty=2, col="firebrick", lwd=2)


legend("topright", legend=c("Percent of observed errors", "Number of observed errors", "Support distribution"), bty="n", fill=c("orange3", "black", "firebrick"))


######### Error by qCal (cumulative) -- Read 2 ########

cumerror<-c()
sumerror<-0
for (i in 1:length(rel)) {
        sumerror<-sumerror+qCal_2$V2[i]
        cumerror[i]<-sumerror
}
cumerror<-cumerror/sumerror

cumcorr<-c()
sumcorr<-0
for (i in 1:length(rel)) {
        sumcorr<-sumcorr+qCal_2$V3[i]
        cumcorr[i]<-sumcorr
}
cumcorr<-cumcorr/sumcorr


plot(qCal_2$V1, cumerror, type="h", lwd=3, col="orange3", axes=F, ylab="Cumulative percent errors in observed positions[%]", xlab="qCal (Calibrated Qs)", main="", ylim=c(0, max(cumerror)+0.1))
lines(qCal_2$V1+0.2, cumcorr, type="h", col="darkolivegreen4", lwd=3)
axis(1, at=seq(-5, 40, by=5))
axis(2, at=c(0, max(cumerror)/20, max(cumerror)/10, max(cumerror)/2, max(cumerror)/10*8, max(cumerror)/10*9, max(cumerror)), labels=c(0, round(max(cumerror)/20*100, 2), round(max(cumerror)/10*100, 2), round(max(cumerror)/2*100, 2), round(max(cumerror)/10*8*100, 2), round(max(cumerror)/10*9*100, 2), round(max(cumerror)*100, 2)))

text(qCal_2$V1, cumerror, adj=c(-.5,.7), labels=qCal_2$V2, cex=0.65, srt=90)
text(qCal_2$V1+0.2, cumcorr, adj=c(-.5,.7), labels=qCal_2$V3, cex=0.65, srt=90)

abline(h=0.8, lty=3, col="gray")
abline(h=0.9, lty=3, col="gray")
abline(h=1, lty=3, col="gray")
abline(h=0.1, lty=2, col="gray")
abline(h=0.05, lty=2, col="gray")

legend("topleft", legend=c("Cumulative percent of observed errors", "Number of observed errors"), bty="n", fill=c("orange3", "black"))





################ Layout for the second slide ##################################################################################################

############# Error by Chastity -- Read 1 ###############

supp<-(chas$V3+chas$V2)
rel<-chas$V2/supp
rel[is.nan(rel)]<-0

plot(chas$V1, rel, type="h", lwd=3, col="orange3", axes=F, ylab="Percent errors of observed positions [%]", xlab="Chastity values", main="Errors in mapped reads at uniquely mapped positions, excluding ambiguous calls", ylim=c(0, max(rel)+0.05))
axis(1, at=seq(50, 100, by=5))
axis(2, at=c(0, max(rel)/2, max(rel)), labels=c(0, round(max(rel)/2*100, 2), round(max(rel)*100, 2)))
axis(4, at=c(0, max(rel)/2, max(rel)), labels=c("0", max(supp[1:51])/2, max(supp[1:51])))

text(chas$V1, rel, adj=c(-.5,.7), labels=chas$V2, cex=0.65, srt=90, col="black")

lines(seq(50,100), supp/max(supp[1:51])*max(rel),lty=2, col="firebrick", lwd=2)


legend("topright", legend=c("Percent of observed errors", "Number of observed errors", "Support distribution"), bty="n", fill=c("orange3", "black", "firebrick"))


######### Error by Chastity (cumulative) -- Read 1 ########

cumerror<-c()
sumerror<-0
for (i in 1:length(rel)) {
        sumerror<-sumerror+chas$V2[i]
        cumerror[i]<-sumerror
}
cumerror<-cumerror/sumerror

cumcorr<-c()
sumcorr<-0
for (i in 1:length(rel)) {
        sumcorr<-sumcorr+chas$V3[i]
        cumcorr[i]<-sumcorr
}
cumcorr<-cumcorr/sumcorr


plot(chas$V1, cumerror, type="h", lwd=3, col="orange3", axes=F, ylab="Cumulative percent errors in observed positions[%]", xlab="Chastity values", main="", ylim=c(0, max(cumerror)+0.1))
lines(chas$V1+0.2, cumcorr, type="h", col="darkolivegreen4", lwd=3)
axis(1, at=seq(50, 100, by=5))
axis(2, at=c(0, max(cumerror)/20, max(cumerror)/10, max(cumerror)/2, max(cumerror)/10*8, max(cumerror)/10*9, max(cumerror)), labels=c(0, round(max(cumerror)/20*100, 2), round(max(cumerror)/10*100, 2), round(max(cumerror)/2*100, 2), round(max(cumerror)/10*8*100, 2), round(max(cumerror)/10*9*100, 2), round(max(cumerror)*100, 2)))

text(chas$V1, cumerror, adj=c(-.5,.7), labels=chas$V2, cex=0.65, srt=90)
text(chas$V1+0.2, cumcorr, adj=c(-.5,.7), labels=chas$V3, cex=0.65, srt=90)

abline(h=0.8, lty=3, col="gray")
abline(h=0.9, lty=3, col="gray")
abline(h=1, lty=3, col="gray")
abline(h=0.1, lty=2, col="gray")
abline(h=0.05, lty=2, col="gray")

legend("topleft", legend=c("Cumulative percent of observed errors", "Number of observed errors"), bty="n", fill=c("orange3", "black"))


############# Error by Chastity -- Read 2 ###############

supp<-(chas_2$V3+chas_2$V2)
rel<-chas_2$V2/supp
rel[is.nan(rel)]<-0

plot(chas_2$V1, rel, type="h", lwd=3, col="orange3", axes=F, ylab="Percent errors of observed positions [%]", xlab="Chastity values", main="Errors in mapped reads at uniquely mapped positions, excluding ambiguous calls", ylim=c(0, max(rel)+0.05))
axis(1, at=seq(50, 100, by=5))
axis(2, at=c(0, max(rel)/2, max(rel)), labels=c(0, round(max(rel)/2*100, 2), round(max(rel)*100, 2)))
axis(4, at=c(0, max(rel)/2, max(rel)), labels=c("0", max(supp[1:51])/2, max(supp[1:51])))

text(chas_2$V1, rel, adj=c(-.5,.7), labels=chas_2$V2, cex=0.65, srt=90, col="black")

lines(seq(50,100), supp/max(supp[1:51])*max(rel),lty=2, col="firebrick", lwd=2)


legend("topright", legend=c("Percent of observed errors", "Number of observed errors", "Support distribution"), bty="n", fill=c("orange3", "black", "firebrick"))


######### Error by Chastity (cumulative) -- Read 2 ########

cumerror<-c()
sumerror<-0
for (i in 1:length(rel)) {
        sumerror<-sumerror+chas_2$V2[i]
        cumerror[i]<-sumerror
}
cumerror<-cumerror/sumerror

cumcorr<-c()
sumcorr<-0
for (i in 1:length(rel)) {
        sumcorr<-sumcorr+chas_2$V3[i]
        cumcorr[i]<-sumcorr
}
cumcorr<-cumcorr/sumcorr


plot(chas_2$V1, cumerror, type="h", lwd=3, col="orange3", axes=F, ylab="Cumulative percent errors in observed positions[%]", xlab="Chastity values", main="", ylim=c(0, max(cumerror)+0.1))
lines(chas_2$V1+0.2, cumcorr, type="h", col="darkolivegreen4", lwd=3)
axis(1, at=seq(50, 100, by=5))
axis(2, at=c(0, max(cumerror)/20, max(cumerror)/10, max(cumerror)/2, max(cumerror)/10*8, max(cumerror)/10*9, max(cumerror)), labels=c(0, round(max(cumerror)/20*100, 2), round(max(cumerror)/10*100, 2), round(max(cumerror)/2*100, 2), round(max(cumerror)/10*8*100, 2), round(max(cumerror)/10*9*100, 2), round(max(cumerror)*100, 2)))

text(chas_2$V1, cumerror, adj=c(-.5,.7), labels=chas_2$V2, cex=0.65, srt=90)
text(chas_2$V1+0.2, cumcorr, adj=c(-.5,.7), labels=chas_2$V3, cex=0.65, srt=90)

abline(h=0.8, lty=3, col="gray")
abline(h=0.9, lty=3, col="gray")
abline(h=1, lty=3, col="gray")
abline(h=0.1, lty=2, col="gray")
abline(h=0.05, lty=2, col="gray")

legend("topleft", legend=c("Cumulative percent of observed errors", "Number of observed errors"), bty="n", fill=c("orange3", "black"))



#### Layout for the fourth slide ############################################################################################################################

layoutmat<-matrix(data=c(1,2,3,4,5,6,7,7,7,8,8,8), nrow=2, ncol=6, byrow=TRUE)
layout(layoutmat)

###### Mismatch per read distribution -- Read 1 ###
mismatch_percent<-(mismatch_dist$V2/sum(mismatch_dist$V2)*100)

plot(mismatch_dist$V1, mismatch_percent, type="b", lwd=1, col="orange3", axes=F, ylab="Fraction [%]", xlab="Mismatches per read", main="Mismatches per Read", ylim=c(0, (1.05 * max(mismatch_percent))), xlim=c(0, 4))

axis(1, at=c(0,1,2,3,4), labels=c("0","1","2","3","4"))
axis(2, at=c(0, round(max(mismatch_percent)/2, 0), round(max(mismatch_percent), 0)), labels=c("0", round(max(mismatch_percent)/2, 0), round(max(mismatch_percent), 0)))

text(mismatch_dist$V1, mismatch_percent, adj=c(0.5, -2), labels=mismatch_dist$V2, cex=0.65, srt=0, col="black")


##### Hits per read distribution -- Read 1 #######

real_origins<-(hits_dist$V2/hits_dist$V1)
log_count<-log(real_origins)

plot(hits_dist$V1, log_count, type="h", lwd=1, col="orange3", axes=F, ylab="Log count", xlab="Hits per read", main="Number of hits in genome per read", ylim=c(0, (1.05 * max(log_count))), xlim=c(1, 50))

axis(1, at=seq(1, 51, by=10))
axis(2, at=c(0, max(log_count)/2, max(log_count)), labels=c("0", round(max(log_count)/2, 1), round(max(log_count), 1)))


############# Read length -- Read 1 ##############
read_length_percent<-(read_length_dist$V2/sum(read_length_dist$V2)*100)

pie( read_length_dist$V2, labels=read_length_dist$V1, lwd=3, col="orange3",  main="Read Length Distribution" )


###### Mismatch per read distribution -- Read 2 ###
mismatch_percent<-(mismatch_dist_2$V2/sum(mismatch_dist_2$V2)*100)

plot(mismatch_dist_2$V1, mismatch_percent, type="b", lwd=1, col="orange3", axes=F, ylab="Fraction [%]", xlab="Mismatches", main="Mismatches per Read", ylim=c(0, (1.05 * max(mismatch_percent))), xlim=c(0, 4))

axis(1, at=c(0,1,2,3,4), labels=c("0","1","2","3","4"))
axis(2, at=c(0, round(max(mismatch_percent)/2, 0), round(max(mismatch_percent), 0)), labels=c("0", round(max(mismatch_percent)/2, 0), round(max(mismatch_percent), 0)))

text(mismatch_dist_2$V1, mismatch_percent, adj=c(0.5, -2), labels=mismatch_dist_2$V2, cex=0.65, srt=0, col="black")

##### Hits per read distribution -- Read 2 #######

real_origins<-(hits_dist_2$V2/hits_dist_2$V1)
log_count<-log(real_origins)

plot(hits_dist_2$V1, log_count, type="h", lwd=1, col="orange3", axes=F, ylab="Log count", xlab="Hits per read", main="Number of hits in genome per read", ylim=c(0, (1.05 * max(log_count))), xlim=c(1, 50))

axis(1, at=seq(1, 51, by=10))
axis(2, at=c(0, max(log_count)/2, max(log_count)), labels=c("0", round(max(log_count)/2, 1), round(max(log_count), 1)))


############# Read length -- Read 2 ##############
read_length_percent<-(read_length_dist_2$V2/sum(read_length_dist_2$V2)*100)

pie( read_length_dist_2$V2, labels=read_length_dist_2$V1, lwd=3, col="orange3",  main="Read Length Distribution" )
#######################################################
######## Preprocess Error plot

supp1<-(position$V3+position$V2)
rel1<-position$V2/supp1

supp2<-(position_2$V3+position_2$V2)
rel2<-position_2$V2/supp2

xmax<-max(c(max(rel1),max(rel2)))

########### Error by position -- Read 1 ###############

plot(position$V1, rel1, type="h", lwd=3, col="orange3", axes=F, ylab="Percent errors [%]", xlab="Position in read", main="Errors in read positions", ylim=c(0, xmax+0.005), xlim=c(1, length(position$V1)))

axis(1, at=seq(1, length(position$V1)))
axis(2, at=c(0, xmax/2, xmax), labels=c("0", round(xmax/2*100, 2), round(xmax*100, 2)))
axis(4, at=c(0, xmax/2, xmax), labels=c("0", max(supp1[1:36])/2, max(supp1[1:36])))

text(position$V1, rel1, adj=c(-.5,.7), labels=position$V2, cex=0.65, srt=90, col="black")

lines(seq(1,length(position$V1)), supp1/max(supp1)*max(rel1),lty=2, col="firebrick", lwd=2)

legend("left", legend=c("Percent of observed errors", "Number of observed errors", "Support distribution"), bty="n", fill=c("orange3", "black", "firebrick"))

########### Error by position -- Read 2 ###############

plot(position_2$V1, rel2, type="h", lwd=3, col="orange3", axes=F, ylab="Percent errors [%]", xlab="Position in read", main="Errors in read positions", ylim=c(0, xmax+0.005), xlim=c(1, length(position_2$V1)))
axis(1, at=seq(1, length(position_2$V1)))
axis(2, at=c(0, xmax/2, xmax), labels=c("0", round(xmax/2*100, 2), round(xmax*100, 2)))
axis(4, at=c(0, xmax/2, xmax), labels=c("0", max(supp2[1:36])/2, max(supp2[1:36])))

text(position_2$V1, rel2, adj=c(-.5,.7), labels=position_2$V2, cex=0.65, srt=90, col="black")

lines(seq(1,length(position_2$V1)), supp2/max(supp2)*max(rel2),lty=2, col="firebrick", lwd=2)

legend("left", legend=c("Percent of observed errors", "Number of observed errors", "Support distribution"), bty="n", fill=c("orange3", "black", "firebrick"))


##############################################################################################################

# Thanks for flying with R, enjoy your stay in the world of colorful R plots and good-bye:

dev.off()

