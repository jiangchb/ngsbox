pdf(file=("figure01.pdf"), width=16, height=15)

basecomp<-read.table("error_basecomparison.txt");
errornumread<-read.table("error_num_per_read.txt");
errorpos<-read.table("error_positionwise.txt");
prbchr1<-read.table("R_error_prb_plot.chr1.txt")
prbchr2<-read.table("R_error_prb_plot.chr2.txt")
prbchr3<-read.table("R_error_prb_plot.chr3.txt")
prbchr4<-read.table("R_error_prb_plot.chr4.txt")
prbchr5<-read.table("R_error_prb_plot.chr5.txt")

chaschr1<-read.table("R_error_chas_plot.chr1.txt")
chaschr2<-read.table("R_error_chas_plot.chr2.txt")
chaschr3<-read.table("R_error_chas_plot.chr3.txt")
chaschr4<-read.table("R_error_chas_plot.chr4.txt")
chaschr5<-read.table("R_error_chas_plot.chr5.txt")



prberror1<-prbchr1$V2[prbchr1$V1 > 0] + prbchr2$V2[prbchr2$V1 > 0] + prbchr3$V2[prbchr3$V1 > 0] + prbchr4$V2[prbchr4$V1 > 0] + prbchr5$V2[prbchr5$V1 > 0]
prberror0<-sum(prbchr1$V2[prbchr1$V1 <= 0]) + sum(prbchr2$V2[prbchr2$V1 <= 0]) + sum(prbchr3$V2[prbchr3$V1 <= 0]) + sum(prbchr4$V2[prbchr4$V1 <= 0]) + sum(prbchr5$V2[prbchr5$V1 <= 0])
prberror<-c(prberror0, prberror1)

prbsupport1<-prbchr1$V3[prbchr1$V1 > 0] + prbchr2$V3[prbchr2$V1 > 0] + prbchr3$V3[prbchr3$V1 > 0] + prbchr4$V3[prbchr4$V1 > 0] + prbchr5$V3[prbchr5$V1 > 0]
prbsupport0<-sum(prbchr1$V3[prbchr1$V1 <= 0]) + sum(prbchr2$V3[prbchr2$V1 <= 0]) + sum(prbchr3$V3[prbchr3$V1 <= 0]) + sum(prbchr4$V3[prbchr4$V1 <= 0]) + sum(prbchr5$V3[prbchr5$V1 <= 0])
prbsupport<-c(prbsupport0, prbsupport1)

chaserror <- chaschr1$V2 + chaschr2$V2 + chaschr3$V2 + chaschr4$V2 + chaschr5$V2
chassupport <- chaschr1$V3 + chaschr2$V3 + chaschr3$V3 + chaschr4$V3 + chaschr5$V3

########### Layout ####################################################################

layoutmat<-matrix(data=c(1,2), nrow=2, ncol=1)
layout(layoutmat)

########### Error by position ##########################################################

supp<-errorpos$V2
rel<-errorpos$V4

plot(errorpos$V1, rel, type="h", ylim=c(0, max(rel)+0.01), xlim=c(1, length(errorpos$V1)), axes=F, ylab="Observed errors [%]", xlab="Position in read", main="Errors by read length", col="darkolivegreen4", lwd=3)
#points(errorpos$V1, rel, pch=4, lwd=3, col="darkolivegreen4")
axis(1, at=seq(1, length(errorpos$V1)), cex.axis=0.8)
#axis(2, at=c(0, 0.01, 0.02, 0.03, 0.04, 0.05, max(rel)), labels=c("0", "1", "2", "3", "4", "5", round(max(rel)*100, 2)), las=1)
axis(2, at=c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06), labels=c("0", "1", "2", "3", "4", "5", "6"), las=1)

abline(h=0.01, lty=2, col="darkgrey")
abline(h=0.02, lty=2, col="darkgrey")
abline(h=0.03, lty=2, col="darkgrey")
abline(h=0.04, lty=2, col="darkgrey")
abline(h=0.05, lty=2, col="darkgrey")
abline(h=0.06, lty=2, col="darkgrey")

#abline(h=max(rel), lty=2, col="darkgrey")

#text(errorpos$V1, rel, adj=c(-.5,.7), labels=supp, cex=0.65, srt=90, col="black")

#legend("left", legend=c("Percentage of observed errors"), bty="n", fill=c("firebrick"))


############# Error by PRB ###############

supp<-prbsupport
maxy=0.4

plot(seq(0, 40), prberror/prbsupport, type="h", lwd=3, col="darkolivegreen4", axes=F, ylab="Observed errors [%]", xlab="\"prb\"-quality value\nSupport percent for prb=40 = 39.76%", main="Errors by base quality", ylim=c(0, maxy))

text(seq(0,40), prberror/prbsupport, adj=c(-.5,.7), labels=prberror, cex=0.8, srt=70, col="black")
axis(1, at=c(0, seq(5, 40, by=5)), labels=c("<= 0", seq(5, 40, by=5)))
axis(2, las= 1, at=seq(0, 0.4, by=0.1), labels=seq(0, 40, by=10))

lines(seq(0,40), prbsupport/max(prbsupport[1:40]) * 0.35, col="steelblue4", lwd=2)
# 0.35 in der Scale entspricht 0.020489573 = 2.04..% im Wertebereich
# fac=0.0102447865
fac=17.0818591486
axis(4, las=1, at=c(0, 0.01 * fac, 0.02 * fac), labels=c(0, 1, 2))



############# Error by Chas ###############

supp<-chassupport
maxy=1.0

plot(seq(50, 100), chaserror/chassupport, type="h", lwd=3, col="darkolivegreen4", axes=F, ylab="Observed errors [%]", xlab="\"Chastity\"-quality value", main="Errors by base quality", ylim=c(0, maxy))

text(seq(50,100), chaserror/chassupport, adj=c(-.5,.7), labels=chaserror, cex=0.8, srt=70, col="black")

axis(1, at=c(seq(50, 100, by=5)), labels=c(seq(50, 100, by=5)))
axis(2, las= 1, at=seq(0, 1.0, by=0.1), labels=seq(0, 100, by=10))

lines(seq(50, 100), chassupport/max(chassupport), col="steelblue4", lwd=2)

axis(4, las=1, at=c(0,0.5,1), labels=c(0, max(chassupport)/2, max(chassupport)))




########### Plot error per reads ######################################################

#errornumreadperc<-round((errornumread$V2/sum(errornumread$V2))*100, digits=2)

#plot(errornumread$V1, errornumread$V2, type="b", lwd=3, col="orange3", axes=F, ylab="Fraction [%]", xlab="Mismatches/errors per read", main="Mismatches and errors per Read", ylim=c(0, (1.05 * max(errornumread))), xlim=c(0, 4))

#abline(h=errornumread$V2[1], lty=2, col="darkgrey")
#abline(h=errornumread$V2[2], lty=2, col="darkgrey")
#abline(h=errornumread$V2[3], lty=2, col="darkgrey")
#abline(h=errornumread$V2[4], lty=2, col="darkgrey")
#abline(h=errornumread$V2[5], lty=2, col="darkgrey")

#axis(1, at=c(0,1,2,3,4), labels=c("0","1","2","3","4"))
#axis(2, at=c(0, errornumread$V2), labels=c("0", errornumreadperc), las=2)


dev.off()
