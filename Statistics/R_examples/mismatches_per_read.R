pdf(file=("figure01_C.pdf"), width=8, height=6)

errornumread<-read.table("error_num_per_read.txt");

########### Plot error per reads ######################################################

errornumreadperc<-round((errornumread$V2/sum(errornumread$V2))*100, digits=2)

barplot(errornumread$V2/sum(errornumread$V2)*100, col="steelblue4", width=0.4, xlim=c(0,4), ylim=c(0,80), axes=F, space=0.8, xlab="Number of mismatches", ylab="Fraction of reads [%]")

axis(1)
axis(2, las=1, at=c(0,10,20,30,40,50,60,70,80), labels=c(0,10,20,30,40,50,60,70,80))

text(c(0,1,2,3,4), errornumread$V2/sum(errornumread$V2)*100+4, adj=c(-.5,.7), labels=round(errornumread$V2/sum(errornumread$V2)*100, digits=2), cex=0.8, srt=0, col="black")

dev.off()
