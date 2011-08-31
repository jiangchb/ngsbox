pdf(file="overlap_hist.pdf", width=10, height=6)

pos<-read.table("overlap_hist.txt")

ymax=150
coln="steelblue4"
colu="goldenrod1"

hist<-hist(pos$V1, main="", xlab="ChipSeq overlapping TFBS-Predictions", xlim=c(0, 100), ylim=c(0,ymax), breaks=100, col=coln)
x<-hist$counts

arrows(82, 0, 82, 20, col="darkolivegreen4", lwd=3, code=1)

dev.off()



