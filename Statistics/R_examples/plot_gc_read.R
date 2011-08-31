pdf(file="gc_content_read.pdf", width=10, height=8)

data<-read.table("rinput.txt")

plot(data$V2, type="l", col="firebrick", lwd=2, ylim=c(0.15, 0.35), axes=F, xlab="GC content by base", ylab="GC content[%]")
lines(data$V3, type="l", lwd = 2,  col="darkblue")
lines(data$V4, type="l", lwd = 2,  col="darkgreen")
lines(data$V5, type="l", lwd = 2,  col="orange3")
axis(2)
axis(1, at=c(1,38), labels=c("1","38"))

legend(15, 0.35, legend=c("A", "C", "G", "T"), col=c("firebrick", "darkblue", "darkgreen", "orange3"), lwd=2, bty="n")

dev.off()

