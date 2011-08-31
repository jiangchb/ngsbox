#  wc -l */* | grep -v total | sed -e "s/length_//g" | sed -e "s/\/reads_0.fl//g" | sort -k2n | sed -e "s/^ *//g" | sed -e "s/ /\t/g" > ../length.dist

pdf(file="length.dist.pdf", width=10, height=8)

data2<-read.table("length.dist.2")
data4<-read.table("length.dist.4")
data6<-read.table("length.dist.6")
data8<-read.table("length.dist.8")

plot(data6$V2, data6$V1, type="l", col="red", lwd=2, xlab="sRNA length", ylab="occurrence", main="length distribution sRNAs", axes=F, xlim=c(5, 36))
axis(2)
labels=c(5, 17, 21, 24, 29, 36)
axis(1, labels=labels, at=labels)

lines(data2$V2, data2$V1, col="green", lwd=2)
lines(data4$V2, data4$V1, col="blue", lwd=2)
lines(data8$V2, data8$V1, col="yellow", lwd=2)


dev.off()
