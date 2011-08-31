pdf(file="zellerPR_solexaPR_bp_overlap.pdf", width=10, height=15)


layoutmat<-matrix(data=c(1,2,3), nrow=3, ncol=1)
layout(layoutmat)

real100<-read.table("bur_100.txt")
real250<-read.table("bur_250.txt")
real8<-read.table("bur_8.txt")

d100_1<-read.table("per_bur_100_chr1.txt")
d100_2<-read.table("per_bur_100_chr2.txt")
d100_3<-read.table("per_bur_100_chr3.txt")
d100_4<-read.table("per_bur_100_chr4.txt")
d100_5<-read.table("per_bur_100_chr5.txt")

d250_1<-read.table("per_bur_250_chr1.txt")
d250_2<-read.table("per_bur_250_chr2.txt")
d250_3<-read.table("per_bur_250_chr3.txt")
d250_4<-read.table("per_bur_250_chr4.txt")
d250_5<-read.table("per_bur_250_chr5.txt")

d8_1<-read.table("per_bur_8_chr1.txt")
d8_2<-read.table("per_bur_8_chr2.txt")
d8_3<-read.table("per_bur_8_chr3.txt")
d8_4<-read.table("per_bur_8_chr4.txt")
d8_5<-read.table("per_bur_8_chr5.txt")

d100<-c()
d250<-c()
d8<-c()

for (i in 1:1000) { 
	d100[i] <- d100_1$V2[i] + d100_2$V2[i] + d100_3$V2[i] + d100_4$V2[i] + d100_5$V2[i]
	d250[i] <- d250_1$V2[i] + d250_2$V2[i] + d250_3$V2[i] + d250_4$V2[i] + d250_5$V2[i]
	d8[i] <- d8_1$V2[i] + d8_2$V2[i] + d8_3$V2[i] + d8_4$V2[i] + d8_5$V2[i]
}


d100r <- d100_1$V3[1] + d100_2$V3[1] + d100_3$V3[1] + d100_4$V3[1] + d100_5$V3[1]
d250r <- d250_1$V3[1] + d250_2$V3[1] + d250_3$V3[1] + d250_4$V3[1] + d250_5$V3[1]
d8r <- d8_1$V3[1] + d8_2$V3[1] + d8_3$V3[1] + d8_4$V3[1] + d8_5$V3[1]

ylimit<-120
hist(d100[]/d100r, breaks=50, xlim=c(0,1.0), ylim=c(0,ylimit), xlab="Permutation of SolexaPR and Zeller PR, length >= 100", main="Histogram of PR overlap")
arrows(real100$V3, 0, real100$V3, 20, col="darkolivegreen4", lwd=3, code=1)
hist(d250[]/d250r, breaks=50, xlim=c(0,1.0), ylim=c(0,ylimit), xlab="Permutation of SolexaPR and Zeller PR, length >= 250", main="Histogram of PR overlap")
arrows(real250$V3, 0, real250$V3, 20, col="darkolivegreen4", lwd=3, code=1)
hist(d8[]/d8r, breaks=50, xlim=c(0,1.0), ylim=c(0,ylimit), xlab="Permutation of SolexaPR and Zeller PR, length >= 8", main="Histogram of PR overlap")
arrows(real8$V3, 0, real8$V3, 20, col="darkolivegreen4", lwd=3, code=1)



dev.off()
