pdf(file="zellerPR_solexaPR_bp_overlap_bur_and_tsu.pdf", width=10, height=10)


layoutmat<-matrix(data=c(1,2), nrow=2, ncol=1)
layout(layoutmat)

real8<-read.table("bur_8.txt")
d8_1<-read.table("per_bur_8_chr1.txt")
d8_2<-read.table("per_bur_8_chr2.txt")
d8_3<-read.table("per_bur_8_chr3.txt")
d8_4<-read.table("per_bur_8_chr4.txt")
d8_5<-read.table("per_bur_8_chr5.txt")

d8<-c()

for (i in 1:1000) { 
	d100[i] <- d100_1$V2[i] + d100_2$V2[i] + d100_3$V2[i] + d100_4$V2[i] + d100_5$V2[i]
	d250[i] <- d250_1$V2[i] + d250_2$V2[i] + d250_3$V2[i] + d250_4$V2[i] + d250_5$V2[i]
	d8[i] <- d8_1$V2[i] + d8_2$V2[i] + d8_3$V2[i] + d8_4$V2[i] + d8_5$V2[i]
}


d8r <- d8_1$V3[1] + d8_2$V3[1] + d8_3$V3[1] + d8_4$V3[1] + d8_5$V3[1]

ylimit<-250
burhist<-hist(d8[]/d8r, breaks=35, xlim=c(0,1.0), ylim=c(0,ylimit), xlab="Permutation of SolexaPR and Zeller PR, length >= 8", main="Histogram of PR overlap")
arrows(real8$V3, 0, real8$V3, 20, col="darkolivegreen4", lwd=3, code=1)



real8<-read.table("Tsu_8_x.txt")

d8_1<-read.table("per_Tsu_8_chr1.txt")
d8_2<-read.table("per_Tsu_8_chr2.txt")
d8_3<-read.table("per_Tsu_8_chr3.txt")
d8_4<-read.table("per_Tsu_8_chr4.txt")
d8_5<-read.table("per_Tsu_8_chr5.txt")

d8<-c()

for (i in 1:1000) {
        d8[i] <- d8_1$V2[i] + d8_2$V2[i] + d8_3$V2[i] + d8_4$V2[i] + d8_5$V2[i]
}

d8r <- d8_1$V3[1] + d8_2$V3[1] + d8_3$V3[1] + d8_4$V3[1] + d8_5$V3[1]

tsuhist<-hist(d8[]/d8r, breaks=35, xlim=c(0,1.0), ylim=c(0,ylimit), xlab="Permutation of SolexaPR and Zeller PR, length >= 8", main="Histogram of PR overlap")
arrows(real8$V3, 0, real8$V3, 20, col="darkolivegreen4", lwd=3, code=1)







dev.off()
