pdf(file="rnv_perlegen.pdf", width=18, height=10)

pos0<-read.table("perlegen_rnv_validation_Tsu_Int_600_pos_0.txt")
posm1<-read.table("perlegen_rnv_validation_Tsu_Int_600_pos_-1.txt")
pos1<-read.table("perlegen_rnv_validation_Tsu_Int_600_pos_1.txt")

layoutmat<-matrix(data=c(1,2,3,4,5,6), nrow=2, ncol=3)
layout(layoutmat)

ymax=300
coln="steelblue4"
colu="goldenrod1"

hist2<-hist(posm1$V4, main="", xlab="1bp upstream of RNV position", xlim=c(0.5, 1.0), ylim=c(0,ymax), breaks=20, col=coln)
hist1<-hist(posm1$V3, main="", xlab="", xlim=c(0.5, 1.0), ylim=c(0,ymax), breaks=hist2$breaks, col=coln)
x<-hist1$counts
y<-hist2$counts


hist3<-hist(pos0$V3, main="", xlab="", xlim=c(0.5, 1.0), ylim=c(0,ymax), breaks=20, col=coln)
hist4<-hist(pos0$V4, main="", xlab="RNV position", xlim=c(0.5, 1.0), ylim=c(0,ymax), breaks=hist3$breaks, col=colu)
x<-hist3$counts
y<-hist4$counts



hist5<-hist(pos1$V3, main="", xlab="", xlim=c(0.5, 1.0), ylim=c(0,ymax), breaks=20, col=coln)
hist6<-hist(pos1$V4, main="", xlab="1bp downstream of RNV position", xlim=c(0.5, 1.0), ylim=c(0,ymax), breaks=hist5$breaks, col=coln)
x<-hist5$counts
y<-hist6$counts



chisq.test(y, p=x, rescale.p=T)




dev.off()



