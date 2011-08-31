args <- commandArgs()

########## Set colors ############

black<-"black"
orange<-"peru"
red<-"firebrick"
green<-"darkolivegreen4"
gray<-"gray"

########## Read in files ############

AT_TA<-read.table(args[5])
CG_GC<-read.table(args[6])
AC_CA_GT_TG<-read.table(args[7])
AG_GA_CT_TC<-read.table(args[8])
AT<-read.table(args[9])
CG<-read.table(args[10])
name<-args[11]

pdf(file=paste(name, ".pdf", sep=""), width=15, height=6)


### Layout
#layoutmat<-matrix(data=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18), nrow=3, ncol=6)
layoutmat<-matrix(data=c(1,2,3,4,5,6,7,8,9,10,11,12), nrow=2, ncol=6)
layout(layoutmat)

## get max y 
maxy=max(max(log10(AT_TA$V2)), max(log10(CG_GC$V2)), max(log10(AC_CA_GT_TG$V2)), max(log10(AG_GA_CT_TC$V2)), max(log10(AT$V2)), max(log10(CG$V2)) )


### AT/TA
percent_gapped<-(AT_TA$V3 / AT_TA$V2)
plot(AT_TA$V1, log10(AT_TA$V2), type="h", lwd=3, col=gray, axes=T, ylab="log10(count)", xlab="Repeat instances", main="AT/TA Dipolymer", xlim=c(min(AT_TA$V1),18), ylim=c(0, maxy))
#plot(AT_TA$V1, percent_gapped, type="h", lwd=3, col=orange, axes=T, ylab="Percent gapped [%]", xlab="Repeat instances", main="", xlim=c(min(AT_TA$V1),18), ylim=c(0, max(percent_gapped)+0.05))
plot(AT_TA$V1, AT_TA$V4, type="h", lwd=3, col=green, axes=T, ylab="Average Coverage", xlab="Repeat instances", main="", xlim=c(min(AT_TA$V1),18), ylim=c(0, 60))


percent_gapped<-(CG_GC$V3 / CG_GC$V2)
plot(CG_GC$V1, log10(CG_GC$V2), type="h", lwd=3, col=gray, axes=T, ylab="log10(count)", xlab="Repeat instances", main="CG/GC Dipolymer", xlim=c(min(AT_TA$V1),18), ylim=c(0, maxy))
#plot(CG_GC$V1, percent_gapped, type="h", lwd=3, col=orange, axes=T, ylab="Percent gapped [%]", xlab="Repeat instances", main="", xlim=c(min(AT_TA$V1),18), ylim=c(0, max(percent_gapped)+0.05))
plot(CG_GC$V1, CG_GC$V4, type="h", lwd=3, col=green, axes=T, ylab="Average Coverage", xlab="Repeat instances", main="", xlim=c(min(AT_TA$V1),18), ylim=c(0, 60))

percent_gapped<-(AC_CA_GT_TG$V3 / AC_CA_GT_TG$V2)
plot(AC_CA_GT_TG$V1, log10(AC_CA_GT_TG$V2), type="h", lwd=3, col=gray, axes=T, ylab="log10(count)", xlab="Repeat instances", main="AC/CA/GT/TG Dipolymer", xlim=c(min(AT_TA$V1),18), ylim=c(0, maxy))
#plot(AC_CA_GT_TG$V1, percent_gapped, type="h", lwd=3, col=orange, axes=T, ylab="Percent gapped [%]", xlab="Repeat instances", main="", xlim=c(min(AT_TA$V1),18), ylim=c(0, max(percent_gapped)+0.05))
plot(AC_CA_GT_TG$V1, AC_CA_GT_TG$V4, type="h", lwd=3, col=green, axes=T, ylab="Average Coverage", xlab="Repeat instances", main="", xlim=c(min(AT_TA$V1),18), ylim=c(0, 60))

percent_gapped<-(AG_GA_CT_TC$V3 / AG_GA_CT_TC$V2)
plot(AG_GA_CT_TC$V1, log10(AG_GA_CT_TC$V2), type="h", lwd=3, col=gray, axes=T, ylab="log10(count)", xlab="Repeat instances", main="AG/GA/CT/TC Dipolymer", xlim=c(min(AT_TA$V1),18), ylim=c(0, maxy))
#plot(AG_GA_CT_TC$V1, percent_gapped, type="h", lwd=3, col=orange, axes=T, ylab="Percent gapped [%]", xlab="Repeat instances", main="", xlim=c(min(AT_TA$V1),18), ylim=c(0, max(percent_gapped)+0.05))
plot(AG_GA_CT_TC$V1, AG_GA_CT_TC$V4, type="h", lwd=3, col=green, axes=T, ylab="Average Coverage", xlab="Repeat instances", main="", xlim=c(min(AT_TA$V1),18), ylim=c(0, 60))

percent_gapped<-(AT$V3 / AT$V2)
plot(AT$V1, log10(AT$V2), type="h", lwd=3, col=gray, axes=T, ylab="log10(count)", xlab="Repeat instances", main="A/T Homopolymer", xlim=c(min(AT$V1),18), ylim=c(0, maxy))
#plot(AT$V1, percent_gapped, type="h", lwd=3, col=orange, axes=T, ylab="Percent gapped [%]", xlab="Repeat instances", main="", xlim=c(min(AT$V1),18), ylim=c(0, max(percent_gapped)+0.05))
plot(AT$V1, AT$V4, type="h", lwd=3, col=green, axes=T, ylab="Average Coverage", xlab="Repeat instances", main="", xlim=c(min(AT$V1),18), ylim=c(0, 60))

percent_gapped<-(CG$V3 / CG$V2)
plot(CG$V1, log10(CG$V2), type="h", lwd=3, col=gray, axes=T, ylab="log10(count)", xlab="Repeat instances", main="C/G Homopolymer", xlim=c(min(AT$V1),18), ylim=c(0, maxy))
#plot(CG$V1, percent_gapped, type="h", lwd=3, col=orange, axes=T, ylab="Percent gapped [%]", xlab="Repeat instances", main="", xlim=c(min(AT$V1),18), ylim=c(0, max(percent_gapped)+0.05))
plot(CG$V1, CG$V4, type="h", lwd=3, col=green, axes=T, ylab="Average Coverage", xlab="Repeat instances", main="", xlim=c(min(AT$V1),18), ylim=c(0, 60))


dev.off()

