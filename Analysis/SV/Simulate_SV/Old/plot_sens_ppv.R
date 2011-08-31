
min=c(		4,	25,	49,	73,	101,	151,	201)
max=c(		24,	48,	72,	100,	150,	200,	1000)
at=c(		14,	37,	61,	87,	125,	175,	600)

# Fertisch
insert_tp=c(0, 		5, 	467, 	441, 	9, 	0, 	0)
insert_fn=c(1000, 	995, 	533, 	559, 	991, 	1000, 1000)
insert_fp=c(0, 		0, 	3, 	4, 	0, 	0, 	0)

# Fertisch
deletion_tp=c(0,	596,	1000,	1000,	1000,	1000,	1000)
deletion_fn=c(1000,	404,	0,	0,	0,	0,	0)
deletion_fp=c(0, 	0,	0,	0, 	0, 	0, 	0)

# Fertisch
inversion_tp=c(0,	2,	336,	944,	1000,	1000,	1000)
inversion_fn=c(1000,	998,	644,	56,	0,	0, 	0)
inversion_fp=c(0,	0,	0,	0,	0,	0,	0)

insert_ppv=insert_tp / (insert_tp + insert_fp)
insert_sens=insert_tp / (insert_tp + insert_fn)

deletion_ppv=deletion_tp / (deletion_tp + deletion_fp)
deletion_sens=deletion_tp / (deletion_tp + deletion_fn)

inversion_ppv=inversion_tp / (inversion_tp + inversion_fp)
inversion_sens=inversion_tp / (inversion_tp + inversion_fn)
  
# plotting

pdf(file="SV_simulation.pdf", width=11, height=8)

layoutmat<-matrix(data=c(1,2), nrow=2, ncol=1)
layout(layoutmat)

red="#F90000"
darkred="#CD2626"
green="#86CC00"
darkgreen="#6E8B3D"
yellow="#FFC125"
darkyellow="#FF8B00"

# ppv

plot(at, insert_ppv, type="n", ylim=c(0.95,1), axes=F, ylab="", xlab="Size range in bp", main="Positive predictive value")
axis(1, at=at, labels=c("4-24", "25-48", "49-72", "73-100", "101-150", "151-200", "201-1000"), las=2)
axis(2, las=1)

pch=23

lines(at, inversion_ppv, type="l", lwd=1, col=darkyellow)
points(at, inversion_ppv, pch=24, col=darkyellow, bg=darkyellow)

lines(at, insert_ppv, type="l", lwd=1, col=darkgreen)
points(at, insert_ppv, pch=23, col=darkgreen, bg=darkgreen)


lines(at, deletion_ppv, type="c", lwd=1, col=darkred)
points(at, deletion_ppv, pch=25, col=darkred, bg=darkred)

# sensitivity

plot(at, insert_ppv, type="n", ylim=c(0,1), axes=F, ylab="", xlab="Size range in bp", main="Sensitivity")
axis(1, at=at, labels=c("4-24", "25-48", "49-72", "73-100", "101-150", "151-200", "201-1000"), las=2)
axis(2, las=1)

lines(at, inversion_sens, type="l", lwd=1, col=darkyellow)
points(at, inversion_sens, pch=24, col=darkyellow, bg=darkyellow)

lines(at, insert_sens, type="l", lwd=1, col=darkgreen)
points(at, insert_sens, pch=23, col=darkgreen, bg=darkgreen)

lines(at, deletion_sens, type="c", lwd=1, col=darkred)
points(at, deletion_sens, pch=25, col=darkred, bg=darkred)




dev.off()

