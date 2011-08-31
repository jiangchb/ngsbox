data<-read.table("80.inversion_high_quality.plot")

pdf(file="80.inversions.pdf", width=10, height=12)

layoutmat=matrix(data=c(1,2,3,4,5), ncol=1, nrow=5)
layout(layoutmat)

for (chr in 1:5) {

	plot(0,0,type="n", xlim=c(0,31000000), ylim=c(0,80), axes = F, ylab="", xlab=paste("Chromosome ", chr))
	axis(1)
	for (i in 1:223) { 
		if (data$V3[i] == chr) {
			segments(data$V4[i], data$V1[i], data$V5[i], data$V1[i], lwd = 5, col="darkred")
			#text(data$V4[i]-550000, data$V1[i], labels=(data$V5[i]-data$V4[i]+1))
		} 
	}

}


dev.off()


