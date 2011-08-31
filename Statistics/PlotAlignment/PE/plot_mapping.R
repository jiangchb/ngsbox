args <- commandArgs()
mat<-read.table(args[5])
bridgeraw<-read.table(args[6])
format<-args[7]
readlength<-as.numeric(args[8])

chr<-as.numeric(args[9])
start<-as.numeric(args[10])
stop<-as.numeric(args[11])

print(paste(chr, start, stop))

#seq<-read.table(args[8])

xmax<-max(mat$V1)
ymax<-max(mat$V2)

if (format == "pdf") {
	pdf(file=paste("mapping_vis.", chr,"-", start ,"-",stop ,".pdf", sep=""), width=max(15, ((stop-start+1)/2000) * 15), height=ymax/50)
	#pdf(file=paste("mapping_vis.", chr,"-", start ,"-",stop ,".pdf", sep=""), width=20, height=10)
} else {
	png(filename=paste("mapping_vis.", start ,"-",stop ,".png", sep=""))
}
par(mai=c(1.75,0,0,0)) # c(bottom, left, top, right)'

#########################
## PE plotting
plot(stop+readlength, ymax, type="n", axes=F, xlab=paste("Chromosome:", chr, "Position:", start, "-", stop), ylab="", xlim=c(start, stop), ylim=c(0,ymax))

for (i in 1:length(mat$V1)) {
	if (mat$V1[i]+readlength > start & mat$V1[i] < stop) {
		lines(c(mat$V1[i], mat$V1[i]+readlength), c(mat$V2[i], mat$V2[i]), col=as.character(mat$V4[i]), lwd=3)
		if (mat$V3[i] == "D") {
			#lines(lines(c(mat$V1[i]+readlength-(readlength/8), mat$V1[i]+readlength), c(mat$V2[i]+0.05, mat$V2[i]), col=as.character(mat$V4[i]), lwd=2))
			lines(c(mat$V1[i]+readlength-(readlength/8), mat$V1[i]+readlength), c(mat$V2[i]+1.5, mat$V2[i]), col=as.character(mat$V4[i]), lwd=3)
		}
		else {
			#lines(lines(c(mat$V1[i], mat$V1[i]+(readlength/8)), c(mat$V2[i], mat$V2[i]+1.5), col=as.character(mat$V4[i]), lwd=3))
			lines(c(mat$V1[i], mat$V1[i]+(readlength/8)), c(mat$V2[i], mat$V2[i]+1.5), col=as.character(mat$V4[i]), lwd=3)
		}
	}
}

#########################
# Pair connections 
bridge_x1<-bridgeraw$V1[(bridgeraw$V1[] >= start & bridgeraw$V3[] >= start) & (bridgeraw$V1[] <= stop & bridgeraw$V3[] <= stop)]
bridge_y1<-bridgeraw$V2[(bridgeraw$V1[] >= start & bridgeraw$V3[] >= start) & (bridgeraw$V1[] <= stop & bridgeraw$V3[] <= stop)]
bridge_x2<-bridgeraw$V3[(bridgeraw$V1[] >= start & bridgeraw$V3[] >= start) & (bridgeraw$V1[] <= stop & bridgeraw$V3[] <= stop)]
bridge_y2<-bridgeraw$V4[(bridgeraw$V1[] >= start & bridgeraw$V3[] >= start) & (bridgeraw$V1[] <= stop & bridgeraw$V3[] <= stop)]
bridge_col<-bridgeraw$V5[(bridgeraw$V1[] >= start & bridgeraw$V3[] >= start) & (bridgeraw$V1[] <= stop & bridgeraw$V3[] <= stop)]
segments(bridge_x1, bridge_y1, bridge_x2, bridge_y2, col=as.character(bridge_col))


########################
#
# Backbone
#
axis(1, at=c(start, seq(from=start+100, to=stop-100, by=100), stop), labels=c(start, seq(from=start+100, to=stop-100, by=100), stop), cex.axis=0.7)


dev.off()


