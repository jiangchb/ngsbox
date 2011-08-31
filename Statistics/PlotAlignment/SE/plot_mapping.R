# Written by Korbinian Schneeberger

args <- commandArgs()
chr<-as.numeric(args[12])
start<-as.numeric(args[6])
stop<-as.numeric(args[7])
format<-args[11]
mapping<-read.table(args[5])
seq<-read.table(args[8])
mm<-read.table(args[9], colClasses = c("numeric", "numeric", "character", "numeric"), header=T)
ins<-read.table(args[10], colClasses = c("numeric", "numeric", "character", "numeric"), header=T)

red<-"red1"
blue<-c("steelblue1", "steelblue4", "steelblue4", "steelblue4")
green<-"limegreen"

x<-c(start, stop)
maxy<-max(mapping$V2)
if (maxy < 10) {
	maxy<-10
}
y<-c(1, maxy)

if (format == "pdf") {
	pdf(file=paste("mapping_vis.", chr,"-", start ,"-",stop ,".pdf", sep=""), width=(stop-start+1)/7, height=maxy/3.7)
} else {
	png(filename=paste("mapping_vis.", start ,"-",stop ,".png", sep=""))
}
par(mai=c(0.75,0,0,0)) # c(bottom, left, top, right)'
plot(x, y, type="n", xlim=x, ylim=y, axes=F, ylab="", xlab="")

#########################
# Read Plotting
#
for (i in 1:length(mapping$V1)) {
	level<-mapping$V2[i]
	x1<-mapping$V1[i]
	x2<-mapping$V1[i]+mapping$V3[i]-1
	if (x2 >= start & x1 <= stop) { 
		if (x1 < start) {
			x1 <- start
		}
		if (x2 > stop) {
			x2 <- stop
		}
		segments(x1, level, x2, level, col=paste(blue[mapping$V4[i]]), lwd=2)
	}
}

########################
#
# Backbone
#
axis(1, at=seq(start,stop), labels=as.character(seq$V1), cex.axis=0.4)


#############################
# Mismatch Plotting
#
if (length(mm$i) > 0) {
	mm$prb<-mm$prb+6
	for (i in 1:length(mm$i)) {
		p<-mm$i[i]
		if (p >= start & p <= stop) {
			segments(p-0.5, mm$j[i], p+0.5, mm$j[i], col="white", lwd=2)
			abline(v=p, lwd=1, col="gray", lty=2)
			if (mm$prb[i] <= 11) { # mm has 6 added!
				text(p, mm$j[i], labels=mm$value[i], col="lightgrey", font=2)
			} 
			else {
				text(p, mm$j[i], labels=mm$value[i], col=heat.colors(65)[57-mm$prb[i]], font=2)
			}
		}
	}
}

#############################
# Insertion Plotting
#
ins$prb<-ins$prb+6
if (length(ins$i) > 0) {
	for (i in 1:length(ins$i)) {
        	p<-ins$i[i]
	        if (p+1 >= start & p-1 <= stop) {
			abline(v=p+0.5, lwd=1, col="gray", lty=2)
			points(p+0.5, ins$j[i]+0.27, cex=1, pch=25, col=blue, lwd=1)
	                text(p+0.5, ins$j[i]+0.5, labels=paste(ins$value[i], sep=""), col=heat.colors(65)[ins$V4[i]], font=2)
        	}
	}
}


dev.off()










