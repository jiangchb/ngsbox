
data<-read.table("shore_count.plot.txt")

pdf(file="CNVanno.pdf", width=20, height=10)

layoutmat<-matrix(data=c(1), nrow=1, ncol=1)
layout(layoutmat)


zoom=6



for (chr in 1:5) {

xval<-data$V2[data$V1==chr]
norm5<-(data$V5[data$V1==chr] + 1) / (sum(data$V5[data$V1==chr])/1000000)
norm6<-(data$V6[data$V1==chr] + 1) / (sum(data$V6[data$V1==chr])/1000000)
anno<-data$V8[data$V1==chr]


plot5=(norm5/norm6) - 1
plot6=((-1) * (norm6/norm5)) + 1

plot5[plot5>zoom]=zoom
plot6[plot6<(-zoom)]=-zoom

frame<-data.frame(xval, norm5, norm6, plot5, plot6, anno)


plot(frame$xval[frame$norm5>=frame$norm6], frame$plot5[frame$norm5>=frame$norm6], ylim=c(-zoom,zoom), axes=F, xlab=paste("Chr", chr, sep=""), ylab="Del-10              Col-0", xlim=c(0,31300000))

axis(1)
axis(2, las=2, at=(-zoom:(zoom)), labels=c(paste((zoom+1), "+", sep=""), zoom:2, (1:zoom), paste((zoom+1), "+", sep="")) )

points(frame$xval[frame$norm6>frame$norm5], frame$plot6[frame$norm6>frame$norm5])
points(frame$xval[frame$norm6>frame$norm5 & frame$anno=="NBS_LRR"], frame$plot6[frame$norm6>frame$norm5 & frame$anno=="NBS_LRR"], col="orange")


}


