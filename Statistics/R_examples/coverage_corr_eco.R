###################################ä
####################################
### Contour plot

colbur<-read.table("cov_corr_col_bur.matrix")
coltsu<-read.table("cov_corr_col_tsu.matrix")
burtsu<-read.table("cov_corr_bur_tsu.matrix")

data<-read.table("cov_corr_eco.out")

data_colbur<-read.table("cov_corr_col_bur.uniq_simple")
data_coltsu<-read.table("cov_corr_col_tsu.uniq_simple")
data_burtsu<-read.table("cov_corr_bur_tsu.uniq_simple")


####################################
####################################
### Bur vs Col Contour plot

postscript(file="cov_corr_all_three_1.eps", width=10, height=10, horizontal=F, paper="special")


max<-300

x<-seq(1:max)
y<-seq(1:max)
z=matrix(data=colbur$V1, 300)

filled.contour(x=x, y=y, z=z, color=heat.colors, col=c("white", rev(heat.colors(27))), xlim=c(0,80), ylim=c(0,80), xlab="Col-0", ylab="Bur-0", plot.axes={axis(1); axis(2)})


dev.off();

#####################################
#####################################
### Bur vs Col dot plot

postscript(file="cov_corr_all_three_2.eps", width=10, height=10, horizontal=F, paper="special")


#plot(data$V3, data$V4, xlim=c(0,max), ylim=c(0,max), type="n", axes=F, xlab="Col-0", ylab="Bur-0")

max<-300

plot(data_colbur$V2, data_colbur$V1, xlim=c(0,max), ylim=c(0,max), type="n", axes=F, xlab="Col-0", ylab="Bur-0")
#polygon(c(0, 50, 50, 0), c(75, 75, 300, 300), col="grey55", border = NA)
#polygon(c(5, 50, 50, 5), c(5, 5, 50, 50), col="grey75", border = NA)
points(data_colbur$V2[data_colbur$V1[]<=300 & data_colbur$V2[]<=300], data_colbur$V1[data_colbur$V1[]<=300 & data_colbur$V2[]<=300], cex=0.25, col="steelblue4", pch=21, bg="steelblue4")
#rect(0, 75, 50, 300, border="grey55")
#rect(5, 5, 50, 50, border="grey75")
axis(1)
axis(2, las=1)

pearson<-cor.test(data$V3, data$V4, method="pearson")
spearman<-cor.test(data$V3, data$V4, method="spearman")

text(150, 50, labels=paste("Col-0 vs. Bur-0\n", " Pearson: ", round(pearson$estimate, digits=4), "  p:",pearson$p.value, "\n Spearman: ",round(spearman$estimate, digits=4), " p:", spearman$p.value, sep=""), adj=0)

#text(200, 50, labels=paste("Pearson: ", round(pearson$estimate, digits=4), "  p:",pearson$p.value, sep=""), adj=0)

dev.off()

######################################
######################################
### Tsu vs Col

postscript(file="cov_corr_all_three_3.eps", width=10, height=10, horizontal=F, paper="special")


max<-300

x<-seq(1:max)
y<-seq(1:max)
z=matrix(data=coltsu$V1, 300)

filled.contour(x=x, y=y, z=z, color=heat.colors, col=c("white", rev(heat.colors(18))), xlim=c(0,80), ylim=c(0,80), xlab="Col-0", ylab="Tsu-1", plot.axes={axis(1); axis(2)})

dev.off()
postscript(file="cov_corr_all_three_4.eps", width=10, height=10, horizontal=F, paper="special")

plot(data_coltsu$V2, data_coltsu$V1, xlim=c(0,max), ylim=c(0,max), type="n", axes=F, xlab="Col-0", ylab="Tsu-0")
#polygon(c(0, 50, 50, 0), c(75, 75, 300, 300), col="grey55", border = NA)
#polygon(c(5, 50, 50, 5), c(5, 5, 50, 50), col="grey75", border = NA)
points(data_coltsu$V2[data_coltsu$V1[]<=300 & data_coltsu$V2[]<=300], data_coltsu$V1[data_coltsu$V1[]<=300 & data_coltsu$V2[]<=300], cex=0.25, col="goldenrod1", pch=21, bg="goldenrod1")
#rect(0, 75, 50, 300, border="grey55")
#rect(5, 5, 50, 50, border="grey75")
axis(1)
axis(2, las=1)

pearson<-cor.test(data_coltsu$V1, data_coltsu$V2, method="pearson")
spearman<-cor.test(data_coltsu$V1, data_coltsu$V2, method="spearman")

#text(200, 50, labels=paste("Pearson: ", round(pearson$estimate, digits=4), "  p:",pearson$p.value, sep=""), adj=0)
text(150, 50, labels=paste("Col-0 vs. Tsu-1\n", " Pearson: ", round(pearson$estimate, digits=4), "  p:",pearson$p.value, "\n Spearman: ",round(spearman$estimate, digits=4), " p:", spearman$p.value, sep=""), adj=0)

dev.off()

######################################
######################################
### Tsu vs Bur

postscript(file="cov_corr_all_three_5.eps", width=10, height=10, horizontal=F, paper="special")


max<-300

x<-seq(1:max)
y<-seq(1:max)
z=matrix(data=burtsu$V1, 300)

filled.contour(x=x, y=y, z=z, color=heat.colors, col=c("white", rev(heat.colors(27))), xlim=c(0,80), ylim=c(0,80), ylab="Bur-0", xlab="Tsu-1", plot.axes={axis(1); axis(2)})


dev.off()

postscript(file="cov_corr_all_three_6.eps", width=10, height=10, horizontal=F, paper="special")


plot(data_burtsu$V1, data_burtsu$V2, xlim=c(0,max), ylim=c(0,max), type="n", axes=F, xlab="Bur-0", ylab="Tsu-0")
#polygon(c(0, 50, 50, 0), c(75, 75, 300, 300), col="grey55", border = NA)
#polygon(c(5, 50, 50, 5), c(5, 5, 50, 50), col="grey75", border = NA)
points(data_burtsu$V1[data_burtsu$V1[]<=300 & data_burtsu$V2[]<=300], data_burtsu$V2[data_burtsu$V1[]<=300 & data_burtsu$V2[]<=300], cex=0.25, col="darkolivegreen4", pch=21, bg="darkolivegreen4")
#rect(0, 75, 50, 300, border="grey55")
#rect(5, 5, 50, 50, border="grey75")
axis(1)
axis(2, las=1)

pearson<-cor.test(data_burtsu$V1, data_burtsu$V2, method="pearson")
spearman<-cor.test(data_burtsu$V1, data_burtsu$V2, method="spearman")

#text(200, 50, labels=paste("Pearson: ", round(pearson$estimate, digits=4), "  p:",pearson$p.value, sep=""), adj=0)
text(150, 50, labels=paste("Tsu-1 vs. Bur-0\n", " Pearson: ", round(pearson$estimate, digits=4), "  p:",pearson$p.value, "\n Spearman: ",round(spearman$estimate, digits=4), " p:", spearman$p.value, sep=""), adj=0)

dev.off()



