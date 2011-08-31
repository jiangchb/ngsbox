
pdf(file="CNV_NBS_LRR_active.pdf", width=12, height=6)

data<-read.table("shore_count.ALL_80.Bak-7.notcomplete.frac_back.txt")
nbs<-read.table("shore_count.ALL_80.Bak-7.notcomplete.frac_nbs.txt")


layoutmat=matrix(data=c(1,2), nrow=1, ncol=2)
layout(layoutmat)



hist(data$V1, xlim=c(-7,7), xlab="Col-0  <----- enrichment ---->  Accession", axes=F, freq=F, ylim=c(0,0.5), breaks=10, main="All genes in all accessions", col="orange")
axis(1, at=seq(-7,7), labels=c(seq(7,0), seq(1,7)))
axis(2)

hist(nbs$V1, xlim=c(-7,7), xlab="Col-0  <----- enrichment ---->  Accession", axes = F, freq=F, ylim=c(0, 0.5), main="NBS LRR loci in all accessions", col="red")
axis(1, at=seq(-7,7), labels=c(seq(7,0), seq(1,7)))
axis(2)






dev.off()

