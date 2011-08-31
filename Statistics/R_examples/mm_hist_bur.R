dups_bur<-read.table("avg_mm_of_duplications_bur")
cons_bur<-read.table("avg_mm_of_cons_simple_bur_col")

dups_tsu<-read.table("avg_mm_of_duplications_tsu")
cons_tsu<-read.table("avg_mm_of_cons_simple_tsu_col")

pdf(file="cov_corr_btw_samples_inlay.pdf", width=8, height=8)

layoutmatrix<-matrix(data=c(1,3,2,4), nrow=2, ncol=2)
layout(layoutmatrix)


## BUR

# Duplication
hist(dups_bur$V5, xlim=c(0,4), col="grey55", las=1, xlab="Average MM", main="Bur-0")
hist(dups_bur$V6, xlim=c(0,4), col="grey55", las=1, xlab="Average MM", main="Col-0", breaks=20)

# Uniform
hist(cons_bur$V5, xlim=c(0,4), col="grey75", las=1, xlab="Average MM", main="Bur-0", breaks=20)
hist(cons_bur$V6, xlim=c(0,4), col="grey75", las=1, xlab="Average MM", main="Col-0", breaks=20)



## TSU

# Duplication
hist(dups_tsu$V5, xlim=c(0,4), col="grey55", las=1, xlab="Average MM", main="Tsu-1", breaks=20)
hist(dups_tsu$V6, xlim=c(0,4), col="grey55", las=1, xlab="Average MM", main="Col-0", breaks=20)

# Uniform
hist(cons_tsu$V5, xlim=c(0,4), col="grey75", las=1, xlab="Average MM", main="Tsu-1", breaks=20)
hist(cons_tsu$V6, xlim=c(0,4), col="grey75", las=1, xlab="Average MM", main="Col-0", breaks=20)


dev.off()


