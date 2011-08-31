data<-read.table("major_effects_by_gene_minusCol_bur.txt")

all<-data$V2
snps<-data$V3
del<-data$V4
ins<-data$V5
pr8_25<-data$V6
pr26_50<-data$V7
pr51_100<-data$V8
pr101<-data$V9


layoutmat<-matrix(data=c(1,2,3,4), nrow=2, ncol=2)
layout(layoutmat)

################# For pr of all length ###############

all6<-all
all6[all6[]>5]<-6


lambda<-mean(all)
x<-seq(0,max(all))
ypois<-c()
for (i in x) {
        ypois[i+1]<-(lambda^i * exp(-lambda) / factorial(i)) * length(all6)
}


summ<-0
for (i in 7:max(all)) {
	summ<-summ + ypois[i]
}
ypois6<-c(ypois[1], ypois[2], ypois[3], ypois[4], ypois[5], ypois[6], summ)


bins<-c(length(all6[all6[]==0]), length(all6[all6[]==1]), length(all6[all6[]==2]), length(all6[all6[]==3]), length(all6[all6[]==4]), length(all6[all6[]==5]), length(all6[all6[]==6]))

plot(0:6, bins, ylim=c(0,3000), type="l", xlim=c(0,6), col="red")
lines(0:6, ypois6)

plot(0:6, bins, type="l", xlim=c(0,6), col="red")
lines(0:6, ypois6)

chisq.test(bins, p=ypois6, rescale.p=T)


############### For PR of length > 100 ###########

all<-c()
all<-c(snps[]+del[]+ins[]+pr101[])


lambda<-mean(all)
x<-seq(0,max(all))
ypois<-c()
for (i in x) {
        ypois[i+1]<-(lambda^i * exp(-lambda) / factorial(i)) * length(all6)
}

all6<-all
all6[all6[]>5]<-6

summ<-0
for (i in 7:max(all)) {
        summ<-summ + ypois[i]
}
ypois6<-c(ypois[1], ypois[2], ypois[3], ypois[4], ypois[5], ypois[6], summ)


bins<-c(length(all6[all6[]==0]), length(all6[all6[]==1]), length(all6[all6[]==2]), length(all6[all6[]==3]), length(all6[all6[]==4]), length(all6[all6[]==5]), length(all6[all6[]==6]))

plot(0:6, bins, ylim=c(0,3000), type="l", xlim=c(0,6))
lines(0:6, ypois6)

plot(0:6, bins, type="l", xlim=c(0,6))
lines(0:6, ypois6)



