args <- commandArgs()
file<-args[5]
minsupp<-as.numeric(args[6]) # 1 = all
region<-as.numeric(args[7]) #xlim ylim
#if (minsupp < 2) {
#	minsupp<-2
#}
minsupp<-0

file2<-args[8]

##############################################

pdf(file="expected_coverage.pdf" , width=10, height=16)
layoutmax<-matrix(data=c(1,2), nrow=2, ncol=1)
layout(layoutmax)

####### Plot OBS versus EXP coverage ##########

data<-read.table(file)

obs<-data$V1[data$V3 >= minsupp]
exp<-data$V2[data$V3 >= minsupp]
supp<-data$V3[data$V3 >= minsupp]
sum<-data$V4[data$V3 >= minsupp]
ssum<-data$V5[data$V3 >= minsupp]

# supp * ssum can create overflow!!!
stddev<-sqrt(((supp * ssum) - sum^2) / (supp * (supp-1)))

plot(obs, exp, ylab="expected coverage", xlab="observed coverage", axes=T, ylim=c(0, region), xlim=c(0, region), type="n")
for (i in seq(10, 100, by=10)) {
	abline(v=i, lty=2, col="grey")
	abline(h=i, lty=2, col="grey")
}
for (i in 1:length(obs)) {
	segments(obs[i], exp[i]-stddev[i], obs[i], exp[i]+stddev[i], col="grey", lwd = 2)
}
lines(seq(1:region), seq(1:region), lty=2, col="darkgrey")

all<-sum(supp)
possupp<-(supp/all)*region
maxpossup<-max(possupp)
possupp<-possupp*(region/maxpossup)
lines(seq(1:region)-1, possupp[1:region], col="orange3", lwd=2)

points(obs, exp)

axis(4, at=c(0,region), labels=c(0, paste(round(maxpossup, digits=1), "%", sep="")))


########## Plot EXP versus OBS coverage ###########

data<-read.table(file2)

exp<-data$V1[data$V3 >= minsupp]
obs<-data$V2[data$V3 >= minsupp]
supp<-data$V3[data$V3 >= minsupp]
sum<-data$V4[data$V3 >= minsupp]
ssum<-data$V5[data$V3 >= minsupp]

stddev<-sqrt(((supp * ssum) - sum^2) / (supp * (supp-1)))

plot(exp, obs, xlab="expected coverage", ylab="observed coverage", axes=T, ylim=c(0, region), xlim=c(0, region), type="n")
for (i in seq(10, 100, by=10)) {
        abline(v=i, lty=2, col="grey")
        abline(h=i, lty=2, col="grey")
}
for (i in 1:length(obs)) {
        segments(obs[i], exp[i]-stddev[i], obs[i], exp[i]+stddev[i], col="grey", lwd = 2)
}
lines(seq(1:region), seq(1:region), lty=2, col="darkgrey")

all<-sum(supp)
possupp<-(supp/all)*region
maxpossup<-max(possupp)
possupp<-possupp*(region/maxpossup)
lines(seq(1:region)-1, possupp[1:region], col="orange3", lwd=2)

points(exp, obs)

axis(4, at=c(0,region), labels=c(0, paste(round(maxpossup, digits=1), "%", sep="")))


########### Thanks for flying with R ###############

dev.off()

