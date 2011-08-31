pdf(file=("reference_error_analysis.pdf"), width=16, height=14)

prephasingrefread<-read.table("error_prephasing_ref_vs_read.txt");
prephasingreadread<-read.table("error_prephasing_read_vs_read.txt");
prephasingrefref<-read.table("error_prephasing_ref_vs_ref.txt");
prephasingreadref<-read.table("error_prephasing_read_vs_ref.txt");
phasingrefread<-read.table("error_phasing_ref_vs_read.txt");
phasingreadread<-read.table("error_phasing_read_vs_read.txt");
phasingrefref<-read.table("error_phasing_ref_vs_ref.txt");
phasingreadref<-read.table("error_phasing_read_vs_ref.txt");
basecomp<-read.table("error_basecomparison.txt");
errornumread<-read.table("error_num_per_read.txt");
errorpos<-read.table("error_positionwise.txt");
dist2<-read.table("error_distances_2.txt");
dist3<-read.table("error_distances_3.txt");
dist4<-read.table("error_distances_4.txt");
distmax<-read.table("error_distances_max.txt");
dist2rand<-read.table("error_distances_2_rand.txt");
dist3rand<-read.table("error_distances_3_rand.txt");
dist4rand<-read.table("error_distances_4_rand.txt");
distmaxrand<-read.table("error_distances_max_rand.txt");
fstocc<-read.table("error_first_occurrence.txt")

########### Prephasing detection #########################################################

layoutmat<-matrix(data=c(1,2,3,4), nrow=2, ncol=2)
layout(layoutmat)

# Plot Matrix 
plot(c(0, 160), c(0, 160), type="n", xlab="", ylab="", axes=F, main="Prephasing - MM and its successor\nY-axes: MM   X-axes: successor")
for (i in 2:6) {
# horizontal
lines(c(10, 120), c(i*20,i*20), col="grey")
}
for (i in 1:4) {
# vertical
lines(c(i*20, i*20), c(30,130), col="grey") 
}
# horizontal
text(10, 110, labels="A", col="black")
text(10, 90, labels="C", col="black")
text(10, 70, labels="G", col="black")
text(10, 50, labels="T", col="black")
# vertical
text(30, 130, labels="A", col="black")
text(50, 130, labels="C", col="black")
text(70, 130, labels="G", col="black")
text(90, 130, labels="T", col="black")

text(0, 80, labels="Error (Read)",  col="black", srt=90 )
text(70, 150, labels="Successor (Ref)", col="black")

all_a<-sum(prephasingrefread$V3[prephasingrefread$V1=="A"])
all_c<-sum(prephasingrefread$V3[prephasingrefread$V1=="C"])
all_g<-sum(prephasingrefread$V3[prephasingrefread$V1=="G"])
all_t<-sum(prephasingrefread$V3[prephasingrefread$V1=="T"])

text(30, 110, labels=round(100*prephasingrefread$V3[prephasingrefread$V2=="A" & prephasingrefread$V1=="A"]/all_a, digits=1), col="deepskyblue4")
text(50, 110, labels=round(100*prephasingrefread$V3[prephasingrefread$V2=="C" & prephasingrefread$V1=="A"]/all_a, digits=1), col="deepskyblue4")
text(70, 110, labels=round(100*prephasingrefread$V3[prephasingrefread$V2=="G" & prephasingrefread$V1=="A"]/all_a, digits=1), col="deepskyblue4")
text(90, 110, labels=round(100*prephasingrefread$V3[prephasingrefread$V2=="T" & prephasingrefread$V1=="A"]/all_a, digits=1), col="deepskyblue4")

text(30, 90, labels=round(100*prephasingrefread$V3[prephasingrefread$V2=="A" & prephasingrefread$V1=="C"]/all_c, digits=1), col="deepskyblue4")
text(50, 90, labels=round(100*prephasingrefread$V3[prephasingrefread$V2=="C" & prephasingrefread$V1=="C"]/all_c, digits=1), col="deepskyblue4")
text(70, 90, labels=round(100*prephasingrefread$V3[prephasingrefread$V2=="G" & prephasingrefread$V1=="C"]/all_c, digits=1), col="deepskyblue4")
text(90, 90, labels=round(100*prephasingrefread$V3[prephasingrefread$V2=="T" & prephasingrefread$V1=="C"]/all_c, digits=1), col="deepskyblue4")

text(30, 70, labels=round(100*prephasingrefread$V3[prephasingrefread$V2=="A" & prephasingrefread$V1=="G"]/all_g, digits=1), col="deepskyblue4")
text(50, 70, labels=round(100*prephasingrefread$V3[prephasingrefread$V2=="C" & prephasingrefread$V1=="G"]/all_g, digits=1), col="deepskyblue4")
text(70, 70, labels=round(100*prephasingrefread$V3[prephasingrefread$V2=="G" & prephasingrefread$V1=="G"]/all_g, digits=1), col="deepskyblue4")
text(90, 70, labels=round(100*prephasingrefread$V3[prephasingrefread$V2=="T" & prephasingrefread$V1=="G"]/all_g, digits=1), col="deepskyblue4")

text(30, 50, labels=round(100*prephasingrefread$V3[prephasingrefread$V2=="A" & prephasingrefread$V1=="T"]/all_t, digits=1), col="deepskyblue4")
text(50, 50, labels=round(100*prephasingrefread$V3[prephasingrefread$V2=="C" & prephasingrefread$V1=="T"]/all_t, digits=1), col="deepskyblue4")
text(70, 50, labels=round(100*prephasingrefread$V3[prephasingrefread$V2=="G" & prephasingrefread$V1=="T"]/all_t, digits=1), col="deepskyblue4")
text(90, 50, labels=round(100*prephasingrefread$V3[prephasingrefread$V2=="T" & prephasingrefread$V1=="T"]/all_t, digits=1), col="deepskyblue4")


# Plot Matrix 
plot(c(0, 160), c(0, 160), type="n", xlab="", ylab="", axes=F, main="Phasing - MM and its successor\nY-axes: MM   X-axes: successor")
for (i in 2:6) {
# horizontal
lines(c(10, 120), c(i*20,i*20), col="grey")
}
for (i in 1:4) {
# vertical
lines(c(i*20, i*20), c(30,130), col="grey") 
}
# horizontal
text(10, 110, labels="A", col="black")
text(10, 90, labels="C", col="black")
text(10, 70, labels="G", col="black")
text(10, 50, labels="T", col="black")
# vertical
text(30, 130, labels="A", col="black")
text(50, 130, labels="C", col="black")
text(70, 130, labels="G", col="black")
text(90, 130, labels="T", col="black")

text(0, 80, labels="Error (Read)",  col="black", srt=90 )
text(70, 150, labels="Successor (Read)", col="black")

all_a<-sum(prephasingreadread$V3[prephasingreadread$V1=="A"])
all_c<-sum(prephasingreadread$V3[prephasingreadread$V1=="C"])
all_g<-sum(prephasingreadread$V3[prephasingreadread$V1=="G"])
all_t<-sum(prephasingreadread$V3[prephasingreadread$V1=="T"])

text(30, 110, labels=round(100*prephasingreadread$V3[prephasingreadread$V2=="A" & prephasingreadread$V1=="A"]/all_a, digits=1), col="deepskyblue4")
text(50, 110, labels=round(100*prephasingreadread$V3[prephasingreadread$V2=="C" & prephasingreadread$V1=="A"]/all_a, digits=1), col="deepskyblue4")
text(70, 110, labels=round(100*prephasingreadread$V3[prephasingreadread$V2=="G" & prephasingreadread$V1=="A"]/all_a, digits=1), col="deepskyblue4")
text(90, 110, labels=round(100*prephasingreadread$V3[prephasingreadread$V2=="T" & prephasingreadread$V1=="A"]/all_a, digits=1), col="deepskyblue4")

text(30, 90, labels=round(100*prephasingreadread$V3[prephasingreadread$V2=="A" & prephasingreadread$V1=="C"]/all_c, digits=1), col="deepskyblue4")
text(50, 90, labels=round(100*prephasingreadread$V3[prephasingreadread$V2=="C" & prephasingreadread$V1=="C"]/all_c, digits=1), col="deepskyblue4")
text(70, 90, labels=round(100*prephasingreadread$V3[prephasingreadread$V2=="G" & prephasingreadread$V1=="C"]/all_c, digits=1), col="deepskyblue4")
text(90, 90, labels=round(100*prephasingreadread$V3[prephasingreadread$V2=="T" & prephasingreadread$V1=="C"]/all_c, digits=1), col="deepskyblue4")

text(30, 70, labels=round(100*prephasingreadread$V3[prephasingreadread$V2=="A" & prephasingreadread$V1=="G"]/all_g, digits=1), col="deepskyblue4")
text(50, 70, labels=round(100*prephasingreadread$V3[prephasingreadread$V2=="C" & prephasingreadread$V1=="G"]/all_g, digits=1), col="deepskyblue4")
text(70, 70, labels=round(100*prephasingreadread$V3[prephasingreadread$V2=="G" & prephasingreadread$V1=="G"]/all_g, digits=1), col="deepskyblue4")
text(90, 70, labels=round(100*prephasingreadread$V3[prephasingreadread$V2=="T" & prephasingreadread$V1=="G"]/all_g, digits=1), col="deepskyblue4")

text(30, 50, labels=round(100*prephasingreadread$V3[prephasingreadread$V2=="A" & prephasingreadread$V1=="T"]/all_t, digits=1), col="deepskyblue4")
text(50, 50, labels=round(100*prephasingreadread$V3[prephasingreadread$V2=="C" & prephasingreadread$V1=="T"]/all_t, digits=1), col="deepskyblue4")
text(70, 50, labels=round(100*prephasingreadread$V3[prephasingreadread$V2=="G" & prephasingreadread$V1=="T"]/all_t, digits=1), col="deepskyblue4")
text(90, 50, labels=round(100*prephasingreadread$V3[prephasingreadread$V2=="T" & prephasingreadread$V1=="T"]/all_t, digits=1), col="deepskyblue4")

# Plot Matrix 
plot(c(0, 160), c(0, 160), type="n", xlab="", ylab="", axes=F, main="Phasing - MM and its successor\nY-axes: MM   X-axes: successor")
for (i in 2:6) {
# horizontal
lines(c(10, 120), c(i*20,i*20), col="grey")
}
for (i in 1:4) {
# vertical
lines(c(i*20, i*20), c(30,130), col="grey") 
}
# horizontal
text(10, 110, labels="A", col="black")
text(10, 90, labels="C", col="black")
text(10, 70, labels="G", col="black")
text(10, 50, labels="T", col="black")
# vertical
text(30, 130, labels="A", col="black")
text(50, 130, labels="C", col="black")
text(70, 130, labels="G", col="black")
text(90, 130, labels="T", col="black")

text(0, 80, labels="Error (Ref)",  col="black", srt=90 )
text(70, 150, labels="Successor (Read)", col="black")

all_a<-sum(prephasingreadref$V3[prephasingreadref$V1=="A"])
all_c<-sum(prephasingreadref$V3[prephasingreadref$V1=="C"])
all_g<-sum(prephasingreadref$V3[prephasingreadref$V1=="G"])
all_t<-sum(prephasingreadref$V3[prephasingreadref$V1=="T"])

text(30, 110, labels=round(100*prephasingreadref$V3[prephasingreadref$V2=="A" & prephasingreadref$V1=="A"]/all_a, digits=1), col="deepskyblue4")
text(50, 110, labels=round(100*prephasingreadref$V3[prephasingreadref$V2=="C" & prephasingreadref$V1=="A"]/all_a, digits=1), col="deepskyblue4")
text(70, 110, labels=round(100*prephasingreadref$V3[prephasingreadref$V2=="G" & prephasingreadref$V1=="A"]/all_a, digits=1), col="deepskyblue4")
text(90, 110, labels=round(100*prephasingreadref$V3[prephasingreadref$V2=="T" & prephasingreadref$V1=="A"]/all_a, digits=1), col="deepskyblue4")

text(30, 90, labels=round(100*prephasingreadref$V3[prephasingreadref$V2=="A" & prephasingreadref$V1=="C"]/all_c, digits=1), col="deepskyblue4")
text(50, 90, labels=round(100*prephasingreadref$V3[prephasingreadref$V2=="C" & prephasingreadref$V1=="C"]/all_c, digits=1), col="deepskyblue4")
text(70, 90, labels=round(100*prephasingreadref$V3[prephasingreadref$V2=="G" & prephasingreadref$V1=="C"]/all_c, digits=1), col="deepskyblue4")
text(90, 90, labels=round(100*prephasingreadref$V3[prephasingreadref$V2=="T" & prephasingreadref$V1=="C"]/all_c, digits=1), col="deepskyblue4")

text(30, 70, labels=round(100*prephasingreadref$V3[prephasingreadref$V2=="A" & prephasingreadref$V1=="G"]/all_g, digits=1), col="deepskyblue4")
text(50, 70, labels=round(100*prephasingreadref$V3[prephasingreadref$V2=="C" & prephasingreadref$V1=="G"]/all_g, digits=1), col="deepskyblue4")
text(70, 70, labels=round(100*prephasingreadref$V3[prephasingreadref$V2=="G" & prephasingreadref$V1=="G"]/all_g, digits=1), col="deepskyblue4")
text(90, 70, labels=round(100*prephasingreadref$V3[prephasingreadref$V2=="T" & prephasingreadref$V1=="G"]/all_g, digits=1), col="deepskyblue4")

text(30, 50, labels=round(100*prephasingreadref$V3[prephasingreadref$V2=="A" & prephasingreadref$V1=="T"]/all_t, digits=1), col="deepskyblue4")
text(50, 50, labels=round(100*prephasingreadref$V3[prephasingreadref$V2=="C" & prephasingreadref$V1=="T"]/all_t, digits=1), col="deepskyblue4")
text(70, 50, labels=round(100*prephasingreadref$V3[prephasingreadref$V2=="G" & prephasingreadref$V1=="T"]/all_t, digits=1), col="deepskyblue4")
text(90, 50, labels=round(100*prephasingreadref$V3[prephasingreadref$V2=="T" & prephasingreadref$V1=="T"]/all_t, digits=1), col="deepskyblue4")


# Plot Matrix 
plot(c(0, 160), c(0, 160), type="n", xlab="", ylab="", axes=F, main="Phasing - MM and its successor\nY-axes: MM   X-axes: successor")
for (i in 2:6) {
# horizontal
lines(c(10, 120), c(i*20,i*20), col="grey")
}
for (i in 1:4) {
# vertical
lines(c(i*20, i*20), c(30,130), col="grey") 
}
# horizontal
text(10, 110, labels="A", col="black")
text(10, 90, labels="C", col="black")
text(10, 70, labels="G", col="black")
text(10, 50, labels="T", col="black")
# vertical
text(30, 130, labels="A", col="black")
text(50, 130, labels="C", col="black")
text(70, 130, labels="G", col="black")
text(90, 130, labels="T", col="black")

text(0, 80, labels="Error (Ref)",  col="black", srt=90 )
text(70, 150, labels="Successor (Ref)", col="black")

all_a<-sum(prephasingrefref$V3[prephasingrefref$V1=="A"])
all_c<-sum(prephasingrefref$V3[prephasingrefref$V1=="C"])
all_g<-sum(prephasingrefref$V3[prephasingrefref$V1=="G"])
all_t<-sum(prephasingrefref$V3[prephasingrefref$V1=="T"])

text(30, 110, labels=round(100*prephasingrefref$V3[prephasingrefref$V2=="A" & prephasingrefref$V1=="A"]/all_a, digits=1), col="deepskyblue4")
text(50, 110, labels=round(100*prephasingrefref$V3[prephasingrefref$V2=="C" & prephasingrefref$V1=="A"]/all_a, digits=1), col="deepskyblue4")
text(70, 110, labels=round(100*prephasingrefref$V3[prephasingrefref$V2=="G" & prephasingrefref$V1=="A"]/all_a, digits=1), col="deepskyblue4")
text(90, 110, labels=round(100*prephasingrefref$V3[prephasingrefref$V2=="T" & prephasingrefref$V1=="A"]/all_a, digits=1), col="deepskyblue4")

text(30, 90, labels=round(100*prephasingrefref$V3[prephasingrefref$V2=="A" & prephasingrefref$V1=="C"]/all_c, digits=1), col="deepskyblue4")
text(50, 90, labels=round(100*prephasingrefref$V3[prephasingrefref$V2=="C" & prephasingrefref$V1=="C"]/all_c, digits=1), col="deepskyblue4")
text(70, 90, labels=round(100*prephasingrefref$V3[prephasingrefref$V2=="G" & prephasingrefref$V1=="C"]/all_c, digits=1), col="deepskyblue4")
text(90, 90, labels=round(100*prephasingrefref$V3[prephasingrefref$V2=="T" & prephasingrefref$V1=="C"]/all_c, digits=1), col="deepskyblue4")

text(30, 70, labels=round(100*prephasingrefref$V3[prephasingrefref$V2=="A" & prephasingrefref$V1=="G"]/all_g, digits=1), col="deepskyblue4")
text(50, 70, labels=round(100*prephasingrefref$V3[prephasingrefref$V2=="C" & prephasingrefref$V1=="G"]/all_g, digits=1), col="deepskyblue4")
text(70, 70, labels=round(100*prephasingrefref$V3[prephasingrefref$V2=="G" & prephasingrefref$V1=="G"]/all_g, digits=1), col="deepskyblue4")
text(90, 70, labels=round(100*prephasingrefref$V3[prephasingrefref$V2=="T" & prephasingrefref$V1=="G"]/all_g, digits=1), col="deepskyblue4")

text(30, 50, labels=round(100*prephasingrefref$V3[prephasingrefref$V2=="A" & prephasingrefref$V1=="T"]/all_t, digits=1), col="deepskyblue4")
text(50, 50, labels=round(100*prephasingrefref$V3[prephasingrefref$V2=="C" & prephasingrefref$V1=="T"]/all_t, digits=1), col="deepskyblue4")
text(70, 50, labels=round(100*prephasingrefref$V3[prephasingrefref$V2=="G" & prephasingrefref$V1=="T"]/all_t, digits=1), col="deepskyblue4")
text(90, 50, labels=round(100*prephasingrefref$V3[prephasingrefref$V2=="T" & prephasingrefref$V1=="T"]/all_t, digits=1), col="deepskyblue4")



########### Phasing detection #########################################################

layoutmat<-matrix(data=c(1,2,3,4), nrow=2, ncol=2)
layout(layoutmat)

# Plot Matrix 
plot(c(0, 160), c(0, 160), type="n", xlab="", ylab="", axes=F, main="Phasing - MM and its predecessor\nY-axes: MM   X-axes: predecessor")
for (i in 2:6) {
# horizontal
lines(c(10, 120), c(i*20,i*20), col="grey")
}
for (i in 1:4) {
# vertical
lines(c(i*20, i*20), c(30,130), col="grey") 
}
# horizontal
text(10, 110, labels="A", col="black")
text(10, 90, labels="C", col="black")
text(10, 70, labels="G", col="black")
text(10, 50, labels="T", col="black")
# vertical
text(30, 130, labels="A", col="black")
text(50, 130, labels="C", col="black")
text(70, 130, labels="G", col="black")
text(90, 130, labels="T", col="black")

text(0, 80, labels="Error (Read)",  col="black", srt=90 )
text(70, 150, labels="Predecessor (Ref)", col="black")

all_a<-sum(phasingrefread$V3[phasingrefread$V2=="A"])
all_c<-sum(phasingrefread$V3[phasingrefread$V2=="C"])
all_g<-sum(phasingrefread$V3[phasingrefread$V2=="G"])
all_t<-sum(phasingrefread$V3[phasingrefread$V2=="T"])

text(30, 110, labels=round(100*phasingrefread$V3[phasingrefread$V1=="A" & phasingrefread$V2=="A"]/all_a, digits=1), col="deepskyblue4")
text(50, 110, labels=round(100*phasingrefread$V3[phasingrefread$V1=="C" & phasingrefread$V2=="A"]/all_a, digits=1), col="deepskyblue4")
text(70, 110, labels=round(100*phasingrefread$V3[phasingrefread$V1=="G" & phasingrefread$V2=="A"]/all_a, digits=1), col="deepskyblue4")
text(90, 110, labels=round(100*phasingrefread$V3[phasingrefread$V1=="T" & phasingrefread$V2=="A"]/all_a, digits=1), col="deepskyblue4")

text(30, 90, labels=round(100*phasingrefread$V3[phasingrefread$V1=="A" & phasingrefread$V2=="C"]/all_c, digits=1), col="deepskyblue4")
text(50, 90, labels=round(100*phasingrefread$V3[phasingrefread$V1=="C" & phasingrefread$V2=="C"]/all_c, digits=1), col="deepskyblue4")
text(70, 90, labels=round(100*phasingrefread$V3[phasingrefread$V1=="G" & phasingrefread$V2=="C"]/all_c, digits=1), col="deepskyblue4")
text(90, 90, labels=round(100*phasingrefread$V3[phasingrefread$V1=="T" & phasingrefread$V2=="C"]/all_c, digits=1), col="deepskyblue4")

text(30, 70, labels=round(100*phasingrefread$V3[phasingrefread$V1=="A" & phasingrefread$V2=="G"]/all_g, digits=1), col="deepskyblue4")
text(50, 70, labels=round(100*phasingrefread$V3[phasingrefread$V1=="C" & phasingrefread$V2=="G"]/all_g, digits=1), col="deepskyblue4")
text(70, 70, labels=round(100*phasingrefread$V3[phasingrefread$V1=="G" & phasingrefread$V2=="G"]/all_g, digits=1), col="deepskyblue4")
text(90, 70, labels=round(100*phasingrefread$V3[phasingrefread$V1=="T" & phasingrefread$V2=="G"]/all_g, digits=1), col="deepskyblue4")

text(30, 50, labels=round(100*phasingrefread$V3[phasingrefread$V1=="A" & phasingrefread$V2=="T"]/all_t, digits=1), col="deepskyblue4")
text(50, 50, labels=round(100*phasingrefread$V3[phasingrefread$V1=="C" & phasingrefread$V2=="T"]/all_t, digits=1), col="deepskyblue4")
text(70, 50, labels=round(100*phasingrefread$V3[phasingrefread$V1=="G" & phasingrefread$V2=="T"]/all_t, digits=1), col="deepskyblue4")
text(90, 50, labels=round(100*phasingrefread$V3[phasingrefread$V1=="T" & phasingrefread$V2=="T"]/all_t, digits=1), col="deepskyblue4")


# Plot Matrix 
plot(c(0, 160), c(0, 160), type="n", xlab="", ylab="", axes=F, main="Phasing - MM and its predecessor\nY-axes: MM   X-axes: predecessor")
for (i in 2:6) {
# horizontal
lines(c(10, 120), c(i*20,i*20), col="grey")
}
for (i in 1:4) {
# vertical
lines(c(i*20, i*20), c(30,130), col="grey") 
}
# horizontal
text(10, 110, labels="A", col="black")
text(10, 90, labels="C", col="black")
text(10, 70, labels="G", col="black")
text(10, 50, labels="T", col="black")
# vertical
text(30, 130, labels="A", col="black")
text(50, 130, labels="C", col="black")
text(70, 130, labels="G", col="black")
text(90, 130, labels="T", col="black")

text(0, 80, labels="Error (Read)",  col="black", srt=90 )
text(70, 150, labels="Predecessor (Read)", col="black")

all_a<-sum(phasingreadread$V3[phasingreadread$V2=="A"])
all_c<-sum(phasingreadread$V3[phasingreadread$V2=="C"])
all_g<-sum(phasingreadread$V3[phasingreadread$V2=="G"])
all_t<-sum(phasingreadread$V3[phasingreadread$V2=="T"])

text(30, 110, labels=round(100*phasingreadread$V3[phasingreadread$V1=="A" & phasingreadread$V2=="A"]/all_a, digits=1), col="deepskyblue4")
text(50, 110, labels=round(100*phasingreadread$V3[phasingreadread$V1=="C" & phasingreadread$V2=="A"]/all_a, digits=1), col="deepskyblue4")
text(70, 110, labels=round(100*phasingreadread$V3[phasingreadread$V1=="G" & phasingreadread$V2=="A"]/all_a, digits=1), col="deepskyblue4")
text(90, 110, labels=round(100*phasingreadread$V3[phasingreadread$V1=="T" & phasingreadread$V2=="A"]/all_a, digits=1), col="deepskyblue4")

text(30, 90, labels=round(100*phasingreadread$V3[phasingreadread$V1=="A" & phasingreadread$V2=="C"]/all_c, digits=1), col="deepskyblue4")
text(50, 90, labels=round(100*phasingreadread$V3[phasingreadread$V1=="C" & phasingreadread$V2=="C"]/all_c, digits=1), col="deepskyblue4")
text(70, 90, labels=round(100*phasingreadread$V3[phasingreadread$V1=="G" & phasingreadread$V2=="C"]/all_c, digits=1), col="deepskyblue4")
text(90, 90, labels=round(100*phasingreadread$V3[phasingreadread$V1=="T" & phasingreadread$V2=="C"]/all_c, digits=1), col="deepskyblue4")

text(30, 70, labels=round(100*phasingreadread$V3[phasingreadread$V1=="A" & phasingreadread$V2=="G"]/all_g, digits=1), col="deepskyblue4")
text(50, 70, labels=round(100*phasingreadread$V3[phasingreadread$V1=="C" & phasingreadread$V2=="G"]/all_g, digits=1), col="deepskyblue4")
text(70, 70, labels=round(100*phasingreadread$V3[phasingreadread$V1=="G" & phasingreadread$V2=="G"]/all_g, digits=1), col="deepskyblue4")
text(90, 70, labels=round(100*phasingreadread$V3[phasingreadread$V1=="T" & phasingreadread$V2=="G"]/all_g, digits=1), col="deepskyblue4")

text(30, 50, labels=round(100*phasingreadread$V3[phasingreadread$V1=="A" & phasingreadread$V2=="T"]/all_t, digits=1), col="deepskyblue4")
text(50, 50, labels=round(100*phasingreadread$V3[phasingreadread$V1=="C" & phasingreadread$V2=="T"]/all_t, digits=1), col="deepskyblue4")
text(70, 50, labels=round(100*phasingreadread$V3[phasingreadread$V1=="G" & phasingreadread$V2=="T"]/all_t, digits=1), col="deepskyblue4")
text(90, 50, labels=round(100*phasingreadread$V3[phasingreadread$V1=="T" & phasingreadread$V2=="T"]/all_t, digits=1), col="deepskyblue4")

# Plot Matrix 
plot(c(0, 160), c(0, 160), type="n", xlab="", ylab="", axes=F, main="Phasing - MM and its predecessor\nY-axes: MM   X-axes: predecessor")
for (i in 2:6) {
# horizontal
lines(c(10, 120), c(i*20,i*20), col="grey")
}
for (i in 1:4) {
# vertical
lines(c(i*20, i*20), c(30,130), col="grey") 
}
# horizontal
text(10, 110, labels="A", col="black")
text(10, 90, labels="C", col="black")
text(10, 70, labels="G", col="black")
text(10, 50, labels="T", col="black")
# vertical
text(30, 130, labels="A", col="black")
text(50, 130, labels="C", col="black")
text(70, 130, labels="G", col="black")
text(90, 130, labels="T", col="black")

text(0, 80, labels="Error (Ref)",  col="black", srt=90 )
text(70, 150, labels="Predecessor (Read)", col="black")

all_a<-sum(phasingreadref$V3[phasingreadref$V2=="A"])
all_c<-sum(phasingreadref$V3[phasingreadref$V2=="C"])
all_g<-sum(phasingreadref$V3[phasingreadref$V2=="G"])
all_t<-sum(phasingreadref$V3[phasingreadref$V2=="T"])

text(30, 110, labels=round(100*phasingreadref$V3[phasingreadref$V1=="A" & phasingreadref$V2=="A"]/all_a, digits=1), col="deepskyblue4")
text(50, 110, labels=round(100*phasingreadref$V3[phasingreadref$V1=="C" & phasingreadref$V2=="A"]/all_a, digits=1), col="deepskyblue4")
text(70, 110, labels=round(100*phasingreadref$V3[phasingreadref$V1=="G" & phasingreadref$V2=="A"]/all_a, digits=1), col="deepskyblue4")
text(90, 110, labels=round(100*phasingreadref$V3[phasingreadref$V1=="T" & phasingreadref$V2=="A"]/all_a, digits=1), col="deepskyblue4")

text(30, 90, labels=round(100*phasingreadref$V3[phasingreadref$V1=="A" & phasingreadref$V2=="C"]/all_c, digits=1), col="deepskyblue4")
text(50, 90, labels=round(100*phasingreadref$V3[phasingreadref$V1=="C" & phasingreadref$V2=="C"]/all_c, digits=1), col="deepskyblue4")
text(70, 90, labels=round(100*phasingreadref$V3[phasingreadref$V1=="G" & phasingreadref$V2=="C"]/all_c, digits=1), col="deepskyblue4")
text(90, 90, labels=round(100*phasingreadref$V3[phasingreadref$V1=="T" & phasingreadref$V2=="C"]/all_c, digits=1), col="deepskyblue4")

text(30, 70, labels=round(100*phasingreadref$V3[phasingreadref$V1=="A" & phasingreadref$V2=="G"]/all_g, digits=1), col="deepskyblue4")
text(50, 70, labels=round(100*phasingreadref$V3[phasingreadref$V1=="C" & phasingreadref$V2=="G"]/all_g, digits=1), col="deepskyblue4")
text(70, 70, labels=round(100*phasingreadref$V3[phasingreadref$V1=="G" & phasingreadref$V2=="G"]/all_g, digits=1), col="deepskyblue4")
text(90, 70, labels=round(100*phasingreadref$V3[phasingreadref$V1=="T" & phasingreadref$V2=="G"]/all_g, digits=1), col="deepskyblue4")

text(30, 50, labels=round(100*phasingreadref$V3[phasingreadref$V1=="A" & phasingreadref$V2=="T"]/all_t, digits=1), col="deepskyblue4")
text(50, 50, labels=round(100*phasingreadref$V3[phasingreadref$V1=="C" & phasingreadref$V2=="T"]/all_t, digits=1), col="deepskyblue4")
text(70, 50, labels=round(100*phasingreadref$V3[phasingreadref$V1=="G" & phasingreadref$V2=="T"]/all_t, digits=1), col="deepskyblue4")
text(90, 50, labels=round(100*phasingreadref$V3[phasingreadref$V1=="T" & phasingreadref$V2=="T"]/all_t, digits=1), col="deepskyblue4")


# Plot Matrix 
plot(c(0, 160), c(0, 160), type="n", xlab="", ylab="", axes=F, main="Phasing - MM and its predecessor\nY-axes: MM   X-axes: predecessor")
for (i in 2:6) {
# horizontal
lines(c(10, 120), c(i*20,i*20), col="grey")
}
for (i in 1:4) {
# vertical
lines(c(i*20, i*20), c(30,130), col="grey") 
}
# horizontal
text(10, 110, labels="A", col="black")
text(10, 90, labels="C", col="black")
text(10, 70, labels="G", col="black")
text(10, 50, labels="T", col="black")
# vertical
text(30, 130, labels="A", col="black")
text(50, 130, labels="C", col="black")
text(70, 130, labels="G", col="black")
text(90, 130, labels="T", col="black")

text(0, 80, labels="Error (Ref)",  col="black", srt=90 )
text(70, 150, labels="Predecessor (Ref)", col="black")

all_a<-sum(phasingrefref$V3[phasingrefref$V2=="A"])
all_c<-sum(phasingrefref$V3[phasingrefref$V2=="C"])
all_g<-sum(phasingrefref$V3[phasingrefref$V2=="G"])
all_t<-sum(phasingrefref$V3[phasingrefref$V2=="T"])

text(30, 110, labels=round(100*phasingrefref$V3[phasingrefref$V1=="A" & phasingrefref$V2=="A"]/all_a, digits=1), col="deepskyblue4")
text(50, 110, labels=round(100*phasingrefref$V3[phasingrefref$V1=="C" & phasingrefref$V2=="A"]/all_a, digits=1), col="deepskyblue4")
text(70, 110, labels=round(100*phasingrefref$V3[phasingrefref$V1=="G" & phasingrefref$V2=="A"]/all_a, digits=1), col="deepskyblue4")
text(90, 110, labels=round(100*phasingrefref$V3[phasingrefref$V1=="T" & phasingrefref$V2=="A"]/all_a, digits=1), col="deepskyblue4")

text(30, 90, labels=round(100*phasingrefref$V3[phasingrefref$V1=="A" & phasingrefref$V2=="C"]/all_c, digits=1), col="deepskyblue4")
text(50, 90, labels=round(100*phasingrefref$V3[phasingrefref$V1=="C" & phasingrefref$V2=="C"]/all_c, digits=1), col="deepskyblue4")
text(70, 90, labels=round(100*phasingrefref$V3[phasingrefref$V1=="G" & phasingrefref$V2=="C"]/all_c, digits=1), col="deepskyblue4")
text(90, 90, labels=round(100*phasingrefref$V3[phasingrefref$V1=="T" & phasingrefref$V2=="C"]/all_c, digits=1), col="deepskyblue4")

text(30, 70, labels=round(100*phasingrefref$V3[phasingrefref$V1=="A" & phasingrefref$V2=="G"]/all_g, digits=1), col="deepskyblue4")
text(50, 70, labels=round(100*phasingrefref$V3[phasingrefref$V1=="C" & phasingrefref$V2=="G"]/all_g, digits=1), col="deepskyblue4")
text(70, 70, labels=round(100*phasingrefref$V3[phasingrefref$V1=="G" & phasingrefref$V2=="G"]/all_g, digits=1), col="deepskyblue4")
text(90, 70, labels=round(100*phasingrefref$V3[phasingrefref$V1=="T" & phasingrefref$V2=="G"]/all_g, digits=1), col="deepskyblue4")

text(30, 50, labels=round(100*phasingrefref$V3[phasingrefref$V1=="A" & phasingrefref$V2=="T"]/all_t, digits=1), col="deepskyblue4")
text(50, 50, labels=round(100*phasingrefref$V3[phasingrefref$V1=="C" & phasingrefref$V2=="T"]/all_t, digits=1), col="deepskyblue4")
text(70, 50, labels=round(100*phasingrefref$V3[phasingrefref$V1=="G" & phasingrefref$V2=="T"]/all_t, digits=1), col="deepskyblue4")
text(90, 50, labels=round(100*phasingrefref$V3[phasingrefref$V1=="T" & phasingrefref$V2=="T"]/all_t, digits=1), col="deepskyblue4")


########### Plot base change matrix ###################################################

# Plot Matrix
plot(c(0, 160), c(0, 160), type="n", xlab="", ylab="", axes=F, main="Base comparison")
for (i in 2:6) {
# horizontal
lines(c(10, 140), c(i*20,i*20), col="grey")
}
for (i in 1:6) {
# vertical
lines(c(i*20, i*20), c(30,130), col="grey") 
}
# horizontal
text(10, 110, labels="A", col="black")
text(10, 90, labels="C", col="black")
text(10, 70, labels="G", col="black")
text(10, 50, labels="T", col="black")
# vertical
text(30, 130, labels="A", col="black")
text(50, 130, labels="C", col="black")
text(70, 130, labels="G", col="black")
text(90, 130, labels="T", col="black")
text(110, 130, labels="-", col="black")
text(130, 130, labels="N", col="black")

# fill matrix
# set max vals
all_a<-sum(basecomp$V3[basecomp$V1=="A"])
all_c<-sum(basecomp$V3[basecomp$V1=="C"])
all_g<-sum(basecomp$V3[basecomp$V1=="G"])
all_t<-sum(basecomp$V3[basecomp$V1=="T"])

text(30, 110, labels=round(100*basecomp$V3[basecomp$V1=="A" & basecomp$V2=="A"]/all_a, digits=3), col="deepskyblue4")
text(50, 110, labels=round(100*basecomp$V3[basecomp$V1=="A" & basecomp$V2=="C"]/all_a, digits=3), col="deepskyblue4")
text(70, 110, labels=round(100*basecomp$V3[basecomp$V1=="A" & basecomp$V2=="G"]/all_a, digits=3), col="deepskyblue4")
text(90, 110, labels=round(100*basecomp$V3[basecomp$V1=="A" & basecomp$V2=="T"]/all_a, digits=3), col="deepskyblue4")
text(110, 110, labels=round(100*basecomp$V3[basecomp$V1=="A" & basecomp$V2=="-"]/all_a, digits=3), col="deepskyblue4")
if (length(round(100*basecomp$V3[basecomp$V1=="A" & basecomp$V2=="N"]/all_a, digits=1)) != 0) {
	text(130, 110, labels=round(100*basecomp$V3[basecomp$V1=="A" & basecomp$V2=="N"]/all_a, digits=3), col="deepskyblue4")
} else {
	text(130, 110, labels="0", col="deepskyblue4")
}


text(30, 90, labels=round(100*basecomp$V3[basecomp$V1=="C" & basecomp$V2=="A"]/all_c, digits=3), col="deepskyblue4")
text(50, 90, labels=round(100*basecomp$V3[basecomp$V1=="C" & basecomp$V2=="C"]/all_c, digits=3), col="deepskyblue4")
text(70, 90, labels=round(100*basecomp$V3[basecomp$V1=="C" & basecomp$V2=="G"]/all_c, digits=3), col="deepskyblue4")
text(90, 90, labels=round(100*basecomp$V3[basecomp$V1=="C" & basecomp$V2=="T"]/all_c, digits=3), col="deepskyblue4")
text(110, 90, labels=round(100*basecomp$V3[basecomp$V1=="C" & basecomp$V2=="-"]/all_c, digits=3), col="deepskyblue4")
if (length(round(100*basecomp$V3[basecomp$V1=="C" & basecomp$V2=="N"]/all_c, digits=1)) != 0) {
	text(130, 90, labels=round(100*basecomp$V3[basecomp$V1=="C" & basecomp$V2=="N"]/all_c, digits=3), col="deepskyblue4")
} else {
	text(130, 90, labels="0", col="deepskyblue4")
}

text(30, 70, labels=round(100*basecomp$V3[basecomp$V1=="G" & basecomp$V2=="A"]/all_g, digits=3), col="deepskyblue4")
text(50, 70, labels=round(100*basecomp$V3[basecomp$V1=="G" & basecomp$V2=="C"]/all_g, digits=3), col="deepskyblue4")
text(70, 70, labels=round(100*basecomp$V3[basecomp$V1=="G" & basecomp$V2=="G"]/all_g, digits=3), col="deepskyblue4")
text(90, 70, labels=round(100*basecomp$V3[basecomp$V1=="G" & basecomp$V2=="T"]/all_g, digits=3), col="deepskyblue4")
text(110, 70, labels=round(100*basecomp$V3[basecomp$V1=="G" & basecomp$V2=="-"]/all_g, digits=3), col="deepskyblue4")
if (length(round(100*basecomp$V3[basecomp$V1=="G" & basecomp$V2=="N"]/all_g, digits=1)) != 0) {
	text(130, 70, labels=round(100*basecomp$V3[basecomp$V1=="G" & basecomp$V2=="N"]/all_g, digits=3), col="deepskyblue4")
} else {
	text(130, 70, labels="0", col="deepskyblue4")
}

text(30, 50, labels=round(100*basecomp$V3[basecomp$V1=="T" & basecomp$V2=="A"]/all_t, digits=3), col="deepskyblue4")
text(50, 50, labels=round(100*basecomp$V3[basecomp$V1=="T" & basecomp$V2=="C"]/all_t, digits=3), col="deepskyblue4")
text(70, 50, labels=round(100*basecomp$V3[basecomp$V1=="T" & basecomp$V2=="G"]/all_t, digits=3), col="deepskyblue4")
text(90, 50, labels=round(100*basecomp$V3[basecomp$V1=="T" & basecomp$V2=="T"]/all_t, digits=3), col="deepskyblue4")
text(110, 50, labels=round(100*basecomp$V3[basecomp$V1=="T" & basecomp$V2=="-"]/all_t, digits=3), col="deepskyblue4")
if (length(round(100*basecomp$V3[basecomp$V1=="T" & basecomp$V2=="N"]/all_t, digits=1)) != 0) {
	text(130, 50, labels=round(100*basecomp$V3[basecomp$V1=="T" & basecomp$V2=="N"]/all_t, digits=3), col="deepskyblue4")
} else {
	text(130, 50, labels="0", col="deepskyblue4")
}

########### Plot error per reads ######################################################

errornumreadperc<-round((errornumread$V2/sum(errornumread$V2))*100, digits=2)

plot(errornumread$V1, errornumread$V2, type="b", lwd=3, col="orange3", axes=F, ylab="Fraction [%]", xlab="Mismatches/errors per read", main="Mismatches and errors per Read", ylim=c(0, (1.05 * max(errornumread))), xlim=c(0, 4))

abline(h=errornumread$V2[1], lty=2, col="darkgrey")
abline(h=errornumread$V2[2], lty=2, col="darkgrey")
abline(h=errornumread$V2[3], lty=2, col="darkgrey")
abline(h=errornumread$V2[4], lty=2, col="darkgrey")
abline(h=errornumread$V2[5], lty=2, col="darkgrey")

axis(1, at=c(0,1,2,3,4), labels=c("0","1","2","3","4"))
axis(2, at=c(0, errornumread$V2), labels=c("0", errornumreadperc), las=2)

########### Error by position ##########################################################

supp<-errorpos$V2
rel<-errorpos$V4

plot(errorpos$V1, rel, type="n", ylim=c(0, max(rel)+0.005), xlim=c(1, length(errorpos$V1)), axes=F, ylab="Percent errors [%]", xlab="Position in read", main="Errors in read positions")
points(errorpos$V1, rel, pch=4, lwd=3, col="firebrick")
axis(1, at=seq(1, length(errorpos$V1)))
axis(2, at=c(0, 0.01, 0.02, 0.03, 0.04, 0.05, max(rel)), labels=c("0", "1", "2", "3", "4", "5", round(max(rel)*100, 2)))

abline(h=0.01, lty=2, col="darkgrey")
abline(h=0.02, lty=2, col="darkgrey")
abline(h=0.03, lty=2, col="darkgrey")
abline(h=0.04, lty=2, col="darkgrey")
abline(h=0.05, lty=2, col="darkgrey")
abline(h=max(rel), lty=2, col="darkgrey")

text(errorpos$V1, rel, adj=c(-.5,.7), labels=supp, cex=0.65, srt=90, col="black")

legend("left", legend=c("Percentage of observed errors"), bty="n", fill=c("firebrick"))

############### Plot first error occurrence distribution ################################
occhist<-rep(0,35)
for (i in 1:36) {
occhist[i]<-length(fstocc[fstocc == i])
}

plot(errorpos$V1, occhist, type="b",  xlim=c(1, length(errorpos$V1)), axes=F, xlab="Position in read", ylab="Number", main="Occurrence of the first error", col="orange3", lwd=3)
axis(1, at=seq(1, length(errorpos$V1)))
axis(2)


############### Plot distribution of distances of the errors in the reads ###############

disthist2<-rep(0,35)
for (i in 1:35) {
disthist2[i]<-length(dist2[dist2 == i])
}
disthist3<-rep(0,35)
for (i in 1:35) {
disthist3[i]<-length(dist3[dist3 == i])
}
disthist4<-rep(0,35)
for (i in 1:35) {
disthist4[i]<-length(dist4[dist4 == i])
}
disthist2rand<-rep(0,35)
for (i in 1:35) {
disthist2rand[i]<-length(dist2rand[dist2rand == i])
}
disthist3rand<-rep(0,35)
for (i in 1:35) {
disthist3rand[i]<-length(dist3rand[dist3rand == i])
}
disthist4rand<-rep(0,35)
for (i in 1:35) {
disthist4rand[i]<-length(dist4rand[dist4rand == i])
}
disthistmax<-rep(0,35)
for (i in 1:35) {
disthistmax[i]<-length(distmax[distmax == i])
}
disthistmaxrand<-rep(0,35)
for (i in 1:35) {
disthistmaxrand[i]<-length(distmaxrand[distmaxrand == i])
}

layoutmat<-matrix(data=c(1,2,3,4), nrow=2, ncol=2)
layout(layoutmat)

# Plot 2MM
plot(1:35, disthist2, axes=F, xlab="Basepair distance between errors", main="Error cumulation (2MM)", ylab="Observations", col="darkorange2", type="p", lwd=3, xlim=c(1,35), ylim=c(0, max(max(disthist2), max(disthist2rand))))
axis(1)
axis(2)
points(1:35, disthist2rand, col="darkgrey", lwd=3) 
legend("topright", legend=c("Reads with 2 Errors", "Random distribution"), col=c("darkorange2", "darkgrey"), bty="n", lwd=3)

# Plot 3MM
plot(1:35, disthist3, axes=F, xlab="Basepair distance between errors", main="Error cumulation (3MM)", ylab="Observations", col="darkorange3", type="p", lwd=3, xlim=c(1,35), ylim=c(0, max(max(disthist3), max(disthist3rand))))
axis(1)
axis(2)
points(1:35, disthist3rand, col="darkgrey", lwd=3) 
legend("topright", legend=c("Reads with 3 Errors", "Random distribution"), col=c("darkorange3", "darkgrey"), bty="n", lwd=3)

# Plot 4MM
plot(1:35, disthist4, axes=F, xlab="Basepair distance between errors", main="Error cumulation (4MM)", ylab="Observations", col="darkorange4", type="p", lwd=3, xlim=c(1,35), ylim=c(0, max(max(disthist4), max(disthist4rand))))
axis(1)
axis(2)
points(1:35, disthist4rand, col="darkgrey", lwd=3) 
legend("topright", legend=c("Reads with 4 Errors", "Random distribution"), col=c("darkorange4", "darkgrey"), bty="n", lwd=3)

# Plot Max
plot(1:35, disthistmax, axes=F, xlab="Basepair distance between errors", main="Error cumulation, max distance (2MM, 3MM, 4MM)", ylab="Observations", col="darkred", type="p", xlim=c(1,35), lwd=3, ylim=c(0, max(max(disthistmaxrand), max(disthistmax))))
axis(1)
axis(2)
points(1:35, disthistmaxrand, col="darkgrey", lwd=3)
legend("topright", legend=c("Reads with 4 Errors","Random distribution"), col=c("darkorange4", "darkgrey"), bty="n", lwd=3)

dev.off()
