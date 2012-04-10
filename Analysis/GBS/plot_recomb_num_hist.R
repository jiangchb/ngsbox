#OVERALL:
#* VIS
#* OUTPUT PERC P1 PER WINDOW
#* RECOMB NUM HIST

list = read.table("file_listing.txt")
a=0

for (i in list$V1) { 
	con = read.table(i)
	if (a==0) {
		df = data.frame(con)
		a=1
	}
	else {
		df = data.frame(df, con)
	}
}



