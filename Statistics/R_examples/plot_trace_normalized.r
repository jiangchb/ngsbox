
library(RMySQL)
library(DBI)

mgr <- dbDriver("MySQL")
con <- dbConnect(mgr, user='dummy', password='stupid', dbname='chip', host='pou.crg.es')



############################################################################################################################
args <- commandArgs()
chr<-as.numeric(args[5])
end<-as.numeric(args[7])
start<-as.numeric(args[6]) # 1 = all

start = ((chr*100000000) + start)
end = ((chr * 100000000) + end)

pdf(file=paste("trace_", start, "_", end, ".pdf", sep=""), width=12, height=8)

layoutmat<-matrix(data=c(1,2), nrow=2, ncol=1)
layout(layoutmat)

ecotype = 'hyb_col_normalized'


############## A Calling ##############

query = paste("select A_A_peak, A_C_peak, A_G_peak, A_T_peak from ",ecotype," where id between ",start," and ",end, sep = '')
rs <- dbSendQuery(con, query)
ecotype_data <- fetch(rs, n = -1)

query = paste("select base from seq_ref where id between ",start," and ",end, sep = '')
rs <- dbSendQuery(con, query)
ref <- fetch(rs, n = -1)
ref = unlist(ref)

A_A_peak	= unlist(ecotype_data$A_A_peak ) 
A_C_peak	= unlist(ecotype_data$A_C_peak ) 
A_G_peak	= unlist(ecotype_data$A_G_peak ) 
A_T_peak	= unlist(ecotype_data$A_T_peak ) 

max_value = 1300 #max(c(A_A_peak,A_C_peak,A_G_peak,A_T_peak))

plot(1,1,axes=F,type='n',xlim = c(start-10,end+10), ylim = c(-30,110), xlab="", ylab="")

points(seq(start,end,1), (A_A_peak/max_value)*100, col = 'green',type = 'l')
points(seq(start,end,1), (A_C_peak/max_value)*100, col = 'blue',type = 'l')
points(seq(start,end,1), (A_G_peak/max_value)*100, col = 'black',type = 'l')
points(seq(start,end,1), (A_T_peak/max_value)*100, col = 'red',type = 'l')

base_range = seq(start,end,1)

for (i in 1:length(seq(start,end,1))) {

	text(base_range[i],0,ref[i], cex = 1.4)

}


############## Z Calling ##############

query = paste("SELECT Z_A_peak, Z_C_peak, Z_G_peak, Z_T_peak from ",ecotype," where id between ",start," and ",end, sep = '')
rs <- dbSendQuery(con, query)
ecotype_data <- fetch(rs, n = -1)

query = paste("select base from seq_ref where id between ",start," and ",end, sep = '')
rs <- dbSendQuery(con, query)
ref <- fetch(rs, n = -1)
ref = unlist(ref)

Z_A_peak        = unlist(ecotype_data$Z_A_peak ) 
Z_C_peak        = unlist(ecotype_data$Z_C_peak ) 
Z_G_peak        = unlist(ecotype_data$Z_G_peak ) 
Z_T_peak        = unlist(ecotype_data$Z_T_peak ) 

max_value = 1300 #max(c(Z_A_peak,Z_C_peak,Z_G_peak,Z_T_peak))

plot(1,1,axes=F,type='n',xlim = c(start-10,end+10), ylim = c(-30,110), xlab="", ylab="")

points(seq(start,end,1), (Z_A_peak/max_value)*100, col = 'green',type = 'l')
points(seq(start,end,1), (Z_C_peak/max_value)*100, col = 'blue',type = 'l')
points(seq(start,end,1), (Z_G_peak/max_value)*100, col = 'black',type = 'l')
points(seq(start,end,1), (Z_T_peak/max_value)*100, col = 'red',type = 'l')

base_range = seq(start,end,1)

for (i in 1:length(seq(start,end,1))) {

        text(base_range[i],0,ref[i], cex = 1.4)

}


dev.off()



