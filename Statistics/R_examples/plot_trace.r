
library(RMySQL)
library(DBI)

mgr <- dbDriver("MySQL")
con <- dbConnect(mgr, user='dummy', password='stupid', dbname='chip', host='pou.crg.es')


par(mfcol = c(2,1))

############################################################################################################################
args <- commandArgs()
chr<-as.numeric(args[5])
end<-as.numeric(args[7])
start<-as.numeric(args[6]) # 1 = all

start = ((chr*100000000) + start)
end = ((chr * 100000000) + end)

pdf(file=paste("trace_", start, "_", end, ".pdf", sep=""), width=12, height=4)

ecotype = 'hyb_col'


############## A Calling ##############

query = paste("select A_call, A_qscore, A_A_peak, A_C_peak, A_G_peak, A_T_peak from ",ecotype," where id between ",start," and ",end, sep = '')
rs <- dbSendQuery(con, query)
ecotype_data <- fetch(rs, n = -1)

query = paste("select base from seq_ref where id between ",start," and ",end, sep = '')
rs <- dbSendQuery(con, query)
ref <- fetch(rs, n = -1)
ref = unlist(ref)

A_call		= unlist(ecotype_data$A_call )
A_qscore	= unlist(ecotype_data$A_qscore ) 
A_A_peak	= unlist(ecotype_data$A_A_peak ) 
A_C_peak	= unlist(ecotype_data$A_C_peak ) 
A_G_peak	= unlist(ecotype_data$A_G_peak ) 
A_T_peak	= unlist(ecotype_data$A_T_peak ) 

max_value = max(c(A_A_peak,A_C_peak,A_G_peak,A_T_peak))

plot(1,1,axes=F,type='n',xlim = c(start-10,end+10), ylim = c(-30,110) )

points(seq(start,end,1), (A_A_peak/max_value)*100, col = 'green',type = 'l')
points(seq(start,end,1), (A_C_peak/max_value)*100, col = 'blue',type = 'l')
points(seq(start,end,1), (A_G_peak/max_value)*100, col = 'black',type = 'l')
points(seq(start,end,1), (A_T_peak/max_value)*100, col = 'red',type = 'l')

base_range = seq(start,end,1)

for (i in 1:length(seq(start,end,1))) {

	text(base_range[i],0,ref[i], cex = 1.4)

}

for (i in 1:length(seq(start,end,1))) {
	
	shade = (32-A_qscore[i])*3
	shade = paste('grey',shade,sep='')
	
	rect(base_range[i]-0.5,-30,base_range[i]+0.5,-20, col = shade, border = shade)
}


############## Z Calling ##############

query = paste("select Z_call, Z_qscore, Z_A_peak, Z_C_peak, Z_G_peak, Z_T_peak from ",ecotype," where id between ",start," and ",end, sep = '')
rs <- dbSendQuery(con, query)
ecotype_data <- fetch(rs, n = -1)

query = paste("select base from seq_ref where id between ",start," and ",end, sep = '')
rs <- dbSendQuery(con, query)
ref <- fetch(rs, n = -1)
ref = unlist(ref)

Z_call          = unlist(ecotype_data$Z_call )
Z_qscore        = unlist(ecotype_data$Z_qscore ) 
Z_A_peak        = unlist(ecotype_data$Z_A_peak ) 
Z_C_peak        = unlist(ecotype_data$Z_C_peak ) 
Z_G_peak        = unlist(ecotype_data$Z_G_peak ) 
Z_T_peak        = unlist(ecotype_data$Z_T_peak ) 

max_value = max(c(Z_A_peak,Z_C_peak,Z_G_peak,Z_T_peak))

plot(1,1,axes=F,type='n',xlim = c(start-10,end+10), ylim = c(-30,110) )

points(seq(start,end,1), (Z_A_peak/max_value)*100, col = 'green',type = 'l')
points(seq(start,end,1), (Z_C_peak/max_value)*100, col = 'blue',type = 'l')
points(seq(start,end,1), (Z_G_peak/max_value)*100, col = 'black',type = 'l')
points(seq(start,end,1), (Z_T_peak/max_value)*100, col = 'red',type = 'l')

base_range = seq(start,end,1)

for (i in 1:length(seq(start,end,1))) {

        text(base_range[i],0,ref[i], cex = 1.4)

}

for (i in 1:length(seq(start,end,1))) {

	shade = (32-Z_qscore[i])*3
        shade = paste('grey',shade,sep='')

        rect(base_range[i]-0.5,-30,base_range[i]+0.5,-20, col = shade, border = shade)
}



dev.off()



