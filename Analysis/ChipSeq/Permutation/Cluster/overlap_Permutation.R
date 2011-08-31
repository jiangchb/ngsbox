### Read in command line options ##################################################
args <- commandArgs()			# get args
real_overlap<-as.integer(args[5])	# Read real overlap
permutation<-read.table(args[6])	# Read permutation data file
opt<-read.table(args[7])		# Get permutation parameters
pdf(file=args[8], width=10, height=6)	# Open PDF file




### Modify plot appearance ########################################################
breaks=100							# define breaks in histogram

xmax = max(real_overlap + real_overlap / 2, 20)			# Calculate x-Axis length

hist<-hist(permutation$V1, plot = FALSE, breaks = breaks)	# Calculate y-Axis length
ymax=max(hist$counts)

coln="steelblue4"						# define colors
colu="goldenrod1"




### Plot historgram ################################################################
hist<-hist( permutation$V1, main="",
		xlab="Overlap: ChIPseq Peaks, Differentially Expressed Genes,  TFBS-Predictions", 
		xlim=c(0,xmax), ylim=c(0,ymax), 
		breaks=breaks, col=coln)
x<-hist$counts




### Check if data follows Poisson distribution #####################################
#POISSON <- function(k,m) (m^k) * exp(-m) / factorial(k)			# define Poisson dist
#m = mean(permutation[,1])							# get parameter lambda
#p = format(sum(POISSON(real_overlap:1000, m)), digits = 20)			# calculate pvalue
#possion_values = (POISSON(0:xmax, m))						# set and scale poisson values
#scale = ymax / max(possion_values);
#possion_values = scale * possion_values;
#lines(0:xmax, possion_values, col="darkgrey")					# plot dots for Poisson dist
#y_pos = ymax / 10								# plot pvalue
#text(real_overlap, y_pos + y_pos, paste("Real overlap: ", real_overlap, "\n", "p-value: ", p))


### Check if data follows Poisson distribution #####################################
m = round(mean(permutation[,1]), digits=0)                                      # get parameter lambda
scale = ymax / dpois(m,m)                                                       # get scaling factor for poisson dist plot
lines(0:xmax, scale * dpois(0:xmax, m), col="red")                              # plot poisson ditribution
prob=(0.5-(dpois(m,m)/2)) - sum(dpois(seq(1:(real_overlap-m-1))+m, m))          # calculate pvalue
y_pos = max(ymax / 10, 10)                                                      # plot pvalue
text(real_overlap, y_pos + y_pos, paste("Real overlap: ", real_overlap, "\n", "p-value: ", round(prob, digits=5), "\nlambda: ", m))



### Print arrow for real overlap ###################################################
arrows(real_overlap, 0, real_overlap, 10, col="darkolivegreen4", lwd=3, code=1)




### Print options ##################################################################
options = paste( paste( toString(opt[1,1]), toString(opt[1,2]), sep = "=" ),
		 paste( toString(opt[2,1]), toString(opt[2,2]), sep = "=" ),
		 paste( toString(opt[3,1]), toString(opt[3,2]), sep = "=" ),
		 paste( toString(opt[4,1]), toString(opt[4,2]), sep = "=" ),
		 paste( toString(opt[5,1]), toString(opt[5,2]), sep = "=" ),
		 #paste( toString(opt[6,1]), toString(opt[6,2]), sep = "=" ),
		 #paste( toString(opt[7,1]), toString(opt[7,2]), sep = "=" ),
		 #paste( toString(opt[8,1]), toString(opt[8,2]), sep = "=" ),
	sep = "\n")

mtext(text = options, side = 3, line = -5)


# Show first 50 warnings if any occurred
warnings()

dev.off()
