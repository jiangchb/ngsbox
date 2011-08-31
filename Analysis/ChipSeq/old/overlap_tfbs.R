pos<-read.table("permutations.txt")
res<-read.table("result.txt")
opt<-read.table("options.txt", header = FALSE)

#"range, max_cov, segment_length, permutations, tfbs_tables"
# create pdf

pdf(file=paste(paste("overlap_hist", "r", toString(opt$V2[1]), "m", toString(opt$V2[2]), "s", toString(opt$V2[3]), "t", toString(opt$V2[5]),"c", toString(opt$V2[6]), sep = "_"),".pdf", sep = ""), width=10, height=6)

# overlap ChipSeq and prediction
x_res = res[1,1];

# add 50% for visualization, but display at least 20
xmax = max(x_res + x_res / 2, 20)

# ???
breaks=100
# compute the max value
hist<-hist(pos$V1, plot = FALSE, breaks = breaks)
ymax=max(hist$counts)

# some colors
coln="steelblue4"
colu="goldenrod1"

# draw the histogram
hist<-hist(pos$V1, main="", xlab="ChipSeq overlapping TFBS-Predictions", xlim=c(0, xmax), ylim=c(0, ymax), col=coln, breaks = breaks)
x<-hist$counts
par(new=T)

# data follow poisson distribution?
POISSON <- function(k,m) (m^k) * exp(-m) / factorial(k)
# parameter lambda
m = mean(pos[,1])
# p-value
p = format(sum(POISSON(x_res:1000, m)), digits = 4)
possion_values = (POISSON(0:xmax, m))
scale = ymax / max(possion_values);
possion_values = scale * possion_values;
# print point from poisson
points(0:xmax, possion_values, col="red")

# print at 10% of the max y value
y_pos = ymax / 10
text(x_res, y_pos + y_pos, paste("Result: ", x_res, "\n", "p-value: ", p))

# print the options (name of the parameter) and (value of the parameter)
options = paste(paste(toString(opt[1,1]),toString(opt[1,2]), sep = "="), paste(toString(opt[2,1]),toString(opt[2,2]), sep = "="), paste(toString(opt[3,1]),toString(opt[3,2]), sep = "="), paste(toString(opt[4,1]),toString(opt[4,2]), sep = "="), paste(toString(opt[5,1]),toString(opt[5,2]), sep = "="), paste(toString(opt[6,1]),toString(opt[6,2]), sep = "="), sep = "\n")
mtext(text = options, side = 3, line = -5)
# print an arrow
arrows(x_res, 0, x_res, y_pos, col="darkolivegreen4", lwd=3, code=1)

par(new=F)

dev.off()
