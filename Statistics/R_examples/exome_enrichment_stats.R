enriched<-read.table("exome_enriched.txt")
depleted<-read.table("exome_depleted.txt")
enriched_count<-read.table("exome_count_enriched.txt")
depleted_count<-read.table("exome_count_depleted.txt")
enrichment<-round(sum(enriched_count$V1)/(sum(enriched_count$V1)+sum(depleted_count$V1)) * 100, digits = 2)

total<-sum(enriched$V1 >= -1.0)
one<-sum(enriched$V1 >= 1.0)
five<-sum(enriched$V1 >= 5.0)
ten<-sum(enriched$V1 >= 10.0)
twenty<-sum(enriched$V1 >= 20.0)


pdf(file="SureSelect_Enrichment_Statistics.pdf", width=10, height=12)
layoutmat<-matrix(data=c(1,2), nrow=2, ncol=1, byrow=TRUE)

layout(layoutmat)


hist(round(depleted$V1), xlim=c(0,300), ylim=c(0,4000), breaks=1*0:4000, col = "red", xlab = "Mean coverage", main = paste(
	"SureSelect depleteded regions\n mean coverage: ", 
	round(mean(depleted$V1), digits = 2) ) 
)


hist(round(enriched$V1), xlim=c(0,300), ylim=c(0,4000), breaks=1*0:4000, col = "blue", xlab = "Mean coverage", main = paste(
	"SureSelect enriched regions",
	"\n mean coverage: ", round(mean(enriched$V1), digits = 2), 
	", Specificity: ",  enrichment, "%", 
	"\n Targets: ",  total, 
	", Cov >= 1: ", one, 
	", Cov >= 5: ", five, 
	", Cov >= 10: ", ten, 
	", Cov >= 20: ", twenty ) 
)
dev.off()

