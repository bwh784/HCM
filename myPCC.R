library(WGCNA)
datExpr<-read.table("HCM_Filtered3CountData.txt")
datExpr_hcm<-datExpr[,1:18]
datExpr_ctrl<-datExpr[,19:23]

hcm<-cbind(datExpr_hcm[,1],datExpr_ctrl)

B<-transposeBigData(datExpr_ctrl)
ctrl_A<-adjacency(B,
selectCols = NULL,
type = "unsigned",
power = 1,
corFnc = "cor", corOptions = list(use = "p", method='pearson'),
weights = NULL)

hcm_B<-transposeBigData(hcm)
hcm_A<-adjacency(hcm_B,
selectCols = NULL,
type = "unsigned",
power = 1,
corFnc = "cor", corOptions = list(use = "p", method='pearson'),
weights = NULL)

for (i in 1:16456)
{
	for (j in (i+1):16457)
	{
		delta_p=hcm_A[i,j]-ctrl_A[i,j]
		dominator=(1-ctrl_A[i,j]*ctrl_A[i,j])/4
		x=c()
		p=pnorm(abs(delta_p),mean=0,sd=dominator,lower.tail = FALSE, log.p = FALSE)
		if(!is.na(p)&p<3.69e-10&ctrl_A[i,j]<1)
		{		
			x=cbind(x,i,j,p,delta_p,ctrl_A[i,j])
			write.table(x, file="HCM1_significant_deltaP.txt",sep = "\t",append = TRUE, row.names = FALSE,col.names = FALSE)		
		}
		
	}
}