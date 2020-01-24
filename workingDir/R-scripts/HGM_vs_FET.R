#============================================
#	Script name: HGM_vs_FET.R
#	Author: Jianying LI
#	History: initially coded 06/15/2012
#		   modified on 02/25/2013
#	Comment: for IB pathway short course
##===========================================

##==================================================
#	Hypergeometric distribution -- first example
##==================================================
m <- 50  ## Number of red balls
n <- 120 ## Number of white balls
k <- 36  ## Number of drawn balls
x <- 20  ## Number of red balls in the drawing

##===========Here is the probability to get 20 red balls
dhyper(x,m,n,k) 


##=========To get a density plot================
p.list <- c()
for (i in 0 : k)
{
	p.list[i] <- dhyper(i,m,n,k)
}
plot(p.list, main="Probablity getting # of success!")


##============================================
#	Example in the presentation
##============================================

x <- 5   #num_of_DEG in GO
m <- 20  #num_of_gene on chip in GO
n <- 500 #num_of_gene on chip NOT in GO
k <- 40  #num_of_DEG

sample.DEG.testing <- matrix(c(x, (k-x), k, m -x, n-k +x, n+m-k, m, n, m+n), 
		nrow=3,
		dimnames =  list (GO_info = c("In GO","Not in GO", "Total"),
					MicroarrayInfo = c ("DEG","notDEG", "onChip")))
sample.DEG.testing

#==================================================================================
#			   DEG	Not DEGs 			totals 
#In a GO category	  x(5) 	m - x(15) 			  m (20)
#Not in GO category k -x (35)	n - k + x(465)		  n (500)
#totals 		  k(40)     (n + m) - k(480) 	        m + n (genes on array, 520)
#==================================================================================

##============================================================
#	Hypergeometric test
##============================================================

phyper (x-1, m, n, k, lower.tail=FALSE)


##=================================
#	Fisher Exact Test (FET)
##=================================
(fisher.test(matrix(c(5,35,15,465),2,2), alternative='greater'))$p.value
(fisher.test(matrix(c(x, (k-x), (m-x), (n-k+x)), 2, 2), alternative='greater'))$p.value


##=================================
#	Matches the HMT test
##=================================
phyper (x-1, m, n, k, lower.tail=FALSE)

##==========================================
#	It doesn't matter which why you do...
##==========================================
(fisher.test(matrix(c(x,(m-x),(k-x),(n-k+x)),2,2), alternative='greater'))$p.value
(fisher.test(matrix(c(x,(k-x),(m-x),(n-k+x)),2,2), alternative='greater'))$p.value
chisq.test(matrix(c(x,(k-x),(m-x),(n-k+x)),2,2), correct=FALSE)$p.value
##=========Should be the end of story==================
