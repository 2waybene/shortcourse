#==============================
#	pathway-scripts.R
##=============================

##=====================================
#	Hypergeometric distribution
##=====================================
##======A toy example===============

m <- 50  ## Number of red balls
n <- 120 ## Number of white balls
k <- 36  ## Number of drawn balls
x <- 20  ## Number of red balls in the drawing

##===========Here is the probability...
dhyper(x,m,n,k)


##==========Get a plot================
p.list <- c()
for (i in 0 : k)
{
	p.list[i] <- dhyper(i,m,n,k)
}
plot(p.list, main="Probablity getting # of success!", xlab = "Number of red ball", ylab="Probability")


##===============Do a hypergeometric testing, is it correct??==================
prob.total <- 0
for (i in 0 : x)
{
	prob.total = prob.total + dhyper(i,m,n,k)
}
1- prob.total
 
##How significant to get x SUCCESS??
1-  phyper (x, m, n, k)
phyper (x, m, n, k, lower.tail=FALSE)




##Well, it ONLY gives the P[X>x], not equivalent to a "formal testing" 

#So, we need this instead
phyper (x-1, m, n, k, lower.tail=FALSE)

##================Okay, let's try this=================================

(fisher.test(matrix(c(20,16,30,104),2,2), alternative='greater'))$p.value
(fisher.test(matrix(c(x,(k-x),(m-x),(n-k+x)),2,2), alternative='greater'))$p.value
phyper ((x -1), m, n, k, lower.tail=FALSE)

##=====Yeh, they are matched!!===========

##============================================
#	Application in microarry experiment
##=================================================



#=======This is normally what we are dealing with=======
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

dhyper (x, m, n, k) # Density functon (PDF)
1-phyper (x, m, n, k) # Distribution function (CDF)

##=====OR===================
phyper (x-1, m, n, k, lower.tail=FALSE)



##=================================
#	Fisher Exact Test (FET)
##=================================
(fisher.test(matrix(c(5,35,15,465),2,2), alternative='greater'))$p.value
(fisher.test(matrix(c(x, (k-x), (m-x), (n-k+x)), 2, 2), alternative='greater'))$p.value

phyper (x-1, m, n, k, lower.tail=FALSE)

##==========================================
#	It doesn't matter which why you do...
##==========================================
(fisher.test(matrix(c(x,(m-x),(k-x),(n-k+x)),2,2), alternative='greater'))$p.value
(fisher.test(matrix(c(x,(k-x),(m-x),(n-k+x)),2,2), alternative='greater'))$p.value

##=========Should be the end of story==================
