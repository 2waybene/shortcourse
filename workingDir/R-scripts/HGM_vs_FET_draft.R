#============================================
#	Script name: HGM_vs_FET_draft.R
#	Author: Jianying LI
#	History: initially coded 06/15/2012
#		   modified on 02/25/2013
#	Comment: for IB pathway short course
##===========================================

##=====================================
#	Hypergeometric distribution
##=====================================
##======A toy example===============

m <- 50  ## Number of red balls
n <- 120 ## Number of white balls
k <- 36  ## Number of drawn balls
x <- 20  ## Number of red balls in the drawing

##===========Here is the probability...
dhyper(x,m,n,k) ## <-- Probability to get 20 red balls


##==========Get a plot================
p.list <- c()
for (i in 0 : k)
{
	p.list[i] <- dhyper(i,m,n,k)
}
plot(p.list, main="Probablity getting # of success!")


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
##=================================
#	Matches the HMT test
##=================================
phyper (x-1, m, n, k, lower.tail=FALSE)

##==========================================
#	It doesn't matter which why you do...
##==========================================
(fisher.test(matrix(c(x,(m-x),(k-x),(n-k+x)),2,2), alternative='greater'))$p.value
(fisher.test(matrix(c(x,(k-x),(m-x),(n-k+x)),2,2), alternative='greater'))$p.value

##===========================================================================
#	It is going to fail if one used chi-square
#     The warning was not correct, since there was no number smaller than 5
#	What is going on here??
##==========================================================================
chisq.test(matrix(c(x,(k-x), (m-x), (n-k +x)), 2, 2), correct = FALSE)$p.value
chisq.test(matrix(c(x,(k-x), (m-x), (n-k +x)), 2, 2))$p.value


##================================================
## Now, let's take a look at Grace's example
##================================================

x <- 10   ## Find a tumor was a "success!!"
m <- 50  
n <- 50 
k <- 13  

sample.DEG.testing <- matrix(c(x, (k-x), k, m -x, n-k +x, n+m-k, m, n, m+n), 
		nrow=3,
		dimnames =  list (Treatment = c("Control","Treated", "Total"),
					TumorStatus = c ("Tumor","No Tumor", "onChip")))
sample.DEG.testing

#==================================================================================
#		  Tumor	No Tumor 			totals 
#Control	  x(3) 	m - x(47) 			  m (50)
#Treated      k -x (10)	n - k + x(40)		  n (50)
#Totals 	  k(13)     (n + m) - k(87) 	        m + n (100)
#==================================================================================
(fisher.test(matrix(c(x,(m-x),(k-x),(n-k+x)),2,2), alternative='greater'))$p.value 
chisq.test(matrix(c(x,(k-x), (m-x), (n-k +x)), 2, 2), correct = FALSE)$p.value
chisq.test(matrix(c(x,(m-x),(k-x),(n-k+x)),2,2), correct = FALSE)$p.value
phyper (x-1, m, n, k, lower.tail=FALSE)


##==================================================================================
##=======Same answers were obtained from both Fisher Exact Test and ChiSqure test==
##	Maybe need help from Grace??
##==================================================================================

##=====================================
##======A toy example===============

m <- 50  ## Number of red balls
n <- 50  ## Number of white balls
k <- 13  ## Number of drawn balls


##===========Here is the probability...
x <- 10  ## Number of red balls in the drawing
dhyper(x,m,n,k) ## <-- Probability to get 10 red balls


x <- 11  ## Number of red balls in the drawing
dhyper(x,m,n,k) ## <-- Probability to get 10 red balls

x <- 12  ## Number of red balls in the drawing
dhyper(x,m,n,k) ## <-- Probability to get 10 red balls

x <- 13  ## Number of red balls in the drawing
dhyper(x,m,n,k) ## <-- Probability to get 10 red balls


