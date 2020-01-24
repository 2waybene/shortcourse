#================================================================
#	pathway-scripts-2018.R
# author: Jianying Li
# comment: this script accompanies the pathway analysis course
#          of summer 2018
##==============================================================

##================================================##
##  Slide 14 -- hypergenometric distribution      ##
##    A clinical trial study example              ##
##================================================##

m <- 5  ## 5 patients to have elevated C
n <- 5  ## 5 controls to have elevated C
k <- 6  ## Number of drawn individual with elevated C
x <- 4  ## Number of individuals with elevated C who are patients

dhyper(x,m,n,k)



sample.trial.testing <- t(matrix(c(x, (k-x), k, m -x, n-k +x, n+m-k, m, n, m+n),
                                 nrow=3,
                                 dimnames =  list (VolunteerInfo = c ("Patient_w_D","Healthy_Control", "Total"),
                                                   TestResults = c("Elev_C","Norm_C", "Total"))))

sample.trial.testing


##================================================##
##  Slide 16 -- hypergenometric distribution      ##
##================================================##
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


##================================================##
##  Slide 17 -- hypergenometric distribution      ##
##================================================##
m <- 5  ## 5 patients to have elevated C
n <- 5  ## 5 controls to have elevated C
k <- 6  ## Number of drawn individual with elevated C
x <- 4  ## Number of individuals with elevated C who are patients

x <- 5 ## Number of individuals with elevated C who are patients
dhyper(x,m,n,k)

##================================================##
##  Slide 18 -- hypergenometric distribution      ##
##================================================##

## Method 1:
x <- 4
p4 <- dhyper(x,m,n,k)

x <- 5
p5 <- dhyper(x,m,n,k)

p4 + p5

## Method 2:

1- phyper(3,5,5,6)
#[1] 0.2619048

dhyper(4,5,5,6) + dhyper(5,5,5,6)
#[1] 0.2619048


##====================================##
##  slies 19                          ##
##
#	Hypergeometric distribution
# Clinical trial example
##=====================================


x <- 4      #patients with elevated level of compound C
m <- 6      #all with elevated level of compound C
n <- 4      #all with normal level of compound C
k <- 5      #total number of patient volunteers


sample.trial.testing <- t(matrix(c(x, (k-x), k, m -x, n-k +x, n+m-k, m, n, m+n),
                               nrow=3,
                               dimnames =  list (VolunteerInfo = c ("Patient_w_D","Healthy_Control", "Total"),
                                                 TestResults = c("Elev_C","Norm_C", "Total"))))
sample.trial.testing

phyper((x-1), m, n, k, lower.tail=FALSE)
(fisher.test(matrix(c(x,(k-x), (m-x), (n-k+x)),2,2), alternative='greater'))$p.value



##====================================##
##  slies 20                          ##
##====================================##
#	Hypergeometric distribution
# For pathway
##======A toy example===============

m <- 50  ## Number of red balls
n <- 120 ## Number of white balls
k <- 36  ## Number of drawn balls
x <- 20  ## Number of red balls in the drawing

##===========Here is the probability...
dhyper(x,m,n,k)


##===========Here is the test...
phyper((x-1), m, n, k, lower.tail=FALSE)
(fisher.test(matrix(c(x,(k-x), (m-x), (n-k+x)),2,2), alternative='greater'))$p.value


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

dhyper   (x, m, n, k) # Density functon (PDF)
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
#	In fisher exatc test,
# it doesn't matter which why you do...
##==========================================
(fisher.test(matrix(c(x,(m-x),(k-x),(n-k+x)),2,2), alternative='greater'))$p.value
(fisher.test(matrix(c(x,(k-x),(m-x),(n-k+x)),2,2), alternative='greater'))$p.value

##=========Should be the end of story==================
