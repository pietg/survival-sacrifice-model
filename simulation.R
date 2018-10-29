rm(list=ls())
library(Rcpp)
sourceCpp("primal_dual.cpp")
NumIt <-1000
n <- 1000

meanMat <- matrix(0, nrow= NumIt, ncol= 2)
timeMat <- matrix(0, nrow= NumIt, ncol= 1)
colnames(meanMat) <- c("mean1", "mean2")
colnames(timeMat) <- "IP"

t<-rep(NA,n)
delta1<-rep(NA,n)
delta2<-rep(NA,n)

for (j in 1: NumIt)
{	
  sim = 101+j
   
  set.seed(sim)
  
  print(j)
  
  ## Data generation  
  for(i in (1:n))
  {
  	x<-rexp(1,rate=0.5)
  	y<-rexp(1)+x
  	t[i]<-rexp(1,rate=0.4)
  	
  	if ((x>t[i]) && (y>t[i]))
	  {
		  delta1[i]=0;
		  delta2[i]=0;
	  }
	  else
	  {
		  if ((x<=t[i]) && (y>t[i]))
		  {
			  delta1[i]=1;
			  delta2[i]=0;
		  }
		  else
		  {
			  delta1[i]=1;
			  delta2[i]=1;
			  t[i]=y;
		  }
	  }
  }
  
  A = matrix(c(t,delta1,delta2),n,3, byrow = FALSE)

  ### Estimation
  ###########################################
  starter_IP = proc.time()
  IP <- primal_dual(A)
  mean_IP = IP$mean
  time_IP = (proc.time() -starter_IP)[3]
 
 meanMat[j,] <- mean_IP
 timeMat[j,] <- time_IP
 
  write(c(mean_IP,time_IP),file = "primal_dual_results_means.txt",ncol =3,append = TRUE)
 
  ###########################################

}

pdf("BoxPlot_means_and_time.pdf")
boxplot(meanMat, main= "Boxplot of means", las=2)
boxplot(timeMat, main="Run Times", las=2) 
 
dev.off()
