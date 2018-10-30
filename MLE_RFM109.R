	rm(list=ls())
	library(Rcpp)
	sourceCpp("primal_dual.cpp")
	sourceCpp("EM.cpp")

	A<-read.table("RFM109.txt")
	output1 <- primal_dual(A)
	output2 <- EM(A)

	B1 <- output1$MLE1
	C1 <- output1$MLE2
	
	B2 <- output2$MLE1
	C2 <- output2$MLE2

	x<-B1[,1]
	y<-B1[,2]
	x1<-C1[,1]	
	y1<-C1[,2]
	
	t<-B2[,1]
	u<-B2[,2]
	t1<-C2[,1]	
	u1<-C2[,2]	
plot(c(-1000,-1000),xlim=c(min(x,x1),max(x,x1)),ylim=c(0,1), main= "",ylab="",xlab="",bty="n",las=1)

	lines(x,y,lwd=2,col="blue",type='s')
	lines(x1,y1,lwd=2,col="red",type='s')
	
	lines(t,u,lwd=2, lty=1,type='s')
	lines(t1,u1,lwd=2, lty=1,type='s')
 


