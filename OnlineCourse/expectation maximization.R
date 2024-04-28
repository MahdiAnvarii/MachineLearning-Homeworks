# Expectation Maximization
library(MASS)
library(mvtnorm)
library(gplots)
library(ggplot2)
library(reshape2)
s1 <- data.frame(rmvnorm(2000,c(1,10),matrix(c(6,3.5,3.5,5),nrow=2)),Cl="A")
s2 <- data.frame(rmvnorm(3000,c(3,8),matrix(c(8,-6,-6,9),nrow=2)),Cl="B")
s3 <- data.frame(rmvnorm(500,c(6,6),matrix(c(3,0,0,1),nrow=2)),Cl="C")
plot(s1[,1:2])
data <- setNames(rbind(s1,s2,s3),c("X","Y","Cl"))
ggplot(data,aes(X,Y,color=Cl))+geom_point()+scale_color_manual(values = c("black","red","blue"))+theme_bw()
given_data <- as.matrix(data[sample(1:nrow(data)),1:2])
plot(given_data)

nclass <- 3
d <- 2
# parameters
mu <- lapply(seq(nclass),function(i) given_data[sample(seq(nrow(given_data)),1),])
mu
sigma <- lapply(seq(nclass),function(i) diag(d))
sigma
phi <- rep(1/nclass,nclass)
phi

for(iteration in seq(50)){
  # Expectation
  w <- matrix(0,nrow=nrow(given_data),ncol=nclass)
  for(i in seq(nrow(given_data)))
    for(j in seq(nclass))
        w[i,j] <- dmvnorm(given_data[i,],mu[[j]],sigma[[j]]) * phi[j]
  w <- w/ rowSums(w)
  
  # Maximization
  mu <- lapply(seq(nclass),function(Cl) colSums(w[,Cl]*given_data)/sum(w[,Cl]))
  sigma <- lapply(seq(nclass),function(Cl) cov.wt(given_data,w[,Cl])$cov)
  phi <- colMeans(w)
  

}
# Visualization
vis <- setNames(data.frame(given_data,w),c("X","Y",paste0("Cl",1:nclass)))
vis.m <- melt(vis,id.vars=c("X","Y"))
print(ggplot(vis.m,aes(X,Y,color=value))+geom_point(alpha=0.3)+scale_color_gradient2(low="blue",high="red",mid="black",midpoint = 0.5)+facet_grid(.~variable)+theme_bw())


Cls <- vis[,3:5]
adad <- vis[,1:2]
j <- as.matrix(max.col(Cls))
ma <- matrix(0,nrow=nrow(Cls),ncol=1)
ll=1
for(z in seq(nrow(Cls))){
  if(j[ll]==1)
    ma[ll,1] <- "A"
  if(j[ll]==2)
    ma[ll,1] <- "B" 
  if(j[ll]==3)
    ma[ll,1] <- "C"
  ll <- ll+1
}
colnames(ma) <- "Cl"
result <- cbind(adad,ma)

ggplot(result,aes(X,Y,color=Cl))+geom_point()+scale_color_manual(values = c("black","red","blue"))+theme_bw()