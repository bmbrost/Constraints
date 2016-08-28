rm(list=ls())
library(sp)
library(raster)
library(rgeos)
library(mvtnorm)
	
setwd("~/Documents/git/Constraints")
ls()


#######################################################################
### Simulate irregular spatial support
#######################################################################

S <- SpatialLines(list(Lines(list(Line(matrix(runif(20,min=0.2,max=0.8),ncol=2))),1)))
# S <- SpatialPoints(matrix(runif(20,min=0.2,max=0.8),ncol=2))
S <- gBuffer(S,byid=TRUE,width=0.05)
# S <- gUnaryUnion(S)
plot(S)

S <- rasterize(S,raster(xmn=0,xmx=1,ymn=0,ymx=1,res=0.005))
plot(S)


#######################################################################
### Simulate true animal locations
#######################################################################

T <- 1000  # number of locations
idx <- which(values(S)==1)
mu <- SpatialPoints(xyFromCell(S,sample(idx,T,replace=TRUE)))
points(mu,pch=3,cex=0.25,col=rgb(0,0,0,0.25))


#######################################################################
### Simulate observed telemetry locations
#######################################################################

# Define observation model parameters
sigma <- c(0.01,0.05,0.1)
a <- c(0.9,0.75,0.5)
rho <- c(0,0.2,0.5)
nu <- c(30,10,2)

# Simulate complex telemetry error
s <- mu
s$lc <- sample(3,T,replace=TRUE)  # define error class (e.g., Argos location quality)

for(i in unique(s$lc)){
	Sigma <- sigma[i]^2*matrix(c(1,sqrt(a[i])*rho[i],sqrt(a[i])*rho[i],a[i]),2)
	Sigma.tilde <- matrix(c(-1,0,0,1),2)%*%Sigma%*%t(matrix(c(-1,0,0,1),2))
	idx <- which(s$lc==i)
	n <- length(idx)
	t.idx <- rbinom(n,1,prob=0.5)  # index for mixture component
	s@coords[idx[t.idx==0],] <- s@coords[idx[t.idx==0],]+rmvt(n-sum(t.idx),Sigma,nu[i])
	s@coords[idx[t.idx==1],] <- s@coords[idx[t.idx==1],]+rmvt(sum(t.idx),Sigma,nu[i])
}

plot(S,col="gray85")
points(mu,pch=3,cex=0.5,col=s$lc)
points(s,col=s$lc,cex=0.75)
segments(mu@coords[,1],mu@coords[,2],s@coords[,1],s@coords[,2],s$lc,lwd=0.5)


#######################################################################
### Fit model
#######################################################################

priors <-  list(sigma=c(0,1),a=c(0,1),rho=c(0,1),nu=c(0,100))  # prior distribution parameters
tune <- list(sigma=c(0.001,0.005,0.013),rho=c(0.06,0.22,0.13),a=c(0.1,0.1,0.14),nu=c(19,53,0.65),
	mu=c(0.01,0.05,0.08))  # tuning parameters
start <- list(mu=mu@coords,sigma=sigma,a=a,rho=rho,nu=nu)  # starting values

source("constraints.mcmc.R")
out1 <- constraints.mcmc(s,S,priors,start,tune,n.mcmc=10000,adapt=TRUE)	
out2 <- constraints.mcmc(s,S,priors,out1$end,out1$tune,n.mcmc=10000,adapt=TRUE)	
out3 <- constraints.mcmc(s,S,priors,out2$end,out2$tune,n.mcmc=10000,adapt=TRUE)	
out4 <- constraints.mcmc(s,S,priors,out3$end,out3$tune,n.mcmc=10000)	

mod <- out1
mod$keep
mod$tune
matplot(mod$sigma,type="l",lty=1);abline(h=sigma,col=1:5,lty=3)
matplot(mod$a,type="l",lty=1);abline(h=a,col=1:5,lty=3)
matplot(mod$rho,type="l",lty=1);abline(h=rho,col=1:5,lty=3)
matplot(mod$nu,type="l",lty=1);abline(h=nu,col=1:5,lty=3)

idx <- 3
s$lc[idx]
image(S,col="gray85")
points(mod$mu[idx,1,],mod$mu[idx,2,],pch=19,col=rgb(0,0,0,0.05),cex=0.25)
points(mu@coords[idx,1],mu@coords[idx,2],col=2,pch=3)
points(s@coords[idx,1],s@coords[idx,2],col=3,pch=3,xpd=TRUE)

