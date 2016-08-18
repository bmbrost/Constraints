###
### R script file for fitting mixture t model to simulated and harbor seal data
### Brian M. Brost (14 MAR 2015)
###

# Software Disclaimer:
# Although this software program has been used by the U.S. Geological Survey (USGS), no warranty, 
# expressed or implied, is made by the USGS or the U.S. Government as to the accuracy and functioning 
# of the program and related program material nor shall the fact of distribution constitute any 
# such warranty, and no responsibility is assumed by the USGS in connection therewith. 

# Data Disclaimer:
# Although these data have been processed successfully on a computer system at the U.S. Geological 
# Survey (USGS), no warranty expressed or implied is made regarding the display or utility of the 
# data on any other system or for general or scientific purposes, nor shall the act of distribution 
# constitute any such warranty. The USGS or the U.S. Government shall not be held liable for improper
# or incorrect use of the data described and/or contained herein.


###
### Load workspace
###

rm(list=ls())
setwd("~/Dropbox/constraints/supplement/")
load("workspace.Rdata") #Load workspace
source("constraints.MCMC.R") #Load MCMC algorithm
ls()

###
### Libraries and subroutines
###

library(sp)
library(raster)
library(mvtnorm)
library(parallel)
library(doParallel)
library(foreach)	
library(data.table)

# Note: an optional C++ routine can be used to update mu (Appendix B, step 5) more quickly. To use
# C++ for this update, load the following packages, source the C++ file, and use the argument 
# Cpp=TRUE to constraints.MCMC below. Alternatively, mu can be updated in native R using Cpp=FALSE. 
library(Rcpp)
library(inline)
library(RcppArmadillo)
Rcpp::sourceCpp("update.mu.cpp")  # C++ routine for mu update (Appendix B, step 5)

###
### Explore contents of workspace
###

plot(kodiak,col="grey")  # portion of the northeastern end of Kodiak Island, AK (gray)
points(haulout,col=2,pch=19)  # location of haul-out site
plot(S,add=TRUE)  # 500-m resolution raster representing S, the spatial support of mu
plot(d.haulout,add=TRUE)  # 500-m resolution raster for distance to haul-out site covariate
plot(bath,add=TRUE)  # 500-m resolution raster for bathymetry covariate
points(seal,pch=19,cex=0.25)  # Argos telemetry locations of harbor seal used in case study
head(X)  # design matrix: intercept, distance to haul-out, and bathymetry for all cells in S
head(X.scale)  # center and scaled design matrix
lu  # lookup table of swim distances between raster cells in S


###
### Fit model to simulated data
###

T <- 500  # number of locations to simulate

# Set parameters for observation model (corresponding to Argos location classes 3,0, and B)
sigma <- c(2291,2727,13252)
nu <- c(16.74,1.60,1.00)
a <- c(0.70,0.50,0.75)
rho <- c(0.85,0.16,0.30)

# Set parameters for process model
phi <- 424  # m/hr
beta <- c(0,-0.0001511532, -0.0270482365)  # resource selection coefficients; not standardized

# Simulate true locations with resource selection and temporal dependence
mu <- get.mu(T=T,haulout,S,lu,X,beta,phi,lc=3,a=1.1,s=4.5)  # true locations
plot(mu,pch="",axes=TRUE)  # plot true locations mu
plot(S,add=TRUE)
points(mu[,1:2],pch=19,cex=0.5)

# Simulate observed locations
s <- get.s(mu,sigma,a,rho,nu)  # observed locations
points(s[,1:2],col=(s$lc+1),pch=19,cex=0.5)  # add observed locations s

# Plot all observed locations s
plot(s[,1:2],col=(s$lc+1),pch=19,cex=0.5,axes=TRUE)  # add observed locations s

# Specify prior distribution parameters
priors <- list(sigma=c(0,20000),phi=c(0,750),a=c(0,1),rho=c(0,1),nu=c(0,30),beta.mn=rep(0,3),beta.var=diag(3)*10)

# Set tuning parameters
tune <- list(mu=c(950,1200,1250),sigma=c(400,750,2750),rho=c(0.1,0.39,0.44),a=c(0.205,0.25,0.275),
	nu=c(15.0,0.59,0.4),phi=90,beta=0.045)

# Specify starting values
start <- list(mu=mu@coords,sigma=sigma,nu=nu,a=a,rho=rho,phi=phi,beta=beta*X.sd)

# Fit model
n.cores <- 7  # number of cores used for parallel processing
out1 <- constraints.MCMC(s,S,lu,X.scale,priors,start,tune,n.mcmc=500,n.cores=n.cores,Cpp=TRUE)	
out1$keep  # proposal acceptance rate; adjust tuning parameters to achieve values between 0.2 and 0.4

# Examine chains; true parameter values are represented by dashed horizontal lines
matplot(out1$sigma,type="l",lty=1);abline(h=sigma,col=1:5,lty=3)
matplot(out1$a,type="l",lty=1);abline(h=a,col=1:5,lty=3)
matplot(out1$rho,type="l",lty=1);abline(h=rho,col=1:5,lty=3)
matplot(out1$nu,type="l",lty=1);abline(h=nu,col=1:5,lty=3)
plot(out1$phi,type="l",lty=1);abline(h=phi,lty=3)
matplot(out1$beta[,-1],type="l",lty=1);abline(h=beta*X.sd,col=0:5,lty=3)

# Posterior means
apply(out1$sigma,2,mean)
apply(out1$a,2,mean)
apply(out1$rho,2,mean)
apply(out1$nu,2,mean)
mean(out1$phi)
apply(out1$beta,2,mean)

# Posterior mode of mu
mu.mode <- apply(out1$mu,2,get.mode)
mu.hat <- xyFromCell(S,mu.mode)
plot(mu.hat,pch="")
plot(S,add=TRUE)
points(mu.hat,pch=19,cex=0.5)  # estimated mu


###
### Fit model to harbor seal data
###

# Specify prior distribution parameters
priors <- list(sigma=c(0,20000),phi=c(0,750),a=c(0,1),rho=c(0,1),nu=c(0,30),beta.mn=rep(0,3),beta.var=diag(3)*10)

# Set tuning parameters
tune <- list(mu=c(800,700,900,1000,1000,1000),
	sigma=c(950,700,450,800,750,2500),
	rho=c(0.5,0.19,0.35,0.29,0.325,0.375),
	a=c(0.375,0.275,0.175,0.13,0.18,0.21),
	nu=c(21,1.25,0.525,0.29,0.24,0.225),
	phi=75,beta=c(0.025,0.025,0.025))

# Specify starting values
start <- list(mu=mu.start,sigma=c(2300,1000,1400,2700,2700,13000),nu=c(16,2,1.8,1.6,1.2,1.0),
	a=c(0.7,0.5,0.65,0.5,0.9,0.75),rho=c(0.85,0.75,0.38,0.16,0.2,0.3),phi=400,beta=c(0,-2.1,-0.8))

# Fit model
n.cores <- 4  # number of cores used for parallel processing
out2 <- constraints.MCMC(seal,S,lu,X.scale,priors,start,tune,n.mcmc=2000,n.cores=n.cores,Cpp=TRUE)	
out2$keep  # proposal acceptance rate; adjust tuning parameters to achieve values between 0.2 and 0.4

# Examine chains
matplot(out2$sigma,type="l",lty=1)
matplot(out2$a,type="l",lty=1)
matplot(out2$rho,type="l",lty=1)
matplot(out2$nu,type="l",lty=1)
plot(out2$phi,type="l",lty=1)
matplot(out2$beta[,-1],type="l",lty=1)

# Posterior means
apply(out2$sigma,2,mean)
apply(out2$a,2,mean)
apply(out2$rho,2,mean)
apply(out2$nu,2,mean)
mean(out2$phi)
apply(out2$beta,2,mean)

# Posterior mode of mu
mu.mode <- apply(out2$mu,2,get.mode)
mu.hat <- xyFromCell(S,mu.mode)
plot(mu.hat,pch="")
plot(S,add=TRUE)
points(mu.hat,pch=19,cex=0.5)  # estimated mu
