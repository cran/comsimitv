## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo = FALSE------------------------------------------------------
library(comsimitv)

## ----gauss, fig.cap="Shape of Gaussian competirion kernel with different $\\sigma_b$ values", echo =FALSE, fig.height=3.5----
oldpar <- par(no.readonly = TRUE)
d<-seq(-1,1,0.01)
sigma.b=0.01
alpha<-exp(-d^2/sigma.b)
par(mai=c(0.82, 0.82, 0.2, 0.42),cex.lab=1.2)
plot(d,alpha, type="l",col="red",
     ylab=expression(alpha["ij"]),
     xlab=expression("B"["i"]-"B"["j"]),
     ylim=c(0,1.2))
sigma.b=0.1
alpha<-exp(-d^2/sigma.b)
lines(d,alpha, type="l",col="blue")
sigma.b=1
alpha<-exp(-d^2/sigma.b)
lines(d,alpha, type="l",col="green")
legend("topright", c(expression(paste(sigma["b"]," = 0.01")), expression(paste(sigma["b"]," = 0.1")),expression(paste(sigma["b"]," = 1"))),
       col = c("red", "blue","green"),lty=1,bty="n")
par(oldpar)

## ----kisdi1, echo =FALSE, fig.cap="Shape of Kisdi's convex-convave function with different  values of v (C=1,$\\sigma_b$=0.1)", fig.height=3.5----
oldpar <- par(no.readonly = TRUE)
d<-seq(-1,1,0.01)
sigma.b=0.1
C=1
v=0.1
alpha<-C*(1-1/(1+v*exp(-(d)/sigma.b)))
par(mai=c(0.82, 0.82, 0.2, 0.42),cex.lab=1.2)
plot(d,alpha, type="l",col="red",
     ylab=expression(alpha["ij"]),
     xlab=expression("B"["i"]-"B"["j"]))
v=1
alpha<-C*(1-1/(1+v*exp(-(d)/sigma.b)))
lines(d,alpha, type="l",col="blue")
v=10
alpha<-C*(1-1/(1+v*exp(-(d)/sigma.b)))
lines(d,alpha, type="l",col="green")
legend("right", c("v=0.1","v=1","v=10"),
       col = c("red", "blue","green"),lty=1,bty="n")
par(oldpar)

## ----fig3, echo =FALSE, fig.cap="Shape of Kisdi's convex-convave function with different  values of $\\sigma_b$ (C=1,v=1)", fig.height=3.5----
oldpar <- par(no.readonly = TRUE) 
d<-seq(-1,1,0.01)
sigma.b=0.01
C=1
v=1
alpha<-C*(1-1/(1+v*exp(-(d)/sigma.b)))
par(mai=c(0.82, 0.82, 0.2, 0.42),cex.lab=1.2)
plot(d,alpha, type="l",col="red",
     ylab=expression(alpha["ij"]),
     xlab=expression("B"["i"]-"B"["j"]))
sigma.b=0.1
alpha<-C*(1-1/(1+v*exp(-(d)/sigma.b)))
lines(d,alpha, type="l",col="blue")
sigma.b=-0.1
alpha<-C*(1-1/(1+v*exp(-(d)/sigma.b)))
lines(d,alpha, type="l",col="green")
legend("right", c(expression(paste(sigma["b"]," = 0.01")), expression(paste(sigma["b"]," = 0.1")),expression(paste(sigma["b"]," = -0.1"))),
       col = c("red", "blue","green"),lty=1,bty="n")
par(oldpar)

## ----fig4, echo =FALSE, fig.cap="Shape of smooth function by Nattrass et al. with different  values of $\\sigma_b$ (C=1)", fig.height=4----
oldpar <- par(no.readonly = TRUE)
d<-seq(-1,1,0.01)
sigma.b=0.01
C=1
v=1
alpha<-1+C-2*C/(1+exp(-(2*d)/sigma.b))
par(mai=c(0.82, 0.82, 0.2, 0.42),cex.lab=1.2)
plot(d,alpha, type="l",col="red",
     ylab=expression(alpha["ij"]),
     xlab=expression("B"["i"]-"B"["j"]))
sigma.b=0.5
alpha<-1+C-2*C/(1+exp(-(2*d)/sigma.b))
lines(d,alpha, type="l",col="blue")
sigma.b=-0.5
alpha<-1+C-2*C/(1+exp(-(2*d)/sigma.b))
lines(d,alpha, type="l",col="green")
legend("right", c(expression(paste(sigma["b"]," = 0.01")), expression(paste(sigma["b"]," = 0.5")),expression(paste(sigma["b"]," = -0.5"))),
       col = c("red", "blue","green"),lty=1,bty="n")
par(oldpar)

