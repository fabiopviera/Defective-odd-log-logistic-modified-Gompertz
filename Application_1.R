library(survival)
require(bshazard)
library(survminer)



t<-c(49, 49, 30, 10, 7, 37, 57, 17, 15, 33, 60, 42, 36, 17, 39,
     38, 35, 45, 3, 50, 52, 24, 40, 34, 18, 60, 39, 5, 24, 4, 5, 54,
     10, 14, 52, 18, 18, 22, 42, 49, 38, 12, 10, 44, 36, 26, 61, 54,
     34, 31, 33, 18, 54, 38, 29, 46, 39, 11, 59, 61, 33, 59, 49, 32,
     22, 49, 22, 41, 59, 37, 36, 34, 35, 51, 57, 42, 47, 9, 36, 16,
     43, 43, 12, 57, 54, 57, 23, 30, 16, 38, 59, 16, 41, 58, 51, 37,
     32, 34, 53, 33, 41, 53, 30, 41, 51, 50, 52, 34, 32, 39, 59, 35,
     12, 36, 61, 30, 46, 15, 57, 33, 45, 27, 48, 56)
d<-c(0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,
     1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0,
     1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
     0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1,
     0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0)


dat_ <- data.frame(time=t,status=d)

fit.km <- survfit(Surv(time,status)~1 ,data=dat_)

summary(fit.km)

plot(fit.km,lwd=3,col=c("gray2","gray50","gray50"))


require(bshazard)

fit.bshazard <- bshazard(Surv(time, status) ~ 1,data=dat_)
summary(fit.bshazard)

x11()
plot(fit.bshazard$time, fit.bshazard$hazard,
     cex.lab=1.2,xlab="Time (months)",
     ylab='Hazard function h(t)',col="gray65",
     type = "s",ylim=c(0,0.02),xlim=c(0,60),lwd=4)


status<-d


### Estimando MGo - modelo reparametrizado
verosi1_rep <- function(par){
  
  eta  <- par[1]
  beta <- par[2]
  
  lambda <- -log(eta)
  
  logS <- -lambda*(1-exp(-beta*t))
  
  lv <- logS + status*(log(lambda) + log(beta) - beta*t)
  
  return(-sum(lv))
}


valorin1 <- c(0.1,0.1)
valorin1 <- c(estima1$par[1],estima1$par[2])
valorin1 <- c(estima12$par[1],estima12$par[2])
valorin1 <- c(estima13$par[1],estima13$par[2])

(estima1 <- optim(valorin1,verosi1_rep,method = "SANN",hessian=T))
(estima12 <- optim(valorin1,verosi1_rep,method = "BFGS",hessian=T))
(estima13 <- optim(valorin1,verosi1_rep, method = "L-BFGS-B",
                   lower = rep(0.001, 2), upper = rep(1500, 2),hessian=T))


##Refinamento Final
(estima12 <- optim(valorin1,verosi1_rep,method = "BFGS",hessian=T))

(verosi11 <-- verosi1_rep(c(estima12$par[1],estima12$par[2])))

d1 <- 2

n <-length(t)

456.6365/2

(TRV.1<-2*(-456.6365/2-verosi11))
(p_valor.1<-pchisq(TRV.1,1,lower.tail = F))


(AIC_MGO <- -2*verosi11+2*d1)
(BIC_MGO <- -2*verosi11+d1*log(n))
(CAIC_MGO <- AIC_MGO +(2*(d1+2)*(d1+3))/(n-d1-3))
(Inv.he_MGO <-  solve(estima13$hessian))
(e.p_MGO<-c(sqrt(Inv.he_MGO[1,1]),sqrt(Inv.he_MGO[2,2])))

alpha0 = 0.0112 ; beta0 =0.0106
pi00 <- exp(-alpha0/beta0)



###################################3

###Estimando MOLLGO
### Estimando OLLMGo - modelo reparametrizado
verosi2_rep <- function(par){
  
  eta   <- par[1]
  beta  <- par[2]
  gamma <- par[3]
  
  # Espaço paramétrico
  if(eta <= 0 || eta >= 1 || beta <= 0 || gamma <= 0)
    return(Inf)
  
  lambda <- log(1 + ((1-eta)/eta)^(1/gamma))
  alpha  <- lambda*beta
  
  g <- alpha*exp(-beta*t)*
    exp(-lambda*(1-exp(-beta*t)))
  
  G <- 1 - exp(-lambda*(1-exp(-beta*t)))
  
  f <- (gamma*g*(G^(gamma-1))*((1-G)^(gamma-1))) /
    (G^gamma + (1-G)^gamma)^2
  
  S <- 1 - (G^gamma)/(G^gamma + (1-G)^gamma)
  
  lv <- status*log(f) + (1-status)*log(S)
  
  if(any(!is.finite(lv)))
    return(Inf)
  
  return(-sum(lv))
}
#alpha tau k lambda
valorin1 <- c(0.1,0.1,0.1)
valorin1 <- c(estima1$par[1],estima1$par[2],estima1$par[3])
valorin1 <- c(estima12$par[1],estima12$par[2],estima12$par[3])
valorin1 <- c(estima13$par[1],estima13$par[2],estima13$par[3])

(estima1 <- optim(valorin1,verosi2_rep,method = "SANN",hessian=T))
(estima12 <- optim(valorin1,verosi2_rep,method = "BFGS",hessian=T))
(estima13 <- optim(valorin1,verosi2_rep, method = "L-BFGS-B",
                   lower = rep(0.001, 2), upper = rep(1500, 2),hessian=T))


##Refinamento Final
(estima12 <- optim(valorin1,verosi2_rep,method = "BFGS",hessian=T))

(verosi22 <-- verosi2_rep(c(estima12$par[1],estima12$par[2],estima12$par[3])))

d2 <- 3

n <-length(t)


(TRV.1<-2*(verosi22-verosi11))

(p_valor.1<-pchisq(TRV.1,1,lower.tail = F))


(AIC_OLLMGO <- -2*verosi22+2*d2)
(BIC_OLLMGO <- -2*verosi22+d2*log(n))
(CAIC_OLLMGO <- AIC_OLLMGO +(2*(d2+2)*(d2+3))/(n-d2-1))
(Inv.he_OLLMGO <-  solve(estima12$hessian))
(e.p_OLLMGO<-c(sqrt(Inv.he_OLLMGO[1,1]),sqrt(Inv.he_OLLMGO[2,2]),sqrt(Inv.he_OLLMGO[3,3])))



alpha1 <- 0.04036214; beta1 <-  0.06814469; gamma1 <-  2.43648116

pi0 <- exp(-alpha1/beta1)
prop <- (pi0^gamma1)/((1-pi0)^gamma1+pi0^gamma1)

