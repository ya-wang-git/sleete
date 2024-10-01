#################################################################################################### .
#   Program/Function Name: utility functions
#   Author: Zhiwei Zhang
#   Description: utility functions for super learner and machine learning algorithms
#   Change History:
#   Last Modified Date: 09/20/2024
####################################################################################################
# preliminary functions

#-----------------------------#
# truncation
#-----------------------------#
trunc = function(u, lo=-Inf, hi=Inf) pmax(lo,pmin(hi,u))

#-----------------------------#
# logit = log-odds
#-----------------------------#
logit = function(p) log(p/(1-p))

#-----------------------------#
# Agresti definition of h
#-----------------------------#
h0 = function(y1, y0) as.numeric(y1>y0)-as.numeric(y1<y0)

#-----------------------------#
# Mann-Whitney definition of h
#-----------------------------#
h1 = function(y1, y0) as.numeric(y1>y0)+0.5*as.numeric(y1==y0)

#------------------------------------------------------------#
# Kaplan-Meier estimate of a restricted (by tau) distribution
#------------------------------------------------------------#
km.tau = function(x, delta, tau=Inf) {
  dat = Surv(x, delta)
  km.est = survfit(dat~1)
  t.all = c(km.est$time, Inf)
  s.all = c(1, km.est$surv, 0)
  p.all = s.all[-length(s.all)]-s.all[-1]
  t.is.small = (t.all<tau)
  t.tau = c(t.all[t.is.small], tau)
  p.tau = p.all[t.is.small]
  pr.lt.tau = sum(p.tau)
  p.tau = c(p.tau, 1-pr.lt.tau)
  list(t=t.tau, p=p.tau)
}

#-------------------------------------------------------------#
# influence function of the Kaplan-Meier estimator at time tau
# estimated from I, then applied to J
#-------------------------------------------------------------#
km.inf.fct = function(x, delta, tau=Inf, I=1:length(x), J=I) {
  n = length(x); n.I = length(I); n.J = length(J)
  km.rst = km.tau(x[I], delta[I])
  S.hat = 1-sum(km.rst$p[km.rst$t<=tau])
  y.bar = rep(NA, n)
  for (i in 1:n) y.bar[i] = mean(x[I]>=x[i])
  inf.1 = delta[J]*(x[J]<=tau)/y.bar[J]
  inf.2 = rep(NA, n.J)
  for (i in 1:n.J) inf.2[i] = sum(delta[I]*(x[I]<=min(tau,x[J[i]]))/(y.bar[I]^2))/n.I
  inf = S.hat*(inf.2-inf.1)
  ifelse(is.na(inf), 0, inf)
}

#-----------------------------------------------------------#
# influence function of the KM estimator of RMST at time tau
# estimated from I, then applied to J
#-----------------------------------------------------------#
rmst.inf.fct = function(x, delta, tau=Inf, I=1:length(x), J=I) {
  n = length(x); n.I = length(I); n.J = length(J)
  km.rst = km.tau(x[I], delta[I], tau=tau)
  S.int = rep(0, n)
  for (i in 1:n) {
    in.window = (km.rst$t>=x[i])
    if (sum(in.window)>0) S.int[i] = sum((km.rst$t[in.window]-x[i])*km.rst$p[in.window])
  }
  y.bar = rep(NA, n)
  for (i in 1:n) y.bar[i] = mean(x[I]>=x[i])
  inf.1 = delta[J]*(x[J]<=tau)*S.int[J]/y.bar[J]
  inf.2 = rep(NA, n.J)
  for (i in 1:n.J) inf.2[i] = sum(delta[I]*(x[I]<=min(tau,x[J[i]]))*S.int[I]/(y.bar[I]^2))/n.I
  inf = inf.2-inf.1
  ifelse(is.na(inf), 0, inf)
}


# methods

#-------------------------------------------------------#
# method for difference between two means or proportions
#-------------------------------------------------------#
# point estimate
pt.est.mean.diff = function(mData) {
  y = mData[,1]; a = mData[,2]
  mean(y[a>0.5])-mean(y[a<0.5])
}
# influence function estimated from I, applied to J
inf.fct.mean.diff = function(mData, I=1:nrow(mData), J=I, pi=NULL) {
  y = mData[,1]; a = mData[,2]
  if (is.null(pi)) pi = mean(a[I])
  mu1 = mean(y[I][a[I]>0.5]); mu0 = mean(y[I][a[I]<0.5])
  (a[J]*(y[J]-mu1)/pi)-((1-a[J])*(y[J]-mu0)/(1-pi))
}
# the method
mean.diff = list(pt.est=pt.est.mean.diff, inf.fct.avail=TRUE, inf.fct=inf.fct.mean.diff)

#-------------------------------------------------#
# method for log-ratio of two means or proportions
#-------------------------------------------------#
# point estimate
pt.est.log.ratio = function(mData) {
  y = mData[,1]; a = mData[,2]
  log(mean(y[a>0.5])/mean(y[a<0.5]))
}

# influence function estimated from I, applied to J
inf.fct.log.ratio = function(mData, I=1:nrow(mData), J=I, pi=NULL) {
  y = mData[,1]; a = mData[,2]
  if (is.null(pi)) pi = mean(a[I])
  mu1 = mean(y[I][a[I]>0.5]); mu0 = mean(y[I][a[I]<0.5])
  (a[J]*(y[J]-mu1)/(pi*mu1))-((1-a[J])*(y[J]-mu0)/((1-pi)*mu0))
}

# the method
log.ratio = list(pt.est=pt.est.log.ratio, inf.fct.avail=TRUE, inf.fct=inf.fct.log.ratio)

#---------------------------------------------#
# method for log-odds-ratio of two proportions
#---------------------------------------------#
# point estimate
pt.est.log.OR = function(mData) {
  y = mData[,1]; a = mData[,2]
  logit(mean(y[a>0.5]))-logit(mean(y[a<0.5]))
}

# influence function estimated from I, applied to J
inf.fct.log.OR = function(mData, I=1:nrow(mData), J=I, pi=NULL) {
  y = mData[,1]; a = mData[,2]
  if (is.null(pi)) pi = mean(a[I])
  p1 = mean(y[I][a[I]>0.5]); p0 = mean(y[I][a[I]<0.5])
  (a[J]*(y[J]-p1)/(pi*p1*(1-p1)))-((1-a[J])*(y[J]-p0)/((1-pi)*p0*(1-p0)))
}

# the method
log.OR = list(pt.est=pt.est.log.OR, inf.fct.avail=TRUE, inf.fct=inf.fct.log.OR)

#----------------------------------------------------------------------#
# method for Mann-Whitney effect based on an arbitrary h (default = h0)
#----------------------------------------------------------------------#
# point estimate
pt.est.MW = function(mData, h=h0) {
  y = mData[,1]; a = mData[,2]
  mean(outer(y[a>0.5], y[a<0.5], FUN=h))
}

# influence function estimated from I, applied to J
inf.fct.MW = function(mData, I=1:nrow(mData), J=I, pi=NULL, h=h0) {
  y = mData[,1]; a = mData[,2]
  if (is.null(pi)) pi = mean(a[I])
  theta = pt.est.wmw(y[I],a[I],h=h)
  m = length(J); inf = numeric(m)
  for (k in 1:m) {
    if (a[J[k]]>0.5) {
      inf[k] = (mean(h(y[J[k]],y[I]))-theta)/pi
    } else {
      inf[k] = (mean(h(y[I],y[J[k]]))-theta)/(1-pi)
    }
  }
  inf
}

# the method
MW = list(pt.est=pt.est.MW, inf.fct.avail=TRUE, inf.fct=inf.fct.MW)

#----------------------------------------------------------------------#
# method for Mann-Whitney effect based on an arbitrary h (default = h0)
# estimated from randomly right-censored data restricted by tau
#----------------------------------------------------------------------#
# point estimate
pt.est.MW.cens = function(mData, tau=Inf, h=h0) {
  x = mData[,1]; delta = mData[,2]; a = mData[,3]
  a.eq.1 = (a>0.5)
  km.rst.1 = km.tau(x[a.eq.1], delta[a.eq.1], tau=tau)
  t1 = km.rst.1$t; p1 = km.rst.1$p
  km.rst.0 = km.tau(x[!a.eq.1], delta[!a.eq.1], tau=tau)
  t0 = km.rst.0$t; p0 = km.rst.0$p
  H = outer(t1, t0, FUN=h)
  as.vector(t(p1)%*%H%*%p0)
}

# the method
MW.cens = list(pt.est=pt.est.MW.cens, inf.fct.avail=FALSE)

#----------------------------------------------------------------------------#
# method for difference between two survival probabilities at time tau
# estimated from randomly right-censored data using the Kaplan-Meier approach
#----------------------------------------------------------------------------#
# point estimate
pt.est.surv.diff = function(mData, tau=Inf) {
  x = mData[,1]; delta = mData[,2]; a = mData[,3]
  a.eq.1 = (a>0.5)
  km.rst.1 = km.tau(x[a.eq.1], delta[a.eq.1])
  s1 = 1-sum(km.rst.1$p[km.rst.1$t<=tau])
  km.rst.0 = km.tau(x[!a.eq.1], delta[!a.eq.1])
  s0 = 1-sum(km.rst.0$p[km.rst.0$t<=tau])
  s1-s0
}

# influence function
inf.fct.surv.diff = function(mData, tau=Inf, I=1:nrow(mData), J=I, pi=NULL) {
  n = nrow(mData)
  x = mData[,1]; delta = mData[,2]; a = mData[,3]
  if (is.null(pi)) pi = mean(a)
  I.0 = I[a[I]<0.5]; I.1 = I[a[I]>0.5]
  J.0 = J[a[J]<0.5]; J.1 = J[a[J]>0.5]
  inf.all = rep(0, n)
  inf.all[J.0] = -km.inf.fct(x, delta, tau=tau, I=I.0, J=J.0)/(1-pi)
  inf.all[J.1] = km.inf.fct(x, delta, tau=tau, I=I.1, J=J.1)/pi
  inf = inf.all[J]
  ifelse(is.na(inf), 0, inf)
}

# the method
surv.diff = list(pt.est=pt.est.surv.diff, inf.fct.avail=TRUE, inf.fct=inf.fct.surv.diff)

#----------------------------------------------------------------------------#
# method for difference between two mean restricted (by tau) survival times
# estimated from randomly right-censored data using the Kaplan-Meier approach
#----------------------------------------------------------------------------#
# point estimate
pt.est.RMST.diff = function(mData, tau=Inf) {
  x = mData[,1]; delta = mData[,2]; a = mData[,3]
  a.eq.1 = (a>0.5)
  km.rst.1 = km.tau(x[a.eq.1], delta[a.eq.1], tau=tau)
  mu1 = sum(km.rst.1$p*km.rst.1$t)
  km.rst.0 = km.tau(x[!a.eq.1], delta[!a.eq.1], tau=tau)
  mu0 = sum(km.rst.0$p*km.rst.0$t)
  mu1-mu0
}

# influence function
inf.fct.RMST.diff = function(mData, tau=Inf, I=1:nrow(mData), J=I, pi=NULL) {
  n = nrow(mData)
  x = mData[,1]; delta = mData[,2]; a = mData[,3]
  if (is.null(pi)) pi = mean(a)
  I.0 = I[a[I]<0.5]; I.1 = I[a[I]>0.5]
  J.0 = J[a[J]<0.5]; J.1 = J[a[J]>0.5]
  inf.all = rep(0, n)
  inf.all[J.0] = -rmst.inf.fct(x, delta, tau=tau, I=I.0, J=J.0)/(1-pi)
  inf.all[J.1] = rmst.inf.fct(x, delta, tau=tau, I=I.1, J=J.1)/pi
  inf = inf.all[J]
  ifelse(is.na(inf), 0, inf)
}

# the method
RMST.diff = list(pt.est=pt.est.RMST.diff, inf.fct.avail=TRUE, inf.fct=inf.fct.RMST.diff)
MRST.diff = RMST.diff

#-----------------------------------------#
# method for log-hazard-ratio in Cox model
#-----------------------------------------#
# point estimate
pt.est.log.HR = function(mData) {
  x = mData[,1]; delta = mData[,2]; a = mData[,3]
  dat = Surv(x, delta)
  mod = coxph(dat~a)
  coef(mod)
}

# influence function
inf.fct.log.HR = function(mData, I=1:nrow(mData), J=I, pi=NULL) {
  n = nrow(mData)
  x = mData[,1]; delta = mData[,2]; a = mData[,3]
  dat = Surv(x, delta)
  mod = coxph(dat~a, subset=I, robust=FALSE)
  b = coef(mod); v = length(I)*as.vector(vcov(mod))
  e = exp(b*a); ea = e*a
  S0 = S1 = rep(NA, n)
  for (i in 1:n) {
    S0[i] = mean((x[I]>=x[i])*e[I])
    S1[i] = mean((x[I]>=x[i])*ea[I])
  }
  r = S1/S0
  inf.1 = delta[J]*(a[J]-r[J]); inf.2 = rep(NA, n)
  for (j in J) inf.2[j] = mean(delta[I]*(x[j]>=x[I])*e[j]*(a[j]-r[I])/S0[I])
  inf = inf.1-inf.2[J]
  v*ifelse(is.na(inf), 0, inf)
}

# the method
log.HR = list(pt.est=pt.est.log.HR, inf.fct.avail=TRUE, inf.fct=inf.fct.log.HR)

#-------------------------------------------#
# method for log-HR in a stratified PH model
#-------------------------------------------#
# point estimate
pt.est.log.HR.strat = function(mData) {
  x = mData[,1]; delta = mData[,2]; a = mData[,3]; s = mData[,4]
  dat = Surv(x, delta)
  mod = coxph(dat~a+strata(s))
  coef(mod)
}

# the method
log.HR.strat = list(pt.est=pt.est.log.HR.strat, inf.fct.avail=FALSE)

#--------------------------------------------------------------------#
# empirical influence function for a generic point estimator
# estimated from subjects in set I, then applied to subjects in set J
#--------------------------------------------------------------------#
emp.inf.fct = function(est, mData, I=1:nrow(mData), J=I) {
  est0 = est(mData[I,])
  m = length(I); n = length(J); jack.est = numeric(n)
  for (j in 1:n) {
    Ij = c(I, J[j])
    jack.est[j] = est(mData[Ij,])
  }
  (m+1)*(jack.est-est0)
}

