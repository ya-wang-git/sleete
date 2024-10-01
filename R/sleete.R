#################################################################################################### .
#   Program/Function Name: sleete
#   Author: Zhiwei Zhang
#   Description: Super learning for efficient estimation of treatment effects with or without sample splitting
#   Change History:
#   Last Modified Date: 09/20/2024
#################################################################################################### .
#' @name sleete
#' @title sleete
#' @description { Super learning for efficient estimation of treatment effects using multiple algorithms}
#' @param data input dataset (a data frame)
#' @param resp name of response variable (observed event time for a survival outcome)
#' @param event observed event type (1 failure; 0 censoring) for a survival outcome
#' @param trt treatment indicator (numerical or character, with higher value representing experimental treatment)
#' @param stratcov name(s) of stratification variables for a stratified analysis
#' @param basecov.cont name(s) of continuous baseline covariates to adjust for
#' @param basecov.cat name(s) of categorical baseline covariates to adjust for
#' @param pi known probability of receiving treatment 1 (as opposed to 0)
#' @param bounds lower and upper bounds for the treatment effect of interest
#' @param method choice of initial estimator (log.HR, surv.diff, rmst.diff or WMW.cens for a survival outcome)
#' @param ... optional arguments required by the chosen method (e.g., tau for surv.diff and rmst.diff)
#' @param SL.library names of SuperLearner wrapper functions for the multiple algorithms to be combined
#' @param n.folds.cv number of folds in the cross-validation procedure within the super learner
#' @param sample.splitting logical indicator for the use of sample splitting
#' @param n.folds.cf number of folds to use in cross-fitting (i.e., sample splitting)
#' @return a (K+2)-by-2 matrix, where K=length(SL.library), consisting of the unadjusted and augmented estimates and their standard errors
#' @export
#' @examples
#' library(sleete)
#' data <- subset(colon, subset = ((etype == 2)&(rx != "Lev")))
#' data$trt <- as.numeric(data$rx == "Lev+5FU")
#' pi  <- 0.5
#' tau <- 5*365
#'
#' set.seed(12345)
#' fit.sleete <- sleete(data, resp = "time", event = "status", trt = "trt",
#'                      basecov.cont = c("age", "nodes", "differ", "extent"),
#'                      basecov.cat  = c("sex", "obstruct", "perfor", "adhere", "surg", "node4"),
#'                      pi = pi, method = "surv.diff", tau = tau)
#' round(fit.sleete, digits = 3)
#'
sleete <- function(data, resp, event=NULL, trt, stratcov=NULL, basecov.cont=NULL, basecov.cat=NULL,
                   pi=NULL, bounds=c(-Inf, Inf), method="log.HR", ...,
                   SL.library=c("SL.lm", "SL.gam", "SL.rpart", "SL.randomForest"),
                   n.folds.cv=5, sample.splitting=TRUE, n.folds.cf=5){

  method = eval(as.symbol(method))

  mData = data.preprocess(data=data, resp=resp, event=event, trt=trt,
                          basecov.cont=basecov.cont, basecov.cat=basecov.cat)
  mData = na.omit(mData)
  n = nrow(mData); K = length(SL.library)
  if (!is.null(event)) { a = mData[,3]; W = mData[,-(1:3)] } else { a = mData[,2]; W = mData[,-(1:2)] }
  if (is.null(pi)) pi = mean(a)
  est = function(dat) method$pt.est(dat,...)
  if (method$inf.fct.avail) {
    inf.fct = function(dat,I=1:nrow(dat),J=I) method$inf.fct(dat,I=I,J=J,pi=pi,...)
  } else {
    inf.fct = function(dat,I=1:nrow(dat),J=I) emp.inf.fct(est,dat,I=I,J=J)
  }
  est.0 = est(mData)
  inf = inf.fct(mData)
  se.0 = sd(inf)/sqrt(n)
  wt = (a-pi)^2; W.df = data.frame(W)
  if (sample.splitting) {
    inf.ss = numeric(n)
    b.hat.ss = matrix(0,n,K+1)
    cf = n.folds.cf
    fold = sample(1:cf, n, replace=TRUE)
    for (v in 1:cf) {
      val = (fold==v); trn = !val
      I = (1:n)[trn]; J = (1:n)[val]
      inf.ss[val] = inf.fct(mData,I=I,J=J)
      inf.trn = inf.fct(mData,I=I)
      rsp.trn = inf.trn/(a[trn]-pi)
      SL.rst.trn = SuperLearner(Y=rsp.trn, X=W.df[trn,], newX=W.df[val,], family=gaussian(),
                                SL.library=SL.library, cvControl=list(V=n.folds.cv), obsWeights=wt[trn])
      b.hat.ss[val,] = cbind(SL.rst.trn$library.predict, SL.rst.trn$SL.predict)
    }
    aug.mat.ss = (a-pi)*b.hat.ss
    est.aug = est.0-colMeans(aug.mat.ss)
    se.aug = sqrt(diag(var(inf.ss-aug.mat.ss))/n)
  } else {
    rsp = inf/(a-pi)
    SL.rst = SuperLearner(Y=rsp, X=W.df, newX=W.df, family=gaussian(), SL.library=SL.library,
                          cvControl=list(V=n.folds.cv), obsWeights=wt)
    b.hat = cbind(SL.rst$library.predict, SL.rst$SL.predict)
    aug.mat = (a-pi)*b.hat
    est.aug = est.0-colMeans(aug.mat)
    se.aug = sqrt(diag(var(inf-aug.mat))/n)
  }
  pe = trunc(c(est.0, est.aug), lo=bounds[1], hi=bounds[2])
  se = c(se.0, se.aug)
  rst = cbind(pe, se)
  colnames(rst) = c("Pt. Est.", "Std. Err.")
  SL.names = c(SL.library, "SL")
  rownames(rst) = c("Unadjusted", SL.names)
  rst
}

