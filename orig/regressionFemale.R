#!/usr/bin/env Rscript
options("width"=200)

library(gplots)
library(MASS)
library(MatchIt)
library(dummies)


data=read.csv("/home/acquireassets/Desktop/newdata.csv")



data1=data[data$stateicp%in%c("California","New York","Massachusetts","Maryland","Michigan","Ohio",
"Mississippi","Oklahoma","Texas","Florida","North Carolina","Virginia","Georgia","Tennessee",
"Missouri","Wisconsin","Vermont","Washington","Alabama","Arizona","Illinois",
"Idaho","Kansas","Mississippi","Nebraska","Oklahoma","South Carolina","South Dakota","Utah","Wyoming"),]



data1=data1[data1$year==2012|data1$year==2016,]



data1$time=as.numeric(data1$year!=2012)


data1$group=as.numeric(data1$stateicp%in%c("Alabama","Florida","Georgia","Idaho","Kansas","Mississippi","Missouri","Nebraska","North Carolina","Oklahoma","South Carolina","South Dakota","Tennessee","Texas","Utah","Virginia","Wisconsin","Wyoming"))

data1$group1=abs(1-data1$group)



data2=data1[data1$incearn>1,]
data2$outcome=(data2$incearn)

data2=data2[data2$sex=="Female",]

 model1=rlm(outcome~costelec+costgas+costwatr+costfuel+valueh+perwt+pernum+
as.numeric(age)+as.factor(marst)+
as.factor(race)+as.factor(hcovany)+as.factor(empstat)
+as.factor(classwkrd)+as.factor(time)+as.factor(group1)+time:group1,data=data2)




good=summary(model1)
dd = data.frame(good$coefficients) 
dd$p.value =  2*pt(abs(dd$t.value),good$df[2], lower.tail=FALSE)      
 dd
################################################
# Code for Economic Theory B
#
# economictheoryb.com
#
# August 2016
# a.d.


summary.lm <- function (object, correlation = FALSE, 
                        symbolic.cor = FALSE, robust=FALSE,
                        cluster=c(NULL,NULL),...) {
  # add extension for robust standard errors
  if(robust==TRUE){ 
    # save variable that are necessary to calcualte robust sd
    X <- model.matrix(object)
    u2 <- residuals(object)^2
    XDX <- 0
    
    ## One needs to calculate X'DX. But due to the fact that
    ## D is huge (NxN), it is better to do it with a cycle.
    for(i in 1:nrow(X)) {
      XDX <- XDX + u2[i]*X[i,]%*%t(X[i,])
    }
    
    # inverse(X'X)
    XX1 <- solve(t(X)%*%X,tol = 1e-100)
    
    # Sandwich Variance calculation (Bread x meat x Bread)
    varcovar <- XX1 %*% XDX %*% XX1
    
    # adjust degrees of freedom 
    dfc_r <- sqrt(nrow(X))/sqrt(nrow(X)-ncol(X))
    
    # Standard errors of the coefficient estimates are the
    # square roots of the diagonal elements
    rstdh <- dfc_r*sqrt(diag(varcovar))
  }
  # add extension for clustered standard errors
  if(!is.null(cluster)&robust==T){warning("Robust standard errors are calculated. Set robust=F to calculate clustered standard errors.")}
  if(!is.null(cluster)&robust==F){
    if(""%in%cluster){stop("No variable for clustering provided.")}
    if(length(cluster)>2){stop("The function only allows max. 2 clusters. You provided more.")}
    n_coef <- all.vars(object$call$formula)
    if(length(cluster)==1){
      dat <- na.omit(get(paste(object$call$data))[,c(n_coef,cluster)])
      if(nrow(dat)<nrow(object$model)){stop("Not all observation have a cluster.")}
      cluster_vector <- dat[,cluster]
      require(sandwich, quietly = TRUE)
      M <- res_length <- length(unique(cluster_vector))
      N <- length(cluster_vector)
      K <- object$rank
      dfc <- (M/(M-1))*((N-1)/(N-K))
      uj  <- na.omit(apply(estfun(object),2, function(x) tapply(x, cluster_vector, sum)));
      varcovar <- dfc*sandwich(object, meat=crossprod(uj)/N)
      rstdh <- sqrt(diag(varcovar))
    } 
    if(length(cluster)==2){
      dat_1 <- na.omit(get(paste(object$call$data))[,c(n_coef,cluster[1])])
      if(nrow(dat_1)<nrow(object$model)){stop("Not all observation have a cluster.")}
      dat_2 <- na.omit(get(paste(object$call$data))[,c(n_coef,cluster[2])])
      if(nrow(dat_2)<nrow(object$model)){stop("Not all observation have a cluster.")}
      
      dat <- na.omit(get(paste(object$call$data))[,c(n_coef,cluster)])
      library(sandwich,quietly = TRUE)
      cluster1 <- dat[,cluster[1]]
      cluster2 <- dat[,cluster[2]]
      cluster12 = paste(cluster1,cluster2, sep="")
      M1  <- length(unique(cluster1))
      M2  <- length(unique(cluster2))   
      M12 <- res_length <-length(unique(cluster12))
      N   <- length(cluster1)          
      K   <- object$rank             
      dfc1  <- (M1/(M1-1))*((N-1)/(N-K))  
      dfc2  <- (M2/(M2-1))*((N-1)/(N-K))  
      dfc12 <- (M12/(M12-1))*((N-1)/(N-K))  
      u1j   <- apply(estfun(object), 2, function(x) tapply(x, cluster1,  sum)) 
      u2j   <- apply(estfun(object), 2, function(x) tapply(x, cluster2,  sum)) 
      u12j  <- apply(estfun(object), 2, function(x) tapply(x, cluster12, sum)) 
      vc1   <-  dfc1*sandwich(object, meat=crossprod(u1j)/N )
      vc2   <-  dfc2*sandwich(object, meat=crossprod(u2j)/N )
      vc12  <- dfc12*sandwich(object, meat=crossprod(u12j)/N)
      varcovar <- vc1 + vc2 - vc12
      rstdh <- sqrt(diag(varcovar))
    } 
    
  }
  z <- object
  p <- z$rank
  rdf <- 14530-34#n - p #z$df.residual
  if (p == 0) {
    r <- z$residuals
    n <- length(r)
    w <- z$weights
    if (is.null(w)) {
      rss <- sum(r^2)
    }
    else {
      rss <- sum(w * r^2)
      r <- sqrt(w) * r
    }
    resvar <- rss/rdf
    ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
    class(ans) <- "summary.lm"
    ans$aliased <- is.na(coef(object))
    ans$residuals <- r
    ans$df <- c(0L, n, length(ans$aliased))
    ans$coefficients <- matrix(NA, 0L, 4L)
    dimnames(ans$coefficients) <- list(NULL, c("Estimate", 
                                               "Std. Error", "t value", "Pr(>|t|)"))
    ans$sigma <- sqrt(resvar)
    ans$r.squared <- ans$adj.r.squared <- 0
    return(ans)
  }
  if (is.null(z$terms)) 
    stop("invalid 'lm' object:  no 'terms' component")
  if (!inherits(object, "lm")) 
    warning("calling summary.lm(<fake-lm-object>) ...")
  Qr <- stats:::qr.lm(object)
  n <- NROW(Qr$qr)
  print(n)
  print(p)
  print(z$df.residual)
  if (is.na(z$df.residual) || n - p != z$df.residual) 
    warning("residual degrees of freedom in object suggest this is not an \"lm\" fit")
  r <- z$residuals
  f <- z$fitted.values
  w <- z$weights
  if (is.null(w)) {
    mss <- if (attr(z$terms, "intercept")) 
      sum((f - mean(f))^2)
    else sum(f^2)
    rss <- sum(r^2)
  }
  else {
    mss <- if (attr(z$terms, "intercept")) {
      m <- sum(w * f/sum(w))
      sum(w * (f - m)^2)
    }
    else sum(w * f^2)
    rss <- sum(w * r^2)
    r <- sqrt(w) * r
  }
  resvar <- rss/rdf
  if (is.finite(resvar) && resvar < (mean(f)^2 + var(f)) * 
      1e-30) 
    warning("essentially perfect fit: summary may be unreliable")
  p1 <- 1L:p
  R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  se <- sqrt(diag(R) * resvar)
  
  if(robust==T){se <- rstdh}
  if(!is.null(cluster)&robust==F){se <- rstdh}
  est <- z$coefficients[Qr$pivot[p1]]
  tval <- est/se
  ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
  ans$residuals <- r
  pval <- 2 * pt(abs(tval), 
                 rdf, lower.tail = FALSE)
  ans$coefficients <- cbind(est, se, tval, pval)
  dimnames(ans$coefficients) <- list(names(z$coefficients)[Qr$pivot[p1]], 
                                     c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  ans$aliased <- is.na(coef(object))
  ans$sigma <- sqrt(resvar)
  ans$df <- c(p, rdf, NCOL(Qr$qr))
  if (p != attr(z$terms, "intercept")) {
    df.int <- if (attr(z$terms, "intercept")) 
      1L
    else 0L
    ans$r.squared <- mss/(mss + rss)
    ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - 
                                                       df.int)/rdf)
    ans$fstatistic <- c(value = (mss/(p - df.int))/resvar, 
                        numdf = p - df.int, dendf = rdf)
    if(robust==T|(!is.null(cluster))){
      if(!is.null(cluster)){rdf <- res_length -1}
      pos_coef <- match(names(z$coefficients)[-match("(Intercept)",
                                                     names(z$coefficients))],
                        names(z$coefficients))
      
      P_m <- matrix(z$coefficients[pos_coef])
      
      R_m <- diag(1, 
                  length(pos_coef), 
                  length(pos_coef))
      
      ans$fstatistic <- c(value = t(R_m%*%P_m)%*%
                            (solve(varcovar[pos_coef,pos_coef],tol = 1e-100))%*%
                            (R_m%*%P_m)/(p - df.int), 
                          numdf = p - df.int, dendf = rdf)
      
    }
    
  }
  else ans$r.squared <- ans$adj.r.squared <- 0
  ans$cov.unscaled <- R
  dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1, 
                                                             1)]
  if (correlation) {
    ans$correlation <- (R * resvar)/outer(se, se)
    dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
    ans$symbolic.cor <- symbolic.cor
  }
  if (!is.null(z$na.action)) 
    ans$na.action <- z$na.action
  class(ans) <- "summary.lm"
  ans
}
summary.lm(model1, cluster=c("stateicp"))
data1=data[data$stateicp%in%c("California","New York","Massachusetts","Maryland","Michigan","Ohio",
"Mississippi","Oklahoma","Texas","Florida","North Carolina","Virginia","Georgia","Tennessee",
"Missouri","Wisconsin","Vermont","Washington","Alabama","Arizona","Illinois",
"Idaho","Kansas","Mississippi","Nebraska","Oklahoma","South Carolina","South Dakota","Utah","Wyoming"),]


data1=data1[data1$year>=2008,]

data1=data1[data1$year<=2016,]

data1$group=as.numeric(data1$stateicp%in%c("Maine","Alabama","Florida","Georgia","Idaho","Kansas","Mississippi","Missouri","Nebraska","North Carolina","Oklahoma","South Carolina","South Dakota","Tennessee","Texas","Utah","Virginia","Wisconsin","Wyoming"))

data1$group1=abs(1-data1$group)


data2=data1[data1$incearn>1,]
data2$outcome=(data2$incearn)

#model1=lm(outcome~costelec+costgas+costwatr+costfuel+valueh+perwt+pernum+as.factor(sex)+age+as.factor(marst)+
#as.factor(race)+as.factor(hcovany)+as.factor(educ)+as.factor(empstat)
#+as.factor(classwkrd)+incinvst+as.factor(time)+as.factor(group1)+time:group1,data=data2)

data2=data2[data2$sex=="Female",]

data3=data2[data2$group1==0,]

data4=data2[data2$group1==1,]
par(mfrow=c(1,1))
#plotmeans(outcome~year,data=data3)
#plotmeans(outcome~year,data=data4)

comparison=aggregate(incearn ~ year, data3, mean)
comparison1=aggregate(incearn ~ year, data4, mean)

png('regressionFemale.png')
#plot(comparison1[,1],comparison1[,2],type="l",ylim=c(184000,250000),col="blue",xlab="year",ylab="Mean income")
plot(comparison1[,1],comparison1[,2],type="l",col="blue",xlab="year",ylab="Mean income")

lines(comparison[,1],comparison[,2], col=c("darkgreen"))


aggregate(incearn~ group1, data2, mean)

dev.off()



