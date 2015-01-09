#
# try to re-implement DD-lasso as reported in

# Ola ElBakry, M. Omair Ahmad, M. N. S. Swamy, "Inference of Gene
# Regulatory Networks with Variable Time Delay from Time-Series
# Microarray Data," IEEE/ACM Transactions on Computational Biology and
# Bioinformatics, vol. 10, no. 3, pp. 671-687, May-June, 2013
#

# this paper assumes replicates of the expression data (same number of
# genes, same number of time points), so we represent the expression
# data as follows:
#
# each column is gene, each row is time point.  the first n rows are
# for first time point, the second n rows are for the second time
# point, and so on.

library("lars")

dd.avg.reps <- function(data, n.rep) {
  # average the replicates
  ng <- ncol(data);
  nd <- nrow(data);
  n <- floor(nd/n.rep);
  r <- matrix(0, nrow=n, ncol=ng);
  for(j in 1:ng) {
    for(i in 1:n) {
      r[i,j] <- mean(data[(i-1)*n.rep + (1:n.rep),j]);
    }
  }
  #
  r
}

dd.find.delay <- function(data, max.delay) {
  # data is average of replicates, data[t,g] is expression of gene g at time t.
  # returns a ng by ng matrix, where the (i,j)th entry is the delay d
  # for which i --(d)--> j has the maximum abs. correlation.
  # first the replicates are averaged.
  ng <- ncol(data);
  n <- nrow(data);
  r <- matrix(0, nrow=ng, ncol=ng);
  for(i in 1:ng) {
    for(j in 1:ng) {
      tmp.m <- 0;
      for(d in 1:max.delay) {
        tmp.c <- abs(cor(data[1:(n-d),i], data[(1+d):n,j]));
        if(tmp.c > tmp.m) {r[i,j] <- d; tmp.m <- tmp.c;}
      }
    }
  }
  r
}

## test lars
# nr <- 50;
# nc <- 100;
# tmp.X <- matrix(rnorm(nr*nc), nrow=nr, ncol=nc);
# tmp.b <- c(8,-7,3,-9,6, rep(0, nc - 5));
# tmp.y <- tmp.X %*% tmp.b + rnorm(nr);
# 
# tmp2 <- cv.lars(x=tmp.X, y=tmp.y, plot.it=FALSE, mode="step");
# m.idx <- which.min(tmp2$cv); # the step with lowest CV value, resolution may not be too good, but that is what is convenient for lars
# tmp <- lars(x=tmp.X, y=tmp.y);
# tmp.p.b <- tmp$beta[m.idx,];


# idx <- which(tmp.p.b != 0);
# tmp.data <- as.data.frame(cbind(tmp.y,tmp.X[,idx]));
# tmp.lm <- lm(V1 ~ ., data=tmp.data);
# tmp.lm.s <- step(tmp.lm, direction="backward", trace=0, k=log(length(tmp.y))); # BIC
## gets the indices of the variables used in the model after backward elimination
## but the indices are with respect to tmp.data, which includes the y, so need adjustment
# s.idx <- which(names(tmp.data) %in% all.vars(formula(tmp.lm.s)));
# u.idx <- s.idx[-1]-1;

##tmp.lm.s2 <- step(tmp.lm, direction="backward", k=2, trace=5);
## seems to work for moderately many data points
####

cal.threshold.c <- function(N,p, alpha=10^(-N/p)) {
  # the threshold c given
  #   N: the total number of samples (including replicates)
  #   p: the number of genes.
  #   alpha: the control of sparseness, smaller -> more sparse.
  #       The default is the guideline given by the paper.
  #       Another default given by the paper is 0.0001
  A <- log(2*pi*N*alpha*alpha);
  B <- log(N)+2*log(p);
  target.f <- function(x) {x + log(B+x) + A};
  # search from around -B to around -A
  ## test.x <- -1000:10000/100;
  ## plot(x=test.x,y=target.f(test.x), type="l");
  # use the simple root finder, because target.f looks quite straight
  r <- uniroot(target.f,interval=c(-B+1e-6, -A+1e-6));
  r$root
}

mBIC2 <- function(RSS, df, N, p, threshold.c, sigma2) {
  # referring to [22] of the above paper,
  # mBIC2 is -2 Log(L(theta)) + df(log(N) + 2 log(p) + c) - 2 log(df!)
  # When sigma2 is known, -2 Log(L(theta)) can be replaced with
  #     RSS/sigma2 (+ other terms involving only known constants sigma2 and N).
  # When sigma2 is unknown, use it's ML estimate RSS/n, -2 Log(L(theta)) can be replaced with
  #     N Log(RSS) (+ other terms involving only known constants and N).
  # So compute only up to relevant constants.
  LogL <- if(is.null(sigma2)) (N*log(RSS)) else (RSS/sigma2);
  LogL + df*log(N) + df*(2*log(p)+threshold.c) - 2*log(factorial(df))
}

dd.gene.lasso.mBIC2 <- function(y, x, threshold.c, sigma2=NULL) {
  # return the vector of c(a,b) for y ~ a + Xb using LASSO
  # the intercept is fitted, but not returned here
  # for N >= p case in DD-lasso, somehow mBIC2 needs the sigma square of errors.
  # threshold.c is calculated by cal.threshold.c(),
  #   which depends only on N,p and alpha, so can be
  #   pre-calculated for a given dataset.
  # use mBIC2 to choose which step to use
  tmp.l <- lars(x=x, y=y);
  mBIC2.s <- mBIC2(tmp.l$RSS, tmp.l$df, length(y), ncol(x), threshold.c, sigma2);
  m.idx <- which.min(mBIC2.s);
  # the index of the non-zero entries of b can be found by which(b != 0).
  tmp.l$beta[m.idx,]
}

dd.gene.lasso.cv <- function(y, x) {
  # return the vector of c(a,b) for y ~ a + Xb using LASSO
  # the intercept is fitted, but not returned here.
  # used for N < p case in DD-lasso.
  # use cross validation to choose which step to use
  tmp.cv <- cv.lars(x=x,y=y, plot.it=FALSE, mode="step");
  # the step with lowest CV value, resolution may not be too good, but
  # that is what is convenient for lars
  m.idx <- which.min(tmp.cv$cv);
  tmp.l <- lars(x=x, y=y); # fit again, because cv.lars does not give the model
  # the index of the non-zero entries of b can be found by which(b != 0).
  tmp.l$beta[m.idx,]
}

dd.backward.elimination <- function(y, x, b) {
  # b is as returned by dd.gene.lasso for fitting y ~ a + Xb
  # use step() to do backward eliminationg by BIC.
  # b may contain zeros, which are ignored.
  # returned a list of:
  #   coef: an updated b where more entries may be zero, and
  #   idx: the indices of the non-zero entires in coef.
  #   intercept: the intercept of the final selected lm model
  idx <- which(b != 0);
  tmp.cx <- cbind(y,x[,idx]);
  colnames(tmp.cx) <- NULL; # to avoid messing up the names
  tmp.data <- as.data.frame(tmp.cx);
  tmp.lm <- lm(V1 ~ ., data=tmp.data);
  tmp.lm.s <- step(tmp.lm, direction="backward", trace=0, k=log(length(y))); # BIC
  ## gets the indices of the variables used in the model after backward elimination
  ## but the indices are with respect to tmp.data, which includes the y, so need adjustment
  s.idx <- which(names(tmp.data) %in% all.vars(formula(tmp.lm.s)));
  u.idx <- s.idx[-1]-1; # this is with respect to x[,idx];
  #
  out.b <- rep(0, length(b));
  out.idx <- idx[u.idx]; # with respect to x and b
  out.b[out.idx] <- coef(tmp.lm.s)[-1]; # the first is intercept
  list(coef=out.b, idx=out.idx, intercept=coef(tmp.lm.s)[1])
}

dd.delay.shift <- function(data, n.rep, idx, delays) {
  # data contains n.rep replicates (first n.rep rows are one time point,
  # the next n.rep rows are another time point and so on)/
  # idx and delays should have the same length.
  # get the columns indexed by idx, shifted properly by delays,
  # return the shifted submatrix.
  md <- max(delays);
  nr <- nrow(data) - md*n.rep;
  r.idx <- 1:nr;
  r <- matrix(0, nrow=nr, ncol=length(idx));
  for(i in 1:length(idx)) {
    r[,i] <- data[(md - delays[i])*n.rep + r.idx, idx[i]];
  }
  r
}

dd.grn.lasso <- function(data, n.rep, max.delay, sigma2=NULL, alpha=-1) {
  # to infer GRN from expression data in the format as described at the top
  # for mBIC2, need to know the sigma square of the errors.
  # return list of:
  #   links: links[i,j] is effect of gene i --> j
  #   delays: delays[i,j] is delay of gene i --> j if positive. 0 means no link
  
  # for each gene, all other genes (including itself) may be potential
  #   parent with delay, so p is ng
  ng <- ncol(data);
  nr <- nrow(data);
  # delays[i,j] is delay of i --> j
  delays <- dd.find.delay(dd.avg.reps(data, n.rep), max.delay);
  r.links <- matrix(0, nrow=ng, ncol=ng);
  r.delays <- matrix(0, nrow=ng, ncol=ng);

  for(i in 1:ng) {
    # do it gene by gene
    tmp.X <- dd.delay.shift(data, n.rep, 1:ng, delays[,i]);
    # N may be different for different genes, because they may have
    #  different delays, so different amount of truncation.
    N <- nrow(tmp.X);
    tmp.y <- data[(nr - N + 1):nr,i];
    #
    if(N >= ng) { # use mBIC2
      if(alpha < 0) {m.alpha <- 10^(-N/ng);} else {m.alpha <- alpha;}
      threshold.c <- cal.threshold.c(N, ng, m.alpha);

      tmp.b <- dd.gene.lasso.mBIC2(tmp.y, tmp.X, threshold.c, sigma2);
    } else { # use cross-validation
      tmp.b <- dd.gene.lasso.cv(tmp.y, tmp.X);
    }
    #
    p.coef <- dd.backward.elimination(tmp.y, tmp.X, tmp.b);
    r.links[,i] <- p.coef$coef;
    r.delays[p.coef$idx, i] <- delays[p.coef$idx, i];
  }
  #
  list(links=r.links, delays=r.delays)
}
