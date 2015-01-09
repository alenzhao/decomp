###
## To test DD-lasso on synthetic data previously generated.
##

source("dd_lasso.R");

#
print.mat.grn <- function(links,delays) {
  ## all use 0-based indices
  ##   links: links[i,j] is effect of gene i --> j
  ##   delays: delays[i,j] is delay of gene i --> j if positive. 0 means no link
  n <- nrow(links); # should be same as ncol(links)=nrow(delays)=ncol(delays)
  for(j in 1:n) { # to
    for(i in which(delays[,j]!=0)) { # from
      cat("To:",j-1, "From:",i-1, "Delay:",delays[i,j], "Effect:",links[i,j], sep="\t");
      cat("\n");
    }
  }
}

# default parameters
test.maxdelay <- 4;

get.expressions <- function(ns) {
  # ns is possibly a single file name, or a vector of file names
  # return a list of a matrix in a format ready for DD-lasso, and the number of replicates
  n <- length(ns);
  if(n == 1) {
    tmp <- read.table(file=ns[1], header=FALSE);
    list(d=data.matrix(tmp), reps=1)
  } else {
    lx <- list();
    L <- rep(0,n);
    nc <- 0;
    for(i in 1:n) {
      lx[[i]] <- read.table(file=ns[i], header=FALSE);
      L[i] <- nrow(lx[[i]]);
      nc <- ncol(lx[[i]]);
    }
    # truncate all the replicates to the shortest one
    mL <- min(L);
    r <- matrix(0, nrow=n*mL, ncol=nc);
    for(i in 1:n) {
      tmp <- data.matrix(lx[[i]]);
      r[seq(from=i, to=(mL-1)*n + i, by=n),] <- tmp[1:mL,]
    }
    #
    list(d=r, reps=n)
  }
}

test.one.case <- function(exp.filename, out.grn.filename) {
  #sr.t <- read.table(file=exp.filename, header=FALSE);
  #sr <- data.matrix(sr.t);
  sr <- get.expressions(exp.filename);
  pn1 <- dd.grn.lasso(sr$d, sr$reps, test.maxdelay);
  # pn1$links and pn1$delays contain the predicted links coefficients and delays

  sink(out.grn.filename);
  print.mat.grn(pn1$links, pn1$delays);
  sink();
}

# The calling is left to external script
####


