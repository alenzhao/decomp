###
## To test our algorithm on synthetic data previously generated.
##

# assume test.n, test.p, test.m, test.max.n are defined

source("sim.R");
source("infer.R");

#
test1 <- function(tn, sr,  n,st1,st2,max_n, max.delay, method=c("pcor","mi"), pruning=c("all","common")) {
  # test one random gene network of n genes, with p as the p.value threshold
  method <- match.arg(method);
  pruning <- match.arg(pruning);
  if(method == "pcor") {u.test.f <- test.corr; c.test.f <- test.p.corr;}
  if(method == "mi") {u.test.f <- mutual_info; c.test.f <- c_mutual_info;}
  #
  pn1 <- infer_grn1(sr,st1,max.delay, u.test.f); # max.delay is counted in time steps
  pn2 <- filter_pcor(sr,st2,pn1,max_n, c.test.f, pruning);
  # stage 1
  r1 <- compare_grn(tn$links, tn$delays, pn1);
  # stage 2
  r2 <- compare_grn(tn$links, tn$delays, pn2);
  list(stage1=r1, stage2=r2)
}

test.many <- function(indir, m, n,mp,md,n.points,e2,alpha, st1,st2,max_n, method=c("pcor","mi"), pruning=c("all","common")) {
  # syndata is a list containing m synthetic cases
  # call test1 m times on the cases
  # store the performance of each trial, and finally print the summary of the performances of all the trials
  method <- match.arg(method);
  pruning <- match.arg(pruning);
  # for stage 1
  d1.recalls <- rep(0,m);
  d1.precisions <- rep(0,m);
  d1.f.measures <- rep(0,m);
  e1.recalls <- rep(0,m);
  e1.precisions <- rep(0,m);
  e1.f.measures <- rep(0,m);
  l1.recalls <- rep(0,m);
  l1.precisions <- rep(0,m);
  l1.f.measures <- rep(0,m);
  # for stage 2
  d2.recalls <- rep(0,m);
  d2.precisions <- rep(0,m);
  d2.f.measures <- rep(0,m);
  e2.recalls <- rep(0,m);
  e2.precisions <- rep(0,m);
  e2.f.measures <- rep(0,m);
  l2.recalls <- rep(0,m);
  l2.precisions <- rep(0,m);
  l2.f.measures <- rep(0,m);
  #
  for(i in 1:m) {
    # the different nps uses the same file: that of nps200
    fname <- paste(indir,"/n",n,"mp",mp,"md",md,"nps200","e",e2,"a",alpha,"r",i,".txt", sep="");
    cat("Reading in ", fname, "\n");
    source(file=fname, local=TRUE, echo=FALSE);
    # cur.grn is the list containing the grn links and delays
    try({ x <- test1(cur.grn,cur.grn$sim.exp[1:n.points,], n,st1,st2,max_n, md, method, pruning);
          s1 <- x$stage1;
          s2 <- x$stage2;
    # for stage 1
          d1.recalls[i] <- s1$delays.recall;
          d1.precisions[i] <- s1$delays.precision;
          d1.f.measures[i] <- 2*(s1$delays.recall)*(s1$delays.precision)/(s1$delays.recall + s1$delays.precision);
          e1.recalls[i] <- s1$effects.recall;
          e1.precisions[i] <- s1$effects.precision;
          e1.f.measures[i] <- 2*(s1$effects.recall)*(s1$effects.precision)/(s1$effects.recall + s1$effects.precision);
          l1.recalls[i] <- s1$links.recall;
          l1.precisions[i] <- s1$links.precision;
          l1.f.measures[i] <- 2*(s1$links.recall)*(s1$links.precision)/(s1$links.recall + s1$links.precision);
    # for stage 2
          d2.recalls[i] <- s2$delays.recall;
          d2.precisions[i] <- s2$delays.precision;
          d2.f.measures[i] <- 2*(s2$delays.recall)*(s2$delays.precision)/(s2$delays.recall + s2$delays.precision);
          e2.recalls[i] <- s2$effects.recall;
          e2.precisions[i] <- s2$effects.precision;
          e2.f.measures[i] <- 2*(s2$effects.recall)*(s2$effects.precision)/(s2$effects.recall + s2$effects.precision);
          l2.recalls[i] <- s2$links.recall;
          l2.precisions[i] <- s2$links.precision;
          l2.f.measures[i] <- 2*(s2$links.recall)*(s2$links.precision)/(s2$links.recall + s2$links.precision);
        }); # to prevent premature termination of the script
    # print the one line result of the individual trial
    cat("\n");
    cat("====",
      n,n.points,e2,alpha,i,
      st1,st2, max_n,method,pruning,
      e1.recalls[i],e1.precisions[i],e1.f.measures[i],
      e2.recalls[i],e2.precisions[i],e2.f.measures[i],
      l1.recalls[i],l1.precisions[i],l1.f.measures[i],
      l2.recalls[i],l2.precisions[i],l2.f.measures[i],
      d1.recalls[i],d1.precisions[i],d1.f.measures[i],
      d2.recalls[i],d2.precisions[i],d2.f.measures[i], sep=",");
    cat("\n");
  }
  ###
  cat("\n");
  cat("####",
      n,n.points,e2,alpha,m,
      st1,st2, max_n,method,pruning,
      median(e1.recalls, na.rm=TRUE),median(e1.precisions, na.rm=TRUE),median(e1.f.measures, na.rm=TRUE),
      median(e2.recalls, na.rm=TRUE),median(e2.precisions, na.rm=TRUE),median(e2.f.measures, na.rm=TRUE),
      median(l1.recalls, na.rm=TRUE),median(l1.precisions, na.rm=TRUE),median(l1.f.measures, na.rm=TRUE),
      median(l2.recalls, na.rm=TRUE),median(l2.precisions, na.rm=TRUE),median(l2.f.measures, na.rm=TRUE),
      median(d1.recalls, na.rm=TRUE),median(d1.precisions, na.rm=TRUE),median(d1.f.measures, na.rm=TRUE),
      median(d2.recalls, na.rm=TRUE),median(d2.precisions, na.rm=TRUE),median(d2.f.measures, na.rm=TRUE),  sep=",");
  cat("\n");
}

test.all <- function(indir, n, nps, e2, alpha, pruning) {
  max.parents <- 4;
  max.delay <- 4;
  max.n <- 4;
  # 20 replicates
  # for pcor
  for(st in c(1, 1.30103, 2, 3, 4)) {
    # corresponding to p-value of 0.1 0.05 0.01 0.001 0.0001
    test.many(indir, 20, n,max.parents,max.delay, nps,e2,alpha, st,st, max.n, "pcor", pruning);
  }
  # for mi
  for(st in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4)) {
    test.many(indir, 20, n,max.parents,max.delay, nps,e2,alpha, st,st, max.n, "mi", pruning);
  }
}

# The calling is left to external script
#test.all(n10mp4, 10, 4, "all");
####

