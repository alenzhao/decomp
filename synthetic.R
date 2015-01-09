
## modified on 16 Jun, 2014
## written on 14 May, 2014
###
### To generate synthetic network and expression data, for testing on Minimum Description Length (MDL) to learn its structure.

#
# to randomly generate (synthetic) gene network in the form of matrices,
# where each link a_ij has a number representing the
# effect of gene i on gene j: +ve for activation, -ve for inhibition.
# Associated with each a_ij =/= 0 is t_ij > 0, which represents the time delay
# of the effect of gene i on gene j.

var.transition.mat <- function(links, delays) {
  # the matrix A of the Vector Autoregressive Model, if the time series
  # x(t) = a_1*x(t-1) + a_2*x(t-2) + ... + a_d*x(t-d)
  # written as X(t) = A X(t-1), where X(t) = [x'(t) x'(t-1) ... x'(t-d)]',
  #         |a_1 a_2 ... a_d|
  #         | I   0  ...  0 |
  # and A = | 0   I  ...  0 |
  #         |        ...    |
  #         |          I  0 |

  n <- nrow(links);
  md <- max(delays);
  m <- matrix(0, nrow=md*n, ncol=md*n);

  # the a's
  for(i in 1:n) {
    for(j in 1:n) {
      if(delays[i,j] > 0) {
        m[j, i + (delays[i,j]-1)*n] <- links[i,j]; # links[i,j] is effect of gene i on gene j
      }
    }
  }
  # the I's
  for(i in 1:((md-1)*n)) {
    m[n+i,i] <- 1;
  }
  #
  m
}

max.abs.eigenvalue <- function(links, delays) {
  x <- eigen(var.transition.mat(links,delays));
  max(abs(x$values))
}

gen.grn <- function(n,is_acyclic,is_self, max.parents) {
  # n is the number of gene
  # is_acyclic is true iff only acyclic network is to be generated
  # is_self is true iff when is_acyclic is false, to include self loops
  # p_link is the probability of a link
  # Returns a list of two matrices of n by n, one is the links, the other
  # is the delays.
  L <- n*n;
  r <- rep(0,L);
  if(is_acyclic) {
    mp <- sample(0:max.parents,n, replace=TRUE); # up to max.parents
  } else { # more likely to have cycle
    mp <- sample(1:max.parents,n, replace=TRUE); # up to max.parents
  }

  # limits the number of non-zero entries in each column
  for(j in 1:n) {
    if(mp[j] > 0) {
      ps <- sample(1:n, mp[j]);
      for(i in 1:length(ps)) {
        si <- (ps[i]-1)*n + j;
        r[si] <- 1;
      }
    }
  }
  #
  if(is_acyclic) {
    # make it upper triangular
    for(i in 2:n) {
      for(j in 1:(i-1)) {
        r[(i-1)*n + j] <- 0;
      }
    }
  }
  if(!is_self) {
    # make the diagonal zero
    for(i in 1:n) {r[(i-1)*n + i] <- 0;}
  }
  delay <- r*runif(L,0,1); # simulate less delays
  r <- r*(1-2*rbinom(L,1,0.5))*runif(L,0.5,1.5);
  # returns r and delay
  list(links=matrix(data=r, nrow=n, ncol=n, byrow=TRUE),
       delays=matrix(data=delay, nrow=n, ncol=n, byrow=TRUE))
}

gen.stable.grn <- function(n,is_acyclic,is_self, max.parents, max.delay) {
  g <- gen.grn(n, is_acyclic, is_self, max.parents);
  g$delays <- ceiling(g$delays * max.delay);
  # scale it so that the eigen value of VAR model has modulus
  # will not be too large
  e <- max.abs.eigenvalue(g$links, g$delays);
  if(e >= 1) {
    g$links <- g$links * 0.5 / e;
  }
  g
}

permute_grn <- function(links, delays) {
  # to randomly permute the genes so that the position do not give any advantage
  x <- sample(1:nrow(links));
  list(links=links[x,x], delays=delays[x,x])
}

### example of The synthetic network: 
# a_ij is effect of gene i on gene j.
# mdl.links <- matrix(c(0    , 0.7 , 0    , 0.31, 0.8 ,
#                       0    , 0   , -0.85, 0   , 0   ,
#                       -0.55, 0.64, 0    , 0   , 0   ,
#                       0    , 0   , 0    , 0   , -0.6,
#                       0    , 0   , 0.4  , 0   , 0   ),
#                     nrow=5, ncol=5, byrow=TRUE);
# mdl.delay <- matrix(c(0, 1, 0, 4, 2,
#                       0, 0, 1, 0, 0,
#                       2, 1, 0, 0, 0,
#                       0, 0, 0, 0, 1,
#                       0, 0, 1, 0, 0),
#                     nrow=5, ncol=5, byrow=TRUE);
# 
#
sim_grn <- function(links, delays, N, init.mean=1, init.sigma2=1, err.sigma2=1, alpha=1) {
  # links is a matrix encoding the links of a gene network (direction and magnitude)
  # delays contains the time delay (in steps) for each link
  # Returns a matrix a matrix of N by n, where n is the number of genes in the gene network.
  # The expression should be interpreted in log scale, since its value may be negative.
  # The initial expression is from N(init.mean, init.sigma2)
  # alpha controls the gaussianity of the noise, the noise is sign(e)*(abs(e)^alpha) where e is N(0, err.sigma2)
  T <- delays;
  maxd <- max(delays);
  ng <- nrow(links);
  r <- matrix(data=rep(0,ng*(N+maxd)),nrow=(N+maxd),ncol=ng);
  # generate the initial expression
  for(i in 1:maxd) {
    r[i,] <- rnorm(ng, init.mean, init.sigma2)
  }
  #
  for(i in (maxd+1):(maxd+N)) {
    for(j in 1:ng) {
      x <- rnorm(1,0,err.sigma2);
      x <- sign(x)*(abs(x)^alpha);
      for(k in 1:ng) {
        if(T[k,j] != 0) {
          x <- x + r[i-T[k,j],k]*links[k,j];
        }
      }
      r[i,j] <- x;
    }
  }
  r[(maxd+1):(maxd+N),]
}

#
plot_exp <- function(r) {
  # r is a n by g matrix, where n is the number of time points,
  # g is the number of genes
  # The values are the expression of the genes at the different time points
  n <- nrow(r);
  g <- ncol(r);
  legend.r <- if(is.null(colnames(r))) rep("",g) else colnames(r);
  for(i in 1:g) {legend.r[i] <- paste("i=",i,legend.r[i]);}
  plot(x=c(1,n),y=range(r), type="n", main="Expressions",xlab="Time",ylab="Expression");
  for(i in 1:g) {
    lines(x=1:n, y=r[,i], type="b",col=i,pch=i);
  }
  legend(x="topright",legend=legend.r,pch=1:g,col=1:g);
}

##
gen_one_case <- function(n, n.small, max.parents, max.delay, n.points, err.s2, alpha) {
  # generate and simulate one random gene network the large network
  # consists of small networks, each with size n.small, and then they
  # are connected as if they are the nodes of a network.  returns a
  # list
  n.m <- ceiling(n/n.small);
  
  # generate the small networks
  rlinks <- matrix(data=0, nrow=n.m*n.small, ncol=n.m*n.small);
  rdelays <- matrix(data=0, nrow=n.m*n.small, ncol=n.m*n.small);
  for(i in 1:n.m) {
    # allow self loops, but have delays
    tmp.tn <- gen.stable.grn(n.small, FALSE, TRUE, max.parents, max.delay);
    idx <- (i-1)*n.small + 1:n.small;
    rlinks[idx,idx] <- tmp.tn$links;
    rdelays[idx,idx] <- tmp.tn$delays;
  }
  # connect the small networks
  overall <- gen.stable.grn(n.m, FALSE, TRUE, max.parents, max.delay);
  for(i in 1:n.m) {
    for(j in 1:n.m) {
      if(overall$delays[i,j] > 0) {
        f.idx <- sample(1:n.small, 2, replace=TRUE);
        rlinks[(i-1)*n.small + f.idx[1], (j-1)*n.small + f.idx[2]] <- overall$links[i,j];
        rdelays[(i-1)*n.small + f.idx[1], (j-1)*n.small + f.idx[2]] <- overall$delays[i,j];
      }
    }
  }
  #

  tn <- permute_grn(rlinks[1:n,1:n], rdelays[1:n,1:n]);
  sr <- sim_grn(tn$links, tn$delays, n.points, 1,1, err.s2, alpha);

  list(n=n, sn=n.small, mp=max.parents, md=max.delay, nps=n.points, err.sigma2=err.s2, alpha=alpha, true.net=tn, sim.exp=sr)
}

print.name <- function(d, r) {
  # only the base name, no nps
  paste("n",d$n,"mp",d$mp,"md",d$md,"e",10*d$err.sigma2, "a",10*d$alpha, "r",r, sep="");
}
##
##
# to generate an example GRN where there are communities, for preliminary work
# for reference only
eg.grn <- function() {
  max.parents <- 10;
  max.delay <- 4;
  sn <- 10;
  g0 <- gen.stable.grn(sn, FALSE,TRUE, max.parents, max.delay);
  g1 <- gen.stable.grn(sn, FALSE,TRUE, max.parents, max.delay);
  g2 <- gen.stable.grn(sn, FALSE,TRUE, max.parents, max.delay);
  g3 <- gen.stable.grn(sn, FALSE,TRUE, max.parents, max.delay);
  g4 <- gen.stable.grn(sn, FALSE,TRUE, max.parents, max.delay);
  # stitch them together
  rlinks <- matrix(data=0, nrow=sn*5, ncol=sn*5);
  rdelays <- matrix(data=0, nrow=sn*5, ncol=sn*5);

  rlinks[1:sn, 1:sn] <- g0$links;
  rdelays[1:sn, 1:sn] <- g0$delays;

  rlinks[sn+1:sn, sn+1:sn] <- g1$links;
  rdelays[sn+1:sn, sn+1:sn] <- g1$delays;

  rlinks[2*sn+1:sn, 2*sn+1:sn] <- g2$links;
  rdelays[2*sn+1:sn, 2*sn+1:sn] <- g2$delays;

  rlinks[3*sn+1:sn, 3*sn+1:sn] <- g3$links;
  rdelays[3*sn+1:sn, 3*sn+1:sn] <- g3$delays;

  rlinks[4*sn+1:sn, 4*sn+1:sn] <- g4$links;
  rdelays[4*sn+1:sn, 4*sn+1:sn] <- g4$delays;

  #
  # only a few links between communities
  # g0 --(0.7/1)--> g1
  idx <- sample(1:sn,2, replace=TRUE);
  rlinks[0*sn + idx[1], 1*sn + idx[2]] <- 0.7;
  rdelays[0*sn + idx[1], 1*sn + idx[2]] <- 1;
  
  # g0 --(0.31/4)--> g3
  idx <- sample(1:sn,2, replace=TRUE);
  rlinks[0*sn + idx[1], 3*sn + idx[2]] <- 0.31;
  rdelays[0*sn + idx[1], 3*sn + idx[2]] <- 4;
  
  # g0 --(0.8/2)--> g4
  idx <- sample(1:sn,2, replace=TRUE);
  rlinks[0*sn + idx[1], 4*sn + idx[2]] <- 0.8;
  rdelays[0*sn + idx[1], 4*sn + idx[2]] <- 2;
  
  # g1 --(-0.85/1)--> g2
  idx <- sample(1:sn,2, replace=TRUE);
  rlinks[1*sn + idx[1], 2*sn + idx[2]] <- 0.85;
  rdelays[1*sn + idx[1], 2*sn + idx[2]] <- 1;
  
  # g2 --(-0.55/2)--> g0
  idx <- sample(1:sn,2, replace=TRUE);
  rlinks[2*sn + idx[1], 0*sn + idx[2]] <- 0.55;
  rdelays[2*sn + idx[1], 0*sn + idx[2]] <- 2;
  
  # g2 --(0.64/1)--> g1
  idx <- sample(1:sn,2, replace=TRUE);
  rlinks[2*sn + idx[1], 1*sn + idx[2]] <- 0.64;
  rdelays[2*sn + idx[1], 1*sn + idx[2]] <- 1;
  
  # g3 --(-0.6/1)--> g4
  idx <- sample(1:sn,2, replace=TRUE);
  rlinks[3*sn + idx[1], 4*sn + idx[2]] <- 0.6;
  rdelays[3*sn + idx[1], 4*sn + idx[2]] <- 1;
  
  # g4 --(0.4/1)--> g2
  idx <- sample(1:sn,2, replace=TRUE);
  rlinks[4*sn + idx[1], 2*sn + idx[2]] <- 0.4;
  rdelays[4*sn + idx[1], 2*sn + idx[2]] <- 1;
  
  #
  list(links=rlinks, delays=rdelays)
}

#eg.grn1 <- eg.grn();

#eg.grn1.exp <- sim_grn(eg.grn1$links, eg.grn1$delays, 100, alpha=2);

print.exp <- function(d) {
  n <- nrow(d);
  for(i in 1:n) {
    cat(d[i,], sep=" ");
    cat("\n");
  }
}

#sink("test_data.txt");
#print.exp(eg.grn1.exp);
#sink();

print.true.grn <- function(links, delays) {
  # all indices 0-based
  n <- nrow(links);
  # links[i,j] is effect of i to j
  for(i in 1:n) {
    for(j in 1:n) {
      if(delays[i,j] > 0) {
        cat("To:",j-1, "From:",i-1, "Delay:",delays[i,j], "Effect:",links[i,j]);
        cat("\n");
      }
    }
  }
}

#sink("test_grn.txt");
#print.true.grn(eg.grn1$links, eg.grn1$delays);
#sink();

####
pairwise.abs.corr <- function(d, max.delay=4) {
  n <- ncol(d);
  m <- nrow(d);
  r <- matrix(0, nrow=n, ncol=n);
  for(i in 1:n) {
    for(j in 1:n) {
      for(k in 1:max.delay) {
        tmp <- abs(cor(x=d[1:(m-k),i], y=d[(1+k):m,j], method="pearson"));
        if(tmp > r[i,j]) {r[i,j] <- tmp;}
        if(tmp > r[j,i]) {r[j,i] <- tmp;}
      }
    }
  }
  #
  r
}

#eg.pair.cor <- pairwise.abs.corr(eg.grn1.exp, 4);

pairwise.cor.summary <- function(p, dn) {
  n <- nrow(p);
  k <- ceiling(n/dn);
  r <- matrix(0, nrow=k, ncol=k);
  for(i in 1:k) {
    ix <- (i-1)*dn;
    for(j in 1:k) {
      jx <- (j-1)*dn;
      r[i,j] <- mean(p[(1+ix):(dn+ix), (1+jx):(dn+jx)]);
    }
  }
  r
}

#eg.pair.cor.summary <- pairwise.cor.summary(eg.pair.cor, 10);

#eg.pair.cor.summary
###
source("infer.R")
source("pcor.R")

#eg.grn1.inferred <- infer_grn(eg.grn1.exp, 2,2,4, 1, "pcor", no.dup=TRUE);
####
grn.to.mat <- function(g) {
  n <- max(max(g$from), max(g$to));
  r <- matrix(0, nrow=n, ncol=n);
  for(i in 1:nrow(g)) {
    from.idx <- g$from[i];
    to.idx <- g$to[i];
    r[from.idx, to.idx] <- g$test.value[i];
  }
  r
}

#eg.grn1.inferred.mat <- grn.to.mat(eg.grn1.inferred);

#eg.grn1.inferred.mat

#pairwise.cor.summary(abs(eg.grn1.inferred.mat), 10);
####
print.inferred.links <- function(g) {
  for(i in 1:nrow(g)) {
    # 0-based indices
    # absolute test.value
    cat(g$from[i]-1, g$to[i]-1, abs(g$test.value[i]), sep="\t");
    cat("\n");
  }
}

#sink("test_g.txt");
#print.inferred.links(eg.grn1.inferred);
#sink();
####
gen.cases <- function() {
  out <- "test_cases/";
  mp <- 4;
  md <- 4;
  n.small <- 50;
  nps <- 200; # need only generate the large number that we will use, then take fragments if we want to test shorter series
  dir.create(out, mode="0755");
  for(n in c(500,1000)) {
    dir.create(paste(out,"n",n, sep=""), recursive=TRUE, mode="0755");
    for(err.s2 in c(0.5, 2, 8)) {
      for(alpha in c(0.5, 1, 2, 3)) {
        for(r in 1:20) { # replicates
          x <- gen_one_case(n, n.small, mp, md, nps, err.s2, alpha);
          name <- print.name(x,r);
          cat("Processing", name, "\n");
          # the true network
          sink(paste(out,"n",n,"/",name,"_grn.txt", sep=""));
          print.true.grn(x$true.net$links, x$true.net$delays);
          sink();
          # expression data, at different number of time points
          for(np in c(20,50,100,200)) {
            sink(paste(out,"n",n,"/",name,"_nps",np,".txt", sep=""));
            print.exp(x$sim.exp[1:np,]);
            sink();
          }
        }
      }
    }
  }
}

# gen.cases();
####

