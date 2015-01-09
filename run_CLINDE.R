###
## To test our algorithm on synthetic data previously generated.
##

# assume test.n, test.p, test.m, test.max.n are defined

source("sim.R");
source("infer.R");

#
print.grn <- function(g) {
  # all use 0-based indices
  n <- nrow(g);
  for(i in 1:n) {
    cat("To:",g$to[i]-1, "From:",g$from[i]-1, "Delay:",g$delay[i], "Effect:",g$test.value[i], sep="\t");
    cat("\n");
  }
}

# default parameters
test.max.n <- 4;
test.md <- 4;
test.st1 <- 2;
test.st2 <- 2;
test.method <- "pcor";
test.pruning <- "all";

test.one.case <- function(exp.filename, out.grn.filename) {
  sr.t <- read.table(file=exp.filename, header=FALSE);
  sr <- data.matrix(sr.t);
  grn <- infer_grn(sr, test.st1,test.st2,test.md,test.max.n, method=test.method,pruning=test.pruning, no.dup=TRUE);
  sink(out.grn.filename);
  print.grn(grn);
  sink();
}

# The calling is left to external script
####


