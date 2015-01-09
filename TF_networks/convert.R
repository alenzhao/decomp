
# convert the format of GDS4238 data

#
interpolate.d <- function(d, xp, tp) {
  # d is the data frame, xp is x of the input d, tp is the target x
  # returns a list of two matrix where the rows are the interpolated versions in d
  # one is interpolated raw value, the other is interpolated log10 values
  nr <- nrow(d);
  nc <- length(tp);
  r1 <- matrix(0, nrow=nr, ncol=nc);
  r2 <- matrix(0, nrow=nr, ncol=nc);
  for(i in 1:nr) {
    tmp <- spline(x=xp, y=d[i,], xout=tp);
    r1[i,] <- tmp$y;
    tmp2 <- spline(x=xp, y=log10(d[i,]), xout=tp);
    r2[i,] <- tmp2$y;
  }
  list(raw=r1, log10=r2)
}

output.d <- function(d, out, name) {
  # d is a list of raw and log10
  nc <- ncol(d$raw);
  tc <- 32; # truncate up to 8h, start from 0.25, increment by 0.25
  # raw, all
  sink(paste(out,"_",name,"_raw_all.txt", sep=""));
  for(i in 1:nc) {cat(d$raw[,i], sep="\t"); cat("\n");}
  sink();
  # raw, truncated
  sink(paste(out,"_",name,"_raw_cut.txt", sep=""));
  for(i in 1:tc) {cat(d$raw[,i], sep="\t"); cat("\n");}
  sink();
  # log10, all
  sink(paste(out,"_",name,"_log10_all.txt", sep=""));
  for(i in 1:nc) {cat(d$log10[,i], sep="\t"); cat("\n");}
  sink();
  # log10, truncated
  sink(paste(out,"_",name,"_log10_cut.txt", sep=""));
  for(i in 1:tc) {cat(d$log10[,i], sep="\t"); cat("\n");}
  sink();
}

convert.one.data <- function(d, out) {
  # d is the data, out is the output prefix
# interpolation
# truncate or not truncate
# log2 or raw
# save

# time points
tp1 <- c(0.25, 0.5, 1, 1.5, 2, 4, 6, 8, 12, 18);
tp2 <- c(0.25, 0.5, 1, 1.5, 4, 6, 8, 12, 18); ## for LTX, with no data for 2h
tp <- seq(from=0.25, to=18, by=0.25);

# get the columns for each type and rep
# IFNb
# 	rep1
# 		0.25h: GSM528681
# 		0.5h : GSM528683
# 		1h   : GSM528687
# 		1.5h : GSM528685
# 		2h   : GSM528693
# 		4h   : GSM528695
# 		6h   : GSM528697
# 		8h   : GSM528700
# 		12h  : GSM528689
# 		18h  : GSM528691
IFNb_rep1 <- interpolate.d(d[,c("GSM528681","GSM528683","GSM528687","GSM528685","GSM528693","GSM528695","GSM528697","GSM528700","GSM528689","GSM528691")], tp1,tp);
output.d(IFNb_rep1, out, "IFNb_rep1");
# 	rep2
# 		0.25h: GSM528682
# 		0.5h : GSM528684
# 		1h   : GSM528688
# 		1.5h : GSM528686
# 		2h   : GSM528694
# 		4h   : GSM528696
# 		6h   : GSM528698
# 		8h   : GSM528699
# 		12h  : GSM528690
# 		18h  : GSM528692
IFNb_rep2 <- interpolate.d(d[,c("GSM528682","GSM528684","GSM528688","GSM528686","GSM528694","GSM528696","GSM528698","GSM528699","GSM528690","GSM528692")], tp1,tp);
output.d(IFNb_rep2, out, "IFNb_rep2");
# PR8pre
# 	rep1
# 		0.25h: GSM528779
# 		0.5h : GSM528780
# 		1h   : GSM528782
# 		1.5h : GSM528781
# 		2h   : GSM528785
# 		4h   : GSM528786
# 		6h   : GSM528787
# 		8h   : GSM528788
# 		12h  : GSM528783
# 		18h  : GSM528784
PR8pre <- interpolate.d(d[,c("GSM528779","GSM528780","GSM528782","GSM528781","GSM528785","GSM528786","GSM528787","GSM528788","GSM528783","GSM528784")], tp1,tp);
output.d(PR8pre, out, "PR8pre");
# PR8post
# 	rep1
# 		0.25h: GSM528759
# 		0.5h : GSM528761
# 		1h   : GSM528766
# 		1.5h : GSM528763
# 		2h   : GSM528771
# 		4h   : GSM528773
# 		6h   : GSM528775
# 		8h   : GSM528777
# 		12h  : GSM528767
# 		18h  : GSM528769
PR8post_rep1 <- interpolate.d(d[,c("GSM528759","GSM528761","GSM528766","GSM528763","GSM528771","GSM528773","GSM528775","GSM528777","GSM528767","GSM528769")], tp1,tp);
output.d(PR8post_rep1, out, "PR8post_rep1");
# 	rep2
# 		0.25h: GSM528760
# 		0.5h : GSM528762
# 		1h   : GSM528765
# 		1.5h : GSM528764
# 		2h   : GSM528772
# 		4h   : GSM528774
# 		6h   : GSM528776
# 		8h   : GSM528778
# 		12h  : GSM528768
# 		18h  : GSM528770
PR8post_rep2 <- interpolate.d(d[,c("GSM528760","GSM528762","GSM528765","GSM528764","GSM528772","GSM528774","GSM528776","GSM528778","GSM528768","GSM528770")], tp1,tp);
output.d(PR8post_rep2, out, "PR8post_rep2");
# delNS1pre
# 	rep1
# 		0.25h: GSM528671
# 		0.5h : GSM528672
# 		1h   : GSM528674
# 		1.5h : GSM528673
# 		2h   : GSM528677
# 		4h   : GSM528678
# 		6h   : GSM528679
# 		8h   : GSM528680
# 		12h  : GSM528675
# 		18h  : GSM528676
delNS1pre <- interpolate.d(d[,c("GSM528671","GSM528672","GSM528674","GSM528673","GSM528677","GSM528678","GSM528679","GSM528680","GSM528675","GSM528676")], tp1,tp);
output.d(delNS1pre, out, "delNS1pre");
# delNS1post
# 	rep1
# 		0.25h: GSM528651
# 		0.5h : GSM528653
# 		1h   : GSM528657
# 		1.5h : GSM528655
# 		2h   : GSM528663
# 		4h   : GSM528665
# 		6h   : GSM528667
# 		8h   : GSM528669
# 		12h  : GSM528659
# 		18h  : GSM528661
delNS1post_rep1 <- interpolate.d(d[,c("GSM528651","GSM528653","GSM528657","GSM528655","GSM528663","GSM528665","GSM528667","GSM528669","GSM528659","GSM528661")], tp1,tp);
output.d(delNS1post_rep1, out, "delNS1post_rep1");
# 	rep2
# 		0.25h: GSM528652
# 		0.5h : GSM528654
# 		1h   : GSM528658
# 		1.5h : GSM528656
# 		2h   : GSM528664
# 		4h   : GSM528666
# 		6h   : GSM528668
# 		8h   : GSM528670
# 		12h  : GSM528660
# 		18h  : GSM528662
delNS1post_rep2 <- interpolate.d(d[,c("GSM528652","GSM528654","GSM528658","GSM528656","GSM528664","GSM528666","GSM528668","GSM528670","GSM528660","GSM528662")], tp1,tp);
output.d(delNS1post_rep2, out, "delNS1post_rep2");
# vRNA LTX+RNA
# 	rep1
# 		0.25h: GSM528701
# 		0.5h : GSM528703
# 		1h   : GSM528707
# 		1.5h : GSM528705
# 		2h   : GSM528713
# 		4h   : GSM528715
# 		6h   : GSM528717
# 		8h   : GSM528720
# 		12h  : GSM528709
# 		18h  : GSM528711
vRNA_rep1 <- interpolate.d(d[,c("GSM528701","GSM528703","GSM528707","GSM528705","GSM528713","GSM528715","GSM528717","GSM528720","GSM528709","GSM528711")], tp1,tp);
output.d(vRNA_rep1, out, "vRNA_rep1");
# 	rep2
# 		0.25h: GSM528702
# 		0.5h : GSM528704
# 		1h   : GSM528708
# 		1.5h : GSM528706
# 		2h   : GSM528714
# 		4h   : GSM528716
# 		6h   : GSM528718
# 		8h   : GSM528719
# 		12h  : GSM528710
# 		18h  : GSM528712
vRNA_rep2 <- interpolate.d(d[,c("GSM528702","GSM528704","GSM528708","GSM528706","GSM528714","GSM528716","GSM528718","GSM528719","GSM528710","GSM528712")], tp1,tp);
output.d(vRNA_rep2, out, "vRNA_rep2");
# LTX # need special handling
# 	rep1
# 		0.25h: GSM528721
# 		0.5h : GSM528723
# 		1h   : GSM528727
# 		1.5h : GSM528725
# 		2h   :
# 		4h   : GSM528733
# 		6h   : GSM528735
# 		8h   : GSM528737
# 		12h  : GSM528729
# 		18h  : GSM528731
LTX_rep1 <- interpolate.d(d[,c("GSM528721","GSM528723","GSM528727","GSM528725", "GSM528733","GSM528735","GSM528737","GSM528729","GSM528731")], tp2,tp);
output.d(LTX_rep1, out, "LTX_rep1");
# 	rep2
# 		0.25h: GSM528722
# 		0.5h : GSM528724
# 		1h   : GSM528728
# 		1.5h : GSM528726
# 		2h   :
# 		4h   : GSM528734
# 		6h   : GSM528736
# 		8h   : GSM528738
# 		12h  : GSM528730
# 		18h  : GSM528732
LTX_rep2 <- interpolate.d(d[,c("GSM528722","GSM528724","GSM528728","GSM528726", "GSM528734","GSM528736","GSM528738","GSM528730","GSM528732")], tp2,tp);
output.d(LTX_rep2, out, "LTX_rep2");
# MockIFN
# 	rep1
# 		0.25h: GSM528739
# 		0.5h : GSM528741
# 		1h   : GSM528745
# 		1.5h : GSM528743
# 		2h   : GSM528751
# 		4h   : GSM528753
# 		6h   : GSM528756
# 		8h   : GSM528757
# 		12h  : GSM528747
# 		18h  : GSM528749
MockIFN_rep1 <- interpolate.d(d[,c("GSM528739","GSM528741","GSM528745","GSM528743","GSM528751","GSM528753","GSM528756","GSM528757","GSM528747","GSM528749")], tp1,tp);
output.d(MockIFN_rep1, out, "MockIFN_rep1");
# 	rep2
# 		0.25h: GSM528740
# 		0.5h : GSM528742
# 		1h   : GSM528746
# 		1.5h : GSM528744
# 		2h   : GSM528752
# 		4h   : GSM528754
# 		6h   : GSM528755
# 		8h   : GSM528758
# 		12h  : GSM528748
# 		18h  : GSM528750
MockIFN_rep2 <- interpolate.d(d[,c("GSM528740","GSM528742","GSM528746","GSM528744","GSM528752","GSM528754","GSM528755","GSM528758","GSM528748","GSM528750")], tp1,tp);
output.d(MockIFN_rep2, out, "MockIFN_rep2");
#
}


convert.all <- function() {
  for(x in c("AG10803", "AoAF", "CD20+", "CD34+_Mobilized", "fBrain", "fHeart", "fLung", "GM06990", "GM12865", "H7-hESC", "HAEpiC", "HA-h", "HCF", "HCM", "HCPEpiC", "HEEpiC", "HepG2", "HFF", "HIPEpiC", "HMF", "HMVEC-dBl-Ad", "HMVEC-dBl-Neo", "HMVEC-dLy-Neo", "HMVEC_LLy", "HPAF", "HPdLF", "HPF", "HRCEpiC", "HSMM", "HVMF", "IMR90", "K562", "NB4", "NH-A", "NHDF-Ad", "NHDF-neo", "NHLF", "SAEC", "SKMC", "SK-N-SH_RA", "Th1")) {
    # read in a file
    cat(x,"\n");
    tmp.d <- read.delim(paste(x,"/GDS4238_subset2.txt", sep=""));
    convert.one.data(tmp.d, paste(x,"/GDS4238", sep=""));
  }
}
