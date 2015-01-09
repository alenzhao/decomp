# decomp

Infer large causal gene regulatory network from time series data by
decomposing into smaller subnetworks.

====================================================================
Introduction
====================================================================

This is a set of programs for inferring causal large gene regulatory
network (GRN) from time series expression data. To alleviate the
problem of insufficient number of time points relative to the number
of genes, we propose to decompose the large GRN into smaller
overlapping subnetworks without using prior information, and then
infer each subnetwork, and finally merge the subnetworks to give the
final GRN. These programs are written mainly for the purpose of
evaluating the effectiveness of this strategy, so may not be as
user-friendly as desired. The programs have been tested in Debian and
Ubuntu Linux distributions, and should build and work with minimal
tweaking in other Linux distributions. In the following, the main
steps are described briefly.

====================================================================
Step 1. Initial GRN using CLINDE
====================================================================

The initial GRN is inferred using CLINDE/clinde from time series
expression data. 

==== Building

Before buliding CLINDE/clinde, you need to build gsl, which is the GNU
Scientific Library, which contains routines used by CLINDE.

gsl-1.16.tar.gz is included for convenience, and we have tested using
gsl-1.16. First extract it by:

    tar xzvf gsl-1.16.tar.gz

Then to build it, typically you can use:

    cd gsl-1.16/
    ./configure
    make
    make install

But refer to gsl-1.16/README for trouble shooting.

After that, you could build clinde by using the provided makefile in CLINDE:

    cd CLINDE/
    make

You may also build it directly using:

    cd CLINDE/
    gcc -Wall -O2 -static clinde.c parse_option.c tsv.c -o clinde -lgsl -lgslcblas -lm

==== Usage:

The usage of CLINDE/clinde is:

    Usage: ./clinde [-?] -data item1 [item2 ...] [-st1 real] [-st2 real]
    [-max.delay int] [-max.n int] [-method str] [-pruning str] [-one.delay]
    [-no.dup]
    
    Description of the options:
      -?:  Showing the usage.
    
      -data [REQUIRED]:  File name(s) of the input expression data (tab/space
        separated), each row is a time point, each column is a gene. Each file
        should have the same number of columns, but may have different number of
        rows.
    
      -st1:  Score threshold for stage 1. For method=pcor, score is -log10(p),
        where p is the intended p-value. For method=mi, score is the mutual
        information. Default 2.
    
      -st2:  Score threshold for stage 2. Similar to st1. Default 2.
    
      -max.delay:  Maximum delay (lags) to use in the inference, default 4.
    
      -max.n:  Maximum number of parents to condition on in stage 2, default 4.
    
      -method:  Method of testing links in stage 1 and 2, can be pcor (partial
        correlation) or mi (mutual information). Default "pcor"
    
      -pruning:  Pruning strategy in stage 2, can be all (consider all neighbors
        of the two vertices when doing conditional test of the link between the
        two) or common (consider only common neighbors of the two vertices).
        Default "all"
    
      -one.delay:  To keep only one delay (the one with best score, smallest
        delay) for each link after stage 1. Default false.
    
      -no.dup:  Remove duplicate links. To keep only one delay (the one with best
        score, smallest delay) for each link after stage 2. Default false.
    

==== Example usage:

An example use for inferring initial GRN is:

    CLINDE/clinde -data CLINDE/n500_case_study_nps50a.txt -no.dup > output.txt

Now output.txt contains some messages and the GRN after stage 1 and
stage 2 of CLINDE. To extract the GRN after stage 2 for use in the
next step, you may use sed and grep, which are available in Linux:

    sed -n '/^==== Stage 2/,/^==== End Stage 2/p' output.txt | grep "^To:" > outgrn.txt

Now outgrn.txt contains the GRN in a format acceptable for the
following step.

====================================================================
Step 2. Decomposition
====================================================================

The next step is to decompose the intial GRN into overlapping
subnetworks, the program is ovc, which uses "edge betweenness" to
decompose a large network into smaller ones by successively deleting
the edge with the highest "edge betweenness" in a component, until all
the components are small enough.

==== Building

You may build ovc by using the provided makefile by:

    make ovc

Or you may directly build it by:

    gcc -Wall -O3 ovc.c parse_option.c -o ovc -lm

==== Usage:

    Usage: ./ovc [-?] -f str [-s int]
    
    Description of the options:
      -?:  Showing the usage.
    
      -f [REQUIRED]:  File name of the input edges. Each line in the file
        represents an edge, which consists of a 0-based 'from' index, a 0-based
        'to' index, and a weight, separated by space..
    
      -s:  Threshold of component size. If a component is no larger than this, it
        is not further divided. Default 60.
    
==== Example usage:

An example use for decomposition is (continuing the step 1, and using
the default parameter):

    ./ovc -f outgrn.txt > outp.txt

Now outp.txt contains some messages and the subnetworks. There are 3
types of subnetworks: "component" which are disjoint, "partition"
which is component with its parents, and "xpartition" which is same as
"partition" except the parents are separated after a letter 'p' in the
output. The three types of subnetworks are given in lines starting
with "= component", "= partition" and "= xpartition" respectively.

In order to extract them, you may use grep again:

    grep "^= component" outp.txt > outcomponents.txt
    grep "^= partition" outp.txt > outpartition.txt
    grep "^= xpartition" outp.txt > outxpartition.txt

====================================================================
Step 3. Infer Subnetworks using either CLINDE or DD-lasso and merge
into final GRN
====================================================================

The program to infer each subnetwork and merge them into the final GRN
is decomp.py, which is a Python program tested under Python 2.7.

==== Building:

There is no need to build the program, but you need to install python
if it is not already installed. Please refer to your Linux's
instructions on how to install Python 2.7.

==== Usage:

    Usage: python decomp.py expfile_name pafile_name outfile_name outprefix infer_grn ...

or

    Usage: python decomp.py [ f1 f2 ... ] pafile_name outfile_name outprefix infer_grn ...

where

    expfile_name is the time series file if you are using only one.

    f1, f2 ... are the time series files if you are using multiple time series.

    pafile_name is the file of the subnetworks, i.e. that extracted in step 2.

    outfile_name is the desired filename of the output final GRN.

    outprefix is the prefix of temporary files in the process of inferring the subnetworks.

    infer_grn is the method for inferring each subnetwork, which can be "dd_lasso" for DD-lasso in DD_lasso/, or "clinde" for CLINDE/clinde.

    ... means passing any extract parameters to DD-lasso or CLINDE, which is mainly used for CLINDE, e.g. to change the score threshold in inferring the subnetworks.

Note that decomp.py assumes the presence of CLINDE/clinde if clinde is
used, and the following R files in DD_lasso if dd_lasso is used:

      DD_lasso/dd_lasso.R
      DD_lasso/run_DD_lasso.R

Also, R has to be installed, and the R package "lars" has to be
installed if DD-lasso is to be used. Please refer to your Linux's
instructions on how to install R, and instructions of R on how to
install the lars package. We have tested on R 2.15.1 and lars 1.2.

==== Example usage:

An example use is (continuing step 2 above):

    python decomp.py CLINDE/n500_case_study_nps50a.txt outpartition.txt outmerged.txt tmp clinde

====================================================================
Step 4. Comparison of two GRNs (if applicable)
====================================================================

If the true GRN is known, you may assess how close the predicted GRN
is to the true GRN by using either grn_cmp.py or grn_cmp_hcc.

grn_cmp.py is again a Python program, and is slower than grn_cmp_hcc,
but functionally the same.

==== Usage of grn_cmp.py:

    usage: python grn_cmp.py predicted_grn_file true_grn_file

where both predicted_grn_file and true_grn_file are in the same format
as outgrn.txt above.

This program outputs a line of 9 numbers, which are:

    links_recall,links_precision,links_f_measure,
    delay_recall,delay_precision,delay_f_measure,
    effect_recall,effect_precision,effect_f_measure

where a predicted link is correct iff both the end points and
     direction is correct; delays is correct iff both link is correct
     the predicted delay is correct; effect is correct iff both link
     is correct and the sign of the coefficient is correct.

==== Building grn_cmp_hcc:

You may build it by using the provided makefile:

    make grn_cmp_hcc

Or you may build it directly by:

    gcc -Wall -O3 grn_cmp_hcc.c parse_option.c -o grn_cmp_hcc -lm

==== Usage of grn_cmp_hcc:

    Usage: ./grn_cmp_hcc [-?] -p str -t str [-n int] [-v]
    
    Description of the options:
      -?:  Showing the usage.
    
      -p [REQUIRED]:  File name of the predicted GRN. Each line in the file
        represents an edge, which consists of a 0-based 'from' index, a 0-based
        'to' index, a delay, and the effect, separated by space..
    
      -t [REQUIRED]:  File name of the true GRN. Each line in the file represents
        an edge, which consists of a 0-based 'from' index, a 0-based 'to' index, a
        delay, and the effect, separated by space..
    
      -n:  Number of non-hidden genes, unspecified if <= 0. If an index is >=
        this, it is regarded as a hidden node. Default 0.
    
      -v:  Verbose mode.


The parameter -n can be ignored because there are no hidden genes
handled in these programs. The output of grn_cmp_hcc has the same
format as that of grn_cmp.py (if -v is NOT used), but grn_cmp_hcc is
in general much faster.

====================================================================
R Scripts for Generation of Synthetic Data (Optional)
====================================================================

The function "gen.cases()" in synthetic.R can help generate a lot of
synthetic data for testing purposes. Of course R need to be installed
first.

Be warned that a LOT of synthetic data will be generated, and it may
take some time (up to a few hours on a slow computer) because R is not
very fast.

You may tweak "gen.cases()" to generate cases with other parameters.

====================================================================
====================================================================
====================================================================

