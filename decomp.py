#!/usr/bin/python

# to read in expression data, and paritions, decompose the data
# according to the paritions, then call GRN inference tool to infer
# the subnetworks, finally combine the subnetworks.

import sys
import os

# In expression file, each line is data for one time point, and each
# line is separated by white space, each column is data for one gene.
# they are read as string, because we do not need the values.
def read_exp_data(file_name):
    Lns = []
    with open(file_name,'r') as f:
        for Ln in f:
            L = Ln.split()
            if len(L)>1:
                Lns.append(L)
    return Lns

#
def n_genes_from_exp_data(Lns):
    # Lns is as returned by read_exp_data
    return len(Lns[0])

#
def write_exp_data(file_name, data):
    with open(file_name,'w') as f:
        for x in data:
            f.write('\t'.join(x))
            f.write("\n")

# decompose, given a list of 0-based column indices
def exp_data_subset(exp_data, col_indices):
    return [[y[i] for i in col_indices] for y in exp_data]

## in the files containing partitions, we care about only lines
## starting with '= ', and has the relevant vertices information

# = parition 0, vertices: 4 12 19 21 22 23 24 25 26 27 28 29 43 49
# = parition 1, vertices: 13 14 17
# = parition 2, vertices: 11 13 16 18 19 38 39
# = parition 3, vertices: 10 17 19
# = parition 4, vertices: 10 11 12 16 17 18
# = parition 5, vertices: 3 8 30 31 32 33 34 37 38 39 47 49
# = parition 6, vertices: 32 36 37
# = parition 7, vertices: 0 8 16
# = parition 8, vertices: 10 12 13 15 19
# = parition 9, vertices: 1 2 5 7 9 26 27 28 29
# = parition 10, vertices: 3 4 6 7 11 15 24 30 35 39
# = parition 11, vertices: 1 20 23 36
# = parition 12, vertices: 5 21 36 38 40 41 43 44 45 47 48
# = parition 13, vertices: 19 29 34 35
# = parition 14, vertices: 7 9 43
# = parition 15, vertices: 40 42 46 47
# 

# return partition as list of list of 0-based indices
def read_vertex_partitions(file_name):
    Lns = []
    with open(file_name,'r') as f:
        for Ln in f:
            if "= " == Ln[0:2]:
                idx = Ln.rfind(":")
                L = ""
                if(idx != -1):
                    L = Ln[(idx+1):]
                else:
                    L = Ln
                idx1 = L.rfind("p")
                L1 = L
                L2 = ""
                if(idx1 != -1):
                    L1 = L[0:(idx1-1)]
                    L2 = L[(idx1+1):]
                Lns.append([[int(x) for x in L1.split()], [int(x) for x in L2.split()]])
    return Lns

#
# GRN in our MDL format
# Only care lines starting with "To:"

# To: 0 From: 2 Delay: 2 Coef: -0.562024
# To: 1 From: 2 Delay: 1 Coef: 0.72559
# To: 1 From: 2 Delay: 3 Coef: -0.309871
# To: 2 From: 1 Delay: 1 Coef: -0.935253
# To: 2 From: 2 Delay: 1 Coef: -0.143159
# To: 3 From: 4 Delay: 2 Coef: 0.229137
# To: 3 From: 3 Delay: 1 Coef: -0.0990973
# To: 4 From: 0 Delay: 2 Coef: 0.759474
# To: 4 From: 0 Delay: 1 Coef: 0.0332619
# 

# a grn is represented as a list of {'from':from, 'to':to, 'delay':delay, 'effect':effect} of the links
def get_link(L):
    # L is a list of strings, representing the to, from, delay and effect.
    to_id = 0;
    from_id = 0;
    delay = 0;
    effect = 0;

    n = len(L)
    i = 0
    while i < n:
        s = L[i]
        if s == 'To:':
            to_id = int(L[i+1])
        elif s == 'From:':
            from_id = int(L[i+1])
        elif s == 'Delay:':
            delay = int(L[i+1])
        elif s == 'Coef:' or s == 'Effect:':
            effect = float(L[i+1])
        i = i+2

    return {'from':from_id, 'to':to_id, 'delay':delay, 'effect':effect}

def read_grn(file_name):
    Lns = []
    with open(file_name,'r') as f:
        for Ln in f:
            L = Ln.split()
            if len(L)>1 and L[0]=='To:':
                Lns.append(get_link(L))
    return Lns

def write_grn(grn, file_name):
    with open(file_name,'w') as f:
        for x in grn:
            f.write("To: " + str(x['to']) + " From: " + str(x['from']))
            if 'effect' in x:
                f.write(" Effect: " + str(x['effect']))
            if 'delay' in x:
                f.write(" Delay: " + str(x['delay']))
            f.write('\n')

# grn is list of {'from':from, 'to':to, 'delay':delay,
# 'effect':effect} replace the 0-based from and to indices with the
# corresponding one in col_indices.
# But only keep those where the 'to' is < n_to_keep
def replace_idx(grn, col_indices, n_to_keep):
    Lns = []
    for x in grn:
        if(x['to'] < n_to_keep):
            L = x.copy()
            L['from'] = col_indices[L['from']]
            L['to'] = col_indices[L['to']]
            Lns.append(L)
    return Lns

#
# infer_grn(exp_filenames, out_grn_filename, out_prefix) should take
# care of inferring the subnetwork, by using the expression data in
# files exp_filenames, and output the GRN in MDL format in file
# out_grn_filename. out_prefix can be used as prefix for temporary
# files needed.
#
def decomp_grn(expfile_names, pafile_name,outfile_name, out_prefix, infer_grn, extra_args):
    partitions = read_vertex_partitions(pafile_name)
    exp_data = [read_exp_data(x) for x in expfile_names]
    # assume all replicates have the same number of genes
    n_genes = n_genes_from_exp_data(exp_data[0])
    out_grn = []
    i = 0
    for x in partitions:
        # if the size of partition is >= half of the number of genes, skip it
        # if the size is too small, skip it
        pa = x[0] + x[1]
        if len(pa) < (n_genes/2) and len(pa) > 5:
            data = [exp_data_subset(x, pa) for x in exp_data]
            sub_exp_filenames = [out_prefix + "_d" + str(i) + "_" + str(j) + ".txt" for j in range(len(expfile_names))]
            out_grn_filename = out_prefix + "_g" + str(i) + ".txt"
            for (n,d) in zip(sub_exp_filenames, data):
                write_exp_data(n, d)
            infer_grn(sub_exp_filenames, out_grn_filename, out_prefix, extra_args)
            sub_grn = read_grn(out_grn_filename)
            m_sub_grn = replace_idx(sub_grn, pa, len(x[0]))
            out_grn = out_grn + m_sub_grn
            # delete the temporary expression data
            for x in sub_exp_filenames:
                os.remove(x)
        i += 1
    write_grn(out_grn, outfile_name)

# The R version as written here does not support multiple replicates
# # For the R version
# # create tmp R file, and run it
# # the output GRN indices should be 0-based
# def R_clinde_infer_grn(exp_filename, out_grn_filename, out_prefix, extra_args):
#     tmpfile_name = out_prefix + "_tmp.R"
#     os.system("cat run_CLINDE.R > \"" + tmpfile_name + "\"")
#     with open(tmpfile_name,'a') as f:
#         f.write("\ntest.one.case(\"" + exp_filename + "\",\"" + out_grn_filename + "\");\n")
#     os.system("R -q --no-save < \"" + tmpfile_name + "\" >/dev/null")

# For the C version
# the output GRN indices should be 0-based
# the extra args should be list of strings, which allows overriding default parameters
def clinde_infer_grn(exp_filenames, out_grn_filename, out_prefix, extra_args):
    print "In clinde_infer_grn"
    tmpfile_name = out_prefix + "_tmp.txt"
    # roughly follow the parameter setting for large network
    # use only stage 2 results
    ns = ["\"" + x + "\"" for x in exp_filenames]
    os.system("./clinde -data " + (" ".join(ns)) + " -max.delay 4 -max.n 4 -method pcor -pruning all -no.dup " + (" ".join(extra_args)) + " > \"" + tmpfile_name + "\"")
    os.system("sed -n '/^==== Stage 2/,/^==== End Stage 2/p' \"" + tmpfile_name + "\" | grep \"^To:\" > \"" + out_grn_filename + "\"")


# For DD-Lasso
# create tmp R file, and run it
# the output GRN indices should be 0-based
def DD_lasso_infer_grn(exp_filenames, out_grn_filename, out_prefix, extra_args):
    print "In DD_lasso_infer_grn"
    ns = ["\"" + x + "\"" for x in exp_filenames]
    tmpfile_name = out_prefix + "_tmp.R"
    os.system("cat run_DD_lasso.R > \"" + tmpfile_name + "\"")
    with open(tmpfile_name,'a') as f:
        f.write("\ntest.one.case(c(" + (",".join(ns)) + "),\"" + out_grn_filename + "\");\n")
    os.system("R -q --no-save < \"" + tmpfile_name + "\" >/dev/null")

#
# usage: python decomp.py expfile_name pafile_name outfile_name outprefix infer_grn
# or
# usage: python decomp.py [ f1 f2 ... ] pafile_name outfile_name outprefix infer_grn
# note that to specify replicate files, add '[' and ']' between them, all separated by space
# currently only clinde_infer_grn for infer_grn
# extra arguments could be given to individual methods

def usage():
    print "Usage: decomp expfile_name pafile_name outfile_name outprefix infer_grn ..."
    print "or"
    print "Usage: decomp [ f1 f2 ... ] pafile_name outfile_name outprefix infer_grn ..."
    sys.exit(1)

def get_params(argv):
    if len(argv) < 6:
        usage()
    expfile_name = [argv[1]]
    res = argv[2:]
    if argv[1] == "[":
        ei = argv.index("]") # exception if not exists
        expfile_name = argv[2:(ei-1)]
        res = argv[(ei+1):]
    if len(res) < 4:
        usage()
    return {'names':expfile_name, 'pa':res[0], 'out':res[1], 'prefix':res[2], 'method':res[3].lower(), 'extra':res[4:]}

def main(argv):
    pars = get_params(argv)
    print pars
    # infer_grn can be one of clinde, DD_lasso (case insensitive)
    infer_grn = clinde_infer_grn
    s_infer_grn = pars['method']
    if s_infer_grn == "dd_lasso":
        print "Use DD_Lasso"
        infer_grn = DD_lasso_infer_grn
    decomp_grn(pars['names'], pars['pa'], pars['out'],pars['prefix'], infer_grn, pars['extra'])

#

main(sys.argv)
#
