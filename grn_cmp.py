#!/usr/bin/python

# to read in the true grn and the predicted grn in Banjo's report format
# and calculate the link and delay Recall, Precision and F-measure.

# tested in Python 2.7.3, seems to work

import sys

# in our MDL format
# sample

# Pseudorandom seed: 1394442261
# Reading expression data from test_case_nps20/mdl_test_case_nps20e10a5r1.txt ... Done
# Read 20 time points, 5 genes.
# Best found for gene 0
# Score: 332
# Target: 0 Offset: 0.94286 Links: 1
# To: 0 From: 2 Delay: 2 Coef: -0.562024
# Best found for gene 1
# Score: 348
# Target: 1 Offset: 0.029339 Links: 2
# To: 1 From: 2 Delay: 1 Coef: 0.72559
# To: 1 From: 2 Delay: 3 Coef: -0.309871
# Best found for gene 2
# Score: 331
# Target: 2 Offset: -0.000291244 Links: 2
# To: 2 From: 1 Delay: 1 Coef: -0.935253
# To: 2 From: 2 Delay: 1 Coef: -0.143159
# Best found for gene 3
# Score: 330
# Target: 3 Offset: 0.036674 Links: 2
# To: 3 From: 4 Delay: 2 Coef: 0.229137
# To: 3 From: 3 Delay: 1 Coef: -0.0990973
# Best found for gene 4
# Score: 345
# Target: 4 Offset: -0.00017839 Links: 2
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

def sign(a):
    if a < 0:
        return -1
    elif a > 0:
        return 1
    else:
        return 0

def link_correct(L1,L2):
    return L1['to']==L2['to'] and L1['from']==L2['from']

def delay_correct(L1,L2):
    return link_correct(L1,L2) and L1['delay']==L2['delay']

def effect_correct(L1,L2):
    return link_correct(L1,L2) and sign(L1['effect'])==sign(L2['effect'])

def correct_in(L, Lns, correct):
    for x in Lns:
        if correct(L,x):
            return True
    return False

def prop_correct(pred,lns, correct):
    # proportion of links in pred that is same as one in lns, as judged by correct
    # use naive method, since the lists would be short in our case
    i = 0
    for L in pred:
        if correct_in(L, lns, correct):
            i = i+1
    if len(pred)>0:
        return float(i)/len(pred)
    return 0

def cal_fmeasure(R,P):
    if R < 1.0e-8 or P < 1.0e-8:
        return 0
    return 2.0*R*P/(R + P)

def cmp_grn(true_grn, pred_grn):
    # prints the links.recall, links.precision, links.fmeasure, delay.recall, delay.precision, delay.fmeasure
    links_recall = prop_correct(true_grn, pred_grn, link_correct)
    links_precision = prop_correct(pred_grn, true_grn, link_correct)
    links_fmeasure = cal_fmeasure(links_recall, links_precision)
    delay_recall = prop_correct(true_grn, pred_grn, delay_correct)
    delay_precision = prop_correct(pred_grn, true_grn, delay_correct)
    delay_fmeasure = cal_fmeasure(delay_recall, delay_precision)
    effect_recall = prop_correct(true_grn, pred_grn, effect_correct)
    effect_precision = prop_correct(pred_grn, true_grn, effect_correct)
    effect_fmeasure = cal_fmeasure(effect_recall, effect_precision)
    # print a one-line summary
    print str(links_recall) + "," + str(links_precision) + "," + str(links_fmeasure) + "," + str(delay_recall) + "," + str(delay_precision) + "," + str(delay_fmeasure) + "," + str(effect_recall) + "," + str(effect_precision) + "," + str(effect_fmeasure) + "\n"

def main(argv):
    # usage: python grn_cmp.py predicted_grn_file true_grn_file
    # predicted file are true grn are in MDL format
    pred_f = argv[1]
    true_f = argv[2]
    grn_pred = read_grn(pred_f)
    grn_true = read_grn(true_f)
    cmp_grn(grn_true, grn_pred)

main(sys.argv)
