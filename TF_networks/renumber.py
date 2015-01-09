#!/usr/bin/python

# to renumber the TF networks from gene names to number, with respect
# to a name list

import sys

def read_names(file_name):
    ns = {}
    i = 0
    with open(file_name, 'r') as f:
        for Ln in f:
            ns[Ln.strip()] = str(i)
            i += 1
    return ns

def convert_links(file_name, name_ids):
    with open(file_name, 'r') as f:
        for Ln in f:
            L = Ln.split()
            if(len(L))==2:
                # the delay and effect are not in fact known,
                # so just use dummy value of 1.
                # the link is from column 2 pointing to column 1
                n_to = L[0]
                n_fr = L[1]
                if n_to not in name_ids:
                    print "Unknown name: " + n_to
                    sys.exit(1)
                if n_fr not in name_ids:
                    print "Unknown name: " + n_fr
                    sys.exit(1)
                print "To: " + name_ids[n_to] + " From: " + name_ids[n_fr] + " Delay: 1 Effect: 1"

def main(argv):
    # usage: python renumber.py infile names
    # infile has two columns,
    # names has one column of names.
    if len(argv) != 3:
        print "Usage: renumber infile names"
        sys.exit(1)
    convert_links(argv[1], read_names(argv[2]))

if __name__ == '__main__':
    main(sys.argv)

