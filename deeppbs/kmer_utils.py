
from itertools import product
import re

def countKmers(seq, k=5):
    seq = re.sub(r"N*","",seq)
    counts = {}
    num_kmers = len(seq) - k + 1
    for i in range(num_kmers):
        kmer = seq[i:i+k]
        if kmer not in counts:
            counts[kmer] = 0
        counts[kmer] += 1
    return counts



def rc(kmer):
    rc_dict = {
            'A': 'T',
            'C': 'G',
            'G': 'C',
            'T': 'A'
            }
    rc_kmer = ''
    for base in kmer[::-1]:
        rc_kmer += rc_dict[base]

    return rc_kmer


def getLabelDict(k=5):
    all_kmers = [''.join(c) for c in list(product('ACGT', repeat=k))]
    label_dict = dict()
    count = 0
    for kmer in all_kmers:
        if(kmer not in label_dict.keys()):
            rc_kmer = rc(kmer)
            if(rc_kmer not in label_dict.keys()):
                label_dict[kmer] = count
                count +=1
    return label_dict

def getRcSeparatedLabelDict(k=5):
    all_kmers = [''.join(c) for c in list(product('ACGT', repeat=k))]
    label_dict = dict()
    count = 0
    for kmer in all_kmers:
        if(kmer not in label_dict.keys()):
            label_dict[kmer] = count
            count +=1
    return label_dict
