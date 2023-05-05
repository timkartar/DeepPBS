import numpy as np

def seqToOneHot(seq):
    mapping = dict(zip("ACGT", range(4)))
    seq2 = [mapping[i] for i in seq]
    return np.eye(4)[seq2]

def oneHotToSeq(one_hot):
    alphabet = { 0 : "A",
            1 : "C",
            2 : "G",
            3 : "T" }
    arg = np.argmax(one_hot, axis = 1).tolist()
    return "".join([alphabet[i] for i in arg])

def rcSeq(seq):
    mapping = { "A" : "T",
            "T": "A",
            "G": "C",
            "C": "G"}
    rc = ""
    for c in reversed(seq):
        rc = rc + mapping[c]

    return rc

