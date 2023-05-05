# third party modules
import numpy as np

def balancedClassIndices(y, classes, max_percentage=1.0, shuffle=True, return_mask=True):
    """ create a balanced index set from a vector of labels """
    idxs = []
    for c in classes:
        idxs.append(y == c)
    
    # find maximum number of class labels to keep
    nb = int(max_percentage*min([idx.sum() for idx in idxs]))
    if nb == 0:
        # use all vertices
        if return_mask:
            return (y >= 0)
        else:
            return y
    
    for i in range(len(idxs)):
        # exclude any indices above nb
        idx = np.argwhere(idxs[i]).flatten()
        if shuffle :
            np.random.shuffle(idx)
        idxs[i][idx[nb:]] = False
    
    idxb = np.array(idxs, dtype=bool).sum(axis=0, dtype=bool) # a balanced index set
    if return_mask:
        return idxb
    else:
        return np.argwhere(idxb).flatten()
