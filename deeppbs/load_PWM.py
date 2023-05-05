import pickle
import numpy as np
import os
from scipy.stats import entropy

DATA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_data/pwms.pickle")


def loadPWM(key):
    pwm_dict = pickle.load(open(DATA_PATH,"rb"))
    untrimmed = np.array(list(pwm_dict[key].pwm.values())).T 
    start = 0
    for i in range(untrimmed.shape[0]):
        ent = entropy(untrimmed[i,:],[0.25,0.25,0.25,0.25], base=2)
        if ent > 0.5:
            break
        else:
            start = i
    end = untrimmed.shape[0]
    for i in range(untrimmed.shape[0]): 
        ent = entropy(untrimmed[untrimmed.shape[0] - i - 1,:],[0.25,0.25,0.25,0.25], base=2)
        if ent > 0.5:
            break
        else:
            end = untrimmed.shape[0] - i
    return  untrimmed[start:end,:]
    
