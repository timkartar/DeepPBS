import numpy as np

base_dict = {0: 'A',
        1:'C',
        2:'G',
        3:'T'
        }

def sampleFromPwm(pwm):
    likelihood = 1
    seq = ""
    for base_probs in pwm:
        idx = np.random.choice([0,1,2,3],p = base_probs)
        likelihood = likelihood * base_probs[idx]
        seq += base_dict[idx]
    
    return seq, likelihood
        
