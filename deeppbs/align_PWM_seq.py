import numpy as np
from deeppbs.nn.metrics import IC_weighted_PCC

def ungappedAlign(ml,  ms, gt): ## ungapped alignment maximising dot product, needs Length X 4 arrays    
        
    max_score = -9999
    opt_i = 0
    opt_j = 0
    opt_k = 0
    l = ml.shape[0]
    s = ms.shape[0]
    for i in range(0,s):
        for k in range(1, s-i+1):   ### k overlap length
            
            for j in range(0,l - k+1):
                #score = np.sum(ms[i:i+k,:]* ml[j:j+k,:]) # dot product scoring
                
                score = 0

                if np.array_equal(ml, gt): #IC_weighted_PCC scoring
                    idx = j
                else:
                    idx = i
                for col in range(k):
                    col_score, _ = IC_weighted_PCC(ms[i:i+k,:][col,:], ml[j:j+k,:][col,:], gt=gt[idx:idx+k,:][col,:])
                    score += col_score

                
                
                if(score > max_score):
                    max_score = score
                    opt_i = i
                    opt_j = j
                    opt_k = k
    return opt_i, opt_j, opt_k, max_score

def alignPWMSeq(pwm, seq):
    if(pwm.shape[0] > seq.shape[0]):
        seq_start, pwm_start, k, max_score = ungappedAlign(pwm, seq, pwm)
    else:
        pwm_start, seq_start, k, max_score = ungappedAlign(seq, pwm, pwm)
        
    return pwm_start, seq_start, k, max_score

