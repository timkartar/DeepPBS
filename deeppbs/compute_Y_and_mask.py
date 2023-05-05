import numpy as np
from .align_PWM_seq import alignPWMSeq

def computeYAndMask(pwm, dna_seq):
    seq_starts = [None, None]
    pwm_starts = [None,None]
    ks = [None, None]
    scores = [None, None]

    pwm_starts[0], seq_starts[0], ks[0], scores[0] = alignPWMSeq(pwm, dna_seq[0,:,:])
    pwm_starts[1], seq_starts[1], ks[1], scores[1] = alignPWMSeq(pwm, dna_seq[1,:,:])

    pwm_mask = [None, None] #aligned region of PWM
    dna_mask = [None, None] #aligned refion on DNA

    if scores[0] > scores[1]: #tracks which strand pwm aligned to
        primary_index = 0
    else:
        primary_index = 1

    #flip the pwm for the complementary strand both direction and complementarity
    Y_pwm = [pwm, pwm]
    Y_pwm[1 - primary_index] = np.flip(pwm)
    pwm_mask[primary_index] = np.zeros(len(pwm)).astype(bool)
    pwm_mask[primary_index][pwm_starts[primary_index]:pwm_starts[primary_index]+ks[primary_index]] = True
    pwm_mask[1-primary_index] = np.flip(pwm_mask[primary_index])

    dna_mask[primary_index] = np.zeros(dna_seq.shape[1]).astype(bool) 
    dna_mask[primary_index][seq_starts[primary_index]:seq_starts[primary_index]+ks[primary_index]] = True
    dna_mask[1-primary_index] = np.flip(dna_mask[primary_index])

    return Y_pwm, pwm_mask, dna_mask, scores[primary_index]/ks[primary_index]

