import matplotlib.pyplot as plt
from matplotlib import colors, cm
import numpy as np
from Bio import motifs
import sys

def plotPWM(Z, ax, use_mask=False, mask=None, fontsize=10, cmaps=None, xaxis=False):
    #shape of Z : (4,N)
    C = []

    def yellow_map(x):
        target = np.array([0.984313725490196, 0.6627450980392157, 0.13333333333333333, 1])#'#FBA922'
        white = np.array([1,1,1,1])
        return white + x*(target - white)
        #return([1, 1, 1-x, 1])

    cmaps = [cm.get_cmap('Greens'),cm.get_cmap('Blues'),yellow_map,cm.get_cmap('Reds')]
    if cmaps is None:
        cmaps = [cm.get_cmap('Greens'),cm.get_cmap('Blues'),cm.get_cmap('Greys'),cm.get_cmap('Reds')]
    

    for i in range(Z.shape[0]):
        to_append = []
        for j in range(Z.shape[1]):
            if use_mask:
                if not mask[j]:
                    to_append.append([0,0,0,1])
                    continue
            to_append.append(cmaps[i](Z[i,j]))
        C.append(to_append)
    C = np.array(C)
    im = ax.imshow(C)
    if not xaxis:
        ax.set_xticks([])
    ax.set_yticks([0,1,2,3])
    ax.set_yticklabels(['A','C','G','T'], fontsize=fontsize)
    ax.tick_params(axis=u'both', which=u'both',length=0)
    return ax, im

def displayAlignedPWM(pwm, pwm_mask, dna, dna_mask, ax1, ax2, label1 = "Co-crystal Sequence", label2 = "PWM",
        show_mask1=False,
        show_mask2=True
        ): 
        '''
        pass everything as lists, True values in masks 
        should denote respective aligned regions
        ax1 for Seq
        ax2 for PWM
        '''
        dna_al_start = dna_mask.index(True)
        
        pwm_al_start = pwm_mask.index(True)

        dummy = [[0,0,0,0]]
        if(pwm_al_start > dna_al_start):
            for idx in range(pwm_al_start - dna_al_start):
                dna = dummy + dna
                dna_mask = [False] + dna_mask
        else:
            for idx in range(dna_al_start - pwm_al_start):
                pwm = dummy + pwm
                pwm_mask = [False] + pwm_mask
        
        dna_al_end = len(dna_mask) #- dna_mask[::-1].index(True) - 1
        pwm_al_end = len(pwm_mask) #- pwm_mask[::-1].index(True) - 1
        
        if(pwm_al_end > dna_al_end):
            for idx in range(pwm_al_end - dna_al_end):
                dna = dna + dummy
                dna_mask = dna_mask + [False]
        else:
            for idx in range(dna_al_end - pwm_al_end):
                pwm = pwm + dummy
                pwm_mask = pwm_mask + [False]


        pwm = np.array(pwm).T
        dna = np.array(dna).T

        ax1, _ = plotPWM(dna, ax1, use_mask=show_mask1, mask=dna_mask)
        ax1.set_title(label1)
        ax2, _ = plotPWM(pwm, ax2, use_mask=show_mask2, mask=pwm_mask)
        ax2.set_title(label2)
        
        return ax1, ax2, dna, pwm, dna_mask, pwm_mask

if __name__ == "__main__":
    data = dict(np.load("5eg6_C_SUH_MOUSE.H11MO.0.A.npz", allow_pickle=True))
    #m = motifs.read(open("./MA1153.1.jaspar","r"),'jaspar')
    pwm = data['Y_pwm'][0].tolist()
    dna = data['Y_hard'][0].tolist()
    pwm_mask = data['pwm_mask'][0].tolist()
    dna_mask = data['dna_mask'][0].tolist()
    
    
    fig, (ax1, ax2) = plt.subplots(2,1)
    

    ax1, ax2 = displayAlignedPWM(pwm, pwm_mask, dna, dna_mask, ax1, ax2)
    plt.tight_layout()
    plt.show()
