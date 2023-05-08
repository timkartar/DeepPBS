import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from deeppbs import displayAlignedPWM, makeLogo, plotPWM
from os.path import join as ospj
import os
from tqdm import tqdm
import sys
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

root = "config_01.04.2023.13.46_737"
folder = "/predictions/"
dirname = "../output/" + root + folder
try:
    outroot = sys.argv[1]
except:
    outroot= root
out_path = "../figs/" + outroot + folder

if not os.path.exists(out_path):
    os.makedirs(out_path)

for f in tqdm(os.listdir(dirname)):
    if "2geq" not in f:
        continue
    pred_data = dict(np.load(ospj(dirname,f), allow_pickle=True))

    f = "_".join(f.split("_")[:-1])
    gt_data = dict(np.load("../../../dataset/assembly/" + f, allow_pickle=True))
    
    #print(gt_data.keys(), pred_data.keys())
    #print(pred_data['P'].shape, pred_data['Y'].shape, pred_data['idx'], gt_data['Y_pwm'][0].shape)
    #print(pred_data['P'][pred_data['idx'] == 0].shape, gt_data['Y_pwm'][0].shape)

    pwm = gt_data['Y_pwm'][0].tolist()
    dna = gt_data['Y_hard'][0].tolist()
    pwm_mask = gt_data['pwm_mask'][0].tolist()
    dna_mask = gt_data['dna_mask'][0].tolist()
    
    aln_score = gt_data['aln_score'][0]
    
    fig, axes = plt.subplots(2,2, figsize=(20,11))

    ax1 = axes[0][0]
    ax2 = axes[1][0]

    ax5 = axes[0][1]
    ax6 = axes[1][1]


    #ax5 = axes[0][2]
    #ax6 = axes[1][2]
    
    #ax7 = axes[0][3]
    #ax8 = axes[1][3]

    ax1, ax2, _, _, _, _ = displayAlignedPWM(pwm, pwm_mask, dna, dna_mask, ax1, ax2,
            label1="Co-crystal Seq (aln_score: {:.3f})".format(aln_score), label2="PWM", show_mask1=False, show_mask2=True)
    
    pwm = pred_data['P'].tolist()
    dna = pred_data['Y'].tolist()
    pwm_rc = pred_data['P_rc'].tolist()

    pwm_mask = pred_data['P_mask'].tolist()
    pwm_rc_mask = pred_data['P_rc_mask'].tolist()

    dna_mask = pred_data['Y_mask'].tolist()

    gt = np.array(dna)[dna_mask]
    pred = np.array(pwm)[pwm_mask]
    pred_rc = np.array(pwm_rc)[pwm_rc_mask]
    pred = (pred + np.flip(pred_rc,[0,1]))/2
    
    #ax3, ax4, pwm_ret, pred_ret, mask_ret, _ = displayAlignedPWM(pwm, pwm_mask, dna, dna_mask, ax3, ax4, label1="GT", label2="Pred", show_mask1=True, show_mask2=True)
    
    #difference = np.abs(pwm_ret - pred_ret)
    #ax3, _ = plotPWM(difference, ax3, use_mask=True, mask=mask_ret, cmaps=[cm.get_cmap('plasma')]*4)
    
    #ax3.set_title("|PWM - Pred|")
    #pos = ax3.get_position()
    #cbar_pos = [pos.x0 + 0.3, pos.y0 + 0.3,  pos.width / 2.0, pos.height / 2.0] 
    #cax = fig.add_axes(pos)
    #divider = make_axes_locatable(ax3)
    #cax = divider.append_axes("bottom", size="20%", pad=0.05)

    #x_ticks = cax.get_xticklabels()
    #cax.set_xticklabels(x_ticks, rotation=0, fontsize=5)

    #mappable = cm.ScalarMappable(cmap='plasma')
    #fig.colorbar(mappable, cax=cax, orientation='horizontal', fraction=0.046, pad=0.04)

    '''
    gt = np.array(dna)[dna_mask]
    pred = np.array(pwm)[pwm_mask]
    pred_rc = np.array(pwm_rc)[pwm_rc_mask]
    '''
    logo = makeLogo(gt, ax5)
    pred_logo = makeLogo(pred, ax6)

    ax5.set_title("PWM")
    ax6.set_title("Pred")
    ax5.set_xlabel("position")
    ax5.set_ylabel("bits")
    ax6.set_xlabel("position")
    ax6.set_ylabel("bits")
    
    ax5.tick_params(axis='both', which='major', labelsize=15)
    ax6.tick_params(axis='both', which='major', labelsize=15)

    ax5.set_ylim(0,2)
    ax6.set_ylim(0,2)
    
    '''
    logo = makeLogo(np.flip(gt,[1]), ax7)
    pred_logo = makeLogo(np.flip(pred_rc,[0]), ax8)

    ax7.set_title("PWM")
    ax8.set_title("Pred")
    ax7.set_xlabel("position")
    ax7.set_ylabel("bits")
    ax8.set_xlabel("position")
    ax8.set_ylabel("bits")
    
    ax7.set_ylim(0,2)
    ax8.set_ylim(0,2)
    '''
    #plt.title("{} Alignment Score: {:.3f}".format(f, aln_score))
    plt.tight_layout()
    plt.show()
    sys.exit()
    plt.savefig(ospj(out_path, "{}.svg".format(f)), dpi=300)

    plt.close()
