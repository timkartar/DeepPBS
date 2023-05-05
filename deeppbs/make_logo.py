import numpy as np
from scipy.stats import entropy
import pandas as pd

def makeLogo(gt, ax):
    import logomaker as lm
    
    e = entropy(gt,[0.25,0.25,0.25,0.25], axis=1, base=2)
    gt = gt * e[:,np.newaxis]
    gt_dict = {'A': [],
        'C':[],
        'G':[],
        'T':[]
        }
    for item in gt:
        gt_dict['A'].append(item[0])
        gt_dict['C'].append(item[1])
        gt_dict['G'].append(item[2])
        gt_dict['T'].append(item[3])

    gt_df = pd.DataFrame(gt_dict)
    colors = {'A': 'green',
        'C': 'blue',
        'G': '#FBA922',#'black',
        'T': 'red'
        }
    logo = lm.Logo(gt_df, color_scheme=colors, ax=ax, show_spines=False, alpha=0.7, fade_probabilities=False, font_name='DejaVu Sans')
    return logo


