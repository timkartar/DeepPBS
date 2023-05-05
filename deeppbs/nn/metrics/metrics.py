import numpy as np
from sklearn.metrics import roc_curve, precision_recall_curve, auc
from sklearn.metrics import balanced_accuracy_score, recall_score, precision_score, jaccard_score
from sklearn.metrics import f1_score, accuracy_score, brier_score_loss, matthews_corrcoef
from sklearn.metrics import average_precision_score, roc_auc_score
from scipy.stats import pearsonr, entropy

def auprc(y_gt, prob, pi=1, average='binary', **kwargs):
    # determine if binary or multiclass
    if average == 'binary':
        pre_vals, rec_vals, _ = precision_recall_curve(y_gt, prob[:,pi])
        auprc = auc(rec_vals, pre_vals)
    else:
        nc = prob.shape[1]
        auprc = average_precision_score(np.eye(nc)[y_gt], prob, average=average, **kwargs)
        
    return auprc
'''
def auroc(y_gt, prob, pi=1, average='binary', **kwargs):
    # determine if binary or multiclass
    if average == 'binary':
        fpr, tpr, _ = roc_curve(y_gt, prob[:,pi])
        auroc = auc(fpr, tpr)
    else:
        #nc = prob.shape[1]
        #auroc = roc_auc_score(np.eye(nc)[y_gt], prob, average=average, **kwargs)
        auroc = roc_auc_score(y_gt, prob, average=average, **kwargs)
    
    return auroc
'''

def auroc(y_gt, prob, pi=1, average='binary', **kwargs):
    # determine if binary or multiclass
    if average == 'binary':
        fpr, tpr, _ = roc_curve(y_gt, prob[:,pi])
        auroc = auc(fpr, tpr)
    else:
        nc = prob.shape[1]
        one_hot = np.eye(nc)[y_gt]
        mask = (np.sum(one_hot, axis = 0) != 0)
        rescaled_prob = prob[:,mask] / prob[:,mask].sum(axis=1)[:,None]

        try:
            present_class = [i for i in range(nc) if i in np.unique(y_gt)]
            mapp = dict()
            for i in present_class:
                mapp[i] = present_class.index(i)
            y_gt = [mapp[i] for i in y_gt]
            if(len(present_class) <= 2):
                fpr, tpr, _ = roc_curve(y_gt, rescaled_prob[:,pi])
                auroc = auc(fpr, tpr)
            else:
                auroc = roc_auc_score(y_gt, rescaled_prob, average=average, **kwargs)
        except Exception as e:
            auroc = np.nan
            #print("helloooo", e, np.unique(y_gt), prob)#, np.eye(nc)[y_gt])
    return auroc

def specificity(ygt, ypr):
    ni = (ygt == 0)
    tn = (ypr[ni] == 0).sum()
    
    return tn/(ni.sum())

def IC_weighted_PCC(true_column, pred_column, gt=None):
    if( gt is None ):
        gt = true_column
    pcc = pearsonr(true_column, pred_column)[0]
    ic_weight = entropy(gt, [0.25,0.25,0.25,0.25], base=2)/2

    return pcc * ic_weight, ic_weight

def IC_corr(true, pred, bkg=[0.25,0.25,0.25,0.25]):
    ic_true = entropy(true, bkg, base=2, axis=1)
    ic_pred = entropy(pred, bkg, base=2, axis=1)

    return pearsonr(ic_true, ic_pred)[0]

def IC_diff(pwm, pred, bkg=[0.25,0.25,0.25,0.25]):
    columnwise = -1*np.abs(entropy(pwm,bkg, base=2, axis=1) - entropy(pred,bkg, base = 2, axis=1))
    return columnwise.mean()

def brier_multi(targets, probs, root=False):
    if not root:
        return np.mean(np.sum((probs - targets)**2, axis=1))
    else:
        return np.sqrt(np.mean(np.sum((probs - targets)**2, axis=1)))

def mae(targets, probs):
    return np.mean(np.sum(np.abs(probs - targets), axis=1))
