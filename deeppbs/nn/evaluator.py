# built-in modules
import warnings

# third party modules
import torch
import numpy as np

# geobind modules
from deeppbs.nn import processBatch
from deeppbs.nn.metrics import auroc, auprc, balanced_accuracy_score, recall_score, brier_score_loss, specificity
from deeppbs.nn.metrics import precision_score, jaccard_score, f1_score, accuracy_score, matthews_corrcoef
from deeppbs.nn.metrics import reportMetrics, chooseBinaryThreshold
from deeppbs.nn.metrics import IC_weighted_PCC, IC_corr, brier_multi, IC_diff, mae

from scipy.stats import spearmanr, pearsonr #, entropy

METRICS_FN = {
    'accuracy': accuracy_score,
    'balanced_accuracy': balanced_accuracy_score,
    'mean_iou': jaccard_score,
    'auroc': auroc,
    'auprc': auprc,
    'recall': recall_score,
    'precision': precision_score,
    'f1_score': f1_score,
    'brier_score': brier_score_loss,
    'matthews_corrcoef': matthews_corrcoef,
    'specificity': specificity,
    'pearsonr': pearsonr,
    'spearmanr': spearmanr,
    'ic_weighted_pcc': IC_weighted_PCC,
    'ic_corr': IC_corr,
    'brier_multi': brier_multi,
    'ic_diff':IC_diff,
    'mae': mae
}

def registerMetric(name, fn):
    METRICS_FN[name] = fn

class Evaluator(object):
    def __init__(self, model, nc, device="cpu", metrics=None, post_process=None, negative_class=0, labels=None):
        self.model = model # must implement the 'forward' method
        self.device = device
        self.negative_class = negative_class
        
        if post_process is None:
            # identity function
            post_process = lambda x: x
        self.post = post_process
        
        # decide what metrics to use
        self.nc = nc
        if metrics == 'none':
            self.metrics = None
        elif metrics is None:
            if nc == 2:
                # binary classifier
                metrics={
                    'auroc': {'average': 'binary'},
                    'auprc': {'average': 'binary'},
                    'balanced_accuracy': {},
                    'mean_iou': {'average': 'weighted'},
                    'precision': {'average': 'binary', 'zero_division': 0},
                    'recall': {'average': 'binary', 'zero_division': 0},
                    'accuracy': {},
                    'specificity': {},
                    'matthews_corrcoef': {}
                }
                metrics_check = {
                    'auroc': lambda n: (n[0] > 0) and (n[1] > 0),
                    'auprc': lambda n: (n[0] > 0) and (n[1] > 0),
                    'balanced_accuracy': lambda n: (n[0] > 0) and (n[1] > 0),
                    'mean_iou': lambda n: (n[0] > 0) and (n[1] > 0),
                    'precision': lambda n: (n[0] > 0) and (n[1] > 0),
                    'recall': lambda n: (n[1] > 0),
                    'accuracy': lambda n: True,
                    'specificity': lambda n: (n[0] > 0),
                    'matthews_corrcoef': lambda n: (n[0] > 0) and (n[1] > 0)
                }
            elif nc > 2:
                # three or more classes 
                if labels is None:
                    labels = list(range(nc))
                    labels.remove(negative_class)
                metrics={
                    'pearsonr': {},
                    'spearmanr':{},
                    'auroc': {'average': 'macro', 'multi_class' : 'ovo'},
                    'ic_weighted_pcc':{},
                    'ic_corr':{},
                    'brier_multi':{},
                    'ic_diff':{},
                    'mae':{}
                }
                metrics_check = {
                    'auroc': None,
                    'auprc': None,
                    'balanced_accuracy': None,
                    'mean_iou': None,
                    'precision': None,
                    'recall': None,
                    'accuracy': None,
                    'specificity': None
                }
            else:
                # TODO - assume regression
                pass
            self.metrics = metrics
            self.metrics_check = metrics_check
        else:
            if not isinstance(metrics, dict):
                raise ValueError("The argument 'metrics' must be a dictionary of kwargs and metric names or 'none'!")
            self.metrics = metrics
    
    @torch.no_grad()
    def eval(self, dataset, 
            eval_mode=True,
            batchwise=False,
            use_mask=True,
            return_masks=False,
            return_predicted=False,
            return_batches=True,
            xtras=None,
            split_batches=False,
            **kwargs
        ):
        """Returns numpy arrays!!!"""
        
        def _loop(batch, data_items, y_gts, y_prs, outps, masks, indexes, out_masks, batches):
            batch_data = processBatch(self.device, batch, xtras=xtras)
            batch, y_pwm0, y_pwm1, pwm_mask0, pwm_mask1, out_mask0, out_mask1 = batch_data['batch'], batch_data['y_pwm0'], batch_data['y_pwm1'], batch_data['pwm_mask0'], batch_data['pwm_mask1'], batch_data['dna_mask0'], batch_data['dna_mask1']

            
            y = torch.cat([y_pwm0, y_pwm1], dim=0)
            y_mask = torch.cat([pwm_mask0, pwm_mask1], dim=0)
            out_mask = torch.cat([out_mask0, out_mask1], dim=0)
            index = torch.cat([torch.ones_like(pwm_mask0),torch.zeros_like(pwm_mask0)])
            out_index = torch.cat([torch.ones_like(out_mask0),torch.zeros_like(out_mask0)])
            output = self.model(batch)
            if use_mask:
                y = y[y_mask].cpu().numpy()
                out = self.post(output[out_mask]).cpu().numpy()
            else:
                y = y.cpu().numpy()
                out = self.post(output).cpu().numpy()
            
            y_gts.append(y)
            logits.append(output.cpu().numpy())
            outps.append(out)
            if return_masks:
                masks.append(y_mask.cpu().numpy())
                out_masks.append(out_mask.cpu().numpy())
                indexes.append(index.cpu().numpy())
                out_indexes.append(out_index.cpu().numpy())
                
            if return_predicted:
                y_prs.append(self.predictClass(out, y, **kwargs))
            
            if xtras is not None:
                # these items will not be masked even if `use_mask == True`
                for item in xtras:
                    data_items[item].append(batch_data[item].cpu().numpy())
            
            if return_batches:
                if isinstance(batch, list):
                    if batchwise:
                        batches.append(batch)
                    else:    
                        batches += batch
                else:
                    if batchwise:
                        batches.append([batch])
                    else:
                        batches.append(batch.to('cpu'))
        
        # eval or training
        if eval_mode:
            self.model.eval()
        
        # evaluate model on given dataset
        data_items = {}
        y_gts = []
        y_prs = []
        outps = []
        y_masks = []
        out_masks = []
        indexes = []
        batches = []
        out_indexes = []
        logits = []
        if xtras is not None:
            # these items will not be masked even if `use_mask == True`
            for item in xtras:
                data_items[item] = []
        
        # loop over dataset
        for batch in dataset:
            if split_batches:
                dl = batch.to_data_list()
                for d in dl:
                    _loop(d, data_items, y_gts, y_prs, outps, y_masks, indexes, out_masks, batches)
            else:
                _loop(batch, data_items, y_gts, y_prs, outps, y_masks, indexes, out_masks, batches)
        
        # decide what to do with each data item
        data_items['y'] = y_gts
        data_items['output'] = outps

        if return_masks:
            data_items['masks'] = y_masks
            data_items['out_masks'] = out_masks
            data_items['indexes'] = indexes
            data_items['out_idx'] = out_indexes
            data_items['logits'] = logits
        if return_predicted:
            data_items['predicted_y'] = y_prs
        
        # concat batches if not batchwise
        if batchwise:
            data_items['num_batches'] = len(y_gts)
        else:
            for item in data_items:
                data_items[item] = np.concatenate(data_items[item], axis=0)
            data_items['num'] = len(data_items['y'])
        
        # add batches if requested
        if return_batches:
            data_items['batches'] = batches
        
        return data_items
    
    def getMetrics(self, *args, 
            eval_mode=True, 
            metric_values=None, 
            threshold=None, 
            threshold_metric='balanced_accuracy', 
            report_threshold=False,
            metrics_calculation="total",
            split_batches=False,
            use_mask=True,
            label_type="vertex",
            **kwargs
        ):
        if self.metrics is None:
            return {}
        if metric_values is None:
            metric_values = {key: [] for key in self.metrics}
        if 'threshold' in metric_values:
            threshold = metric_values['threshold']
        
        if metrics_calculation == "total":
            batchwise = False
        elif metrics_calculation == "average_batches":
            batchwise = True
        else:
            raise ValueError("Invalid option for `metrics_calculation`: {}".format(metrics_calculation))
        
        # Determine what we were given (a dataset or labels/predictions)
        if len(args) == 1:
            evald = self.eval(args[0], eval_mode=eval_mode, use_mask=False, return_masks=True, return_batches=True, batchwise=batchwise, split_batches=split_batches, **kwargs)
            if batchwise:
                y_gt = evald['y']
                outs = evald['output']
                batches = evald['batches']
                masks = evald['masks']
                out_masks = evald['out_masks']
            else:
                y_gt = [evald['y']]
                outs = [evald['output']]
                batches = [evald['batches']]
                masks = [evald['masks']]
                out_masks = [evald['out_masks']]
        else:
            if len(args) == 4:
                y_gt, outs, masks, out_masks = args
                batches = None
            else:
                y_gt, outs, masks, out_masks, batches = args
                batches = [batches]
            y_gt = [y_gt]
            outs = [outs]
            masks = [masks]
            out_masks = [out_masks]

        # Get predicted class labels
        nan = float('nan')
        for i in range(len(y_gt)):
            if label_type == "graph":
                y_gt[i] = y_gt[i].flatten()
            
            
            if masks is not None and use_mask:    
                y_gt[i] = y_gt[i][masks[i]]
                outs[i] = outs[i][out_masks[i]]
            
            # Compute metrics
            for metric, kw in self.metrics.items():
                if metric == 'auprc' or metric == 'auroc':
                    # AUC metrics
                    metric_values[metric].append(METRICS_FN[metric](np.argmax(y_gt[i], axis=1), outs[i], **kw))
                elif metric == 'ic_weighted_pcc':
                    ic_pccs = []
                    weights = []
                    for index in range(outs[i].shape[0]):
                        ic_pcc, weight = METRICS_FN[metric](y_gt[i][index,:], outs[i][index,:], **kw)
                        ic_pccs.append(ic_pcc)
                        weights.append(weight)
                        metric_values[metric].append(np.average(ic_pccs, weights=weights))
                elif metric in ["ic_corr", "brier_multi", "ic_diff", "mae"]:
                    metric_values[metric].append(float(METRICS_FN[metric](y_gt[i], outs[i])))
                elif metric == "matthews_corrcoef":
                    continue
                elif metric in ['pearsonr','spearmanr']:
                    columnwise = []
                    IC_weights = []
                    for index in range(outs[i].shape[0]):
                        columnwise.append(METRICS_FN[metric](y_gt[i][index,:], outs[i][index,:], **kw)[0])
                    
                    metric_values[metric].append(np.average(columnwise))
                else:
                    #if self.metrics_check[metric](ngt):
                    metric_values[metric].append(METRICS_FN[metric](y_gt[i], y_pr, **kw))
        with warnings.catch_warnings():
            # ignore empty-slice warnings from numpy
            warnings.simplefilter("ignore")
            for key in metric_values:
                metric_values[key] = np.nanmean(metric_values[key])
        return metric_values
    
    def getGraphMetrics(self, batches, y_gt, y_pr, metric_values=None):
        if self.metrics is None:
            return {}
        if metric_values is None:
            metric_values = {}
        
        smooth_gt = []
        smooth_pr = []
        ptr = 0
        if not isinstance(batches, list):
            batches = [batches]
        for batch in batches:
            edge_index = batch.edge_index.cpu().numpy()
            slc = slice(ptr, ptr + batch.num_nodes)
            if "smoothness" in self.metrics:
                smooth_gt.append(METRICS_FN["smoothness"](y_gt[slc], edge_index, **self.metrics['smoothness']))
                smooth_pr.append(METRICS_FN["smoothness"](y_pr[slc], edge_index, **self.metrics['smoothness']))
            ptr += batch.num_nodes
        
        if "smoothness" in self.metrics:
            smooth_pr = np.mean(smooth_pr)
            smooth_gt = np.mean(smooth_gt)
            metric_values['smoothness'].append(smooth_pr)
            metric_values['smoothness_relative'].append((smooth_pr - smooth_gt)/smooth_gt)
        
        return metric_values
    
    def predictClass(self, outs, y_gt=None, metrics_dict=None, threshold=None, threshold_metric='balanced_accuracy', report_threshold=False):
        # Decide how to determine `y_pr` from `outs`
        if self.nc == 2:
            if (threshold is None) and (y_gt is not None):
                # sample n_samples threshold values
                threshold, _  = chooseBinaryThreshold(y_gt, outs[:,1], metric_fn=METRICS_FN[threshold_metric], **self.metrics[threshold_metric])
            elif (threshold is None) and (y_gt is None):
                threshold = 0.5
            y_pr = (outs[:,1] >= threshold)
            if report_threshold and (metrics_dict is not None):
                metrics_dict['threshold'] = threshold
        elif self.nc > 2:
            y_pr = np.argmax(outs, axis=1).flatten()
        
        return y_pr
