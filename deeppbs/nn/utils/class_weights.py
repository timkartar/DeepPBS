# third party modules
import torch

# deeppbs modules
from deeppbs.nn import processBatch

def classWeights(data, nc, device='cpu', use_mask=True):
    if isinstance(data, torch.Tensor):
        # a tensor of class labels
        weight = data.shape[0]/(nc*torch.eye(nc)[data].sum(axis=0))
    else:
        # a dataloader object
        ys = []
        for batch in data:
            batch_data = processBatch(device, batch)
            y, mask = batch_data['y'], batch_data['mask'] 
            if use_mask :
                y = y[mask]
            ys.append(y)
        
        ys = torch.cat(ys, axis=0)
        weight = ys.shape[0]/(nc*torch.eye(nc)[ys].sum(axis=0))
    
    return weight.to(device)
