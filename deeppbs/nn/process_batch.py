# third party modules
import torch

def processBatch(device, batch, xtras=None):
    keys = list(batch[0].to_dict().keys())
    batch_data = {}
    if isinstance(batch, list):
        for key in keys:
            batch_data[key] = torch.cat([data.to_dict()[key] for data in batch]).to(device)
        if xtras is not None:
            for item in xtras:
                batch_data[item] = torch.cat([getattr(data, item) for data in batch]).to(device)
        batch_data['batch'] = batch
    else:
        for key in keys:
            batch_data['batch'] = batch.to(device)
            batch_data[key] = batch.to_dict()[key]
    
        if xtras is not None:
            for item in xtras:
                batch_data[item] = getattr(batch, item)
    
    return batch_data
