from copy import deepcopy
import random

import numpy as np
from torch_geometric.data import Batch

class BalancedDataLoader(object):
    def __init__(self, data_list, nc, batch_size, shuffle=True):
        
        self.data_list = data_list
        self.nc = nc
        self.batch_size = max(batch_size, nc-1)
        self.classes = list(range(1, nc))
        self.shuffle=shuffle
        
        # get class of each data object
        self.data_classes = [[] for _ in range(nc-1)]
        
        for i in range(len(data_list)):
            d = data_list[i]
            mask = (d.y >= 0)
            ys = np.argwhere(np.bincount(d.y[mask], minlength=3) > 0).flatten()
            for y in ys[1:]:
                self.data_classes[y-1].append(i)
    
    def __iter__(self):
        self._data_classes = deepcopy(self.data_classes)
        self.available = set()
        self.counter = 0
        for i in range(self.nc-1):
            if self.shuffle:
                random.shuffle(self._data_classes[i])
            self.available.update(self._data_classes[i])
        
        return self
    
    def __next__(self):
        datas = []
        for i in range(self.batch_size):
            ci = self.counter % (self.nc - 1)
            while True:
                if len(self._data_classes[ci]) == 0:
                    raise StopIteration
                
                di = self._data_classes[ci].pop()
                if di in self.available:
                    datas.append(self.data_list[di])
                    self.available.remove(di)
                    break
                else:
                    continue
            self.counter += 1
        
        return Batch.from_data_list(datas)
