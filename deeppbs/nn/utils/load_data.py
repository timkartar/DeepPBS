# builtin modules
import os.path as osp
import hashlib
from pickle import dump, load
import sys

# third party modules
import numpy as np
import torch
from torch_geometric.data import Data, InMemoryDataset
from sklearn.preprocessing import StandardScaler

# deeppbs modules
from deeppbs.nn.utils import balancedClassIndices

class NodeScaler(object):
    def __init__(self):
        self._data_arrays = []
        self.scaler = StandardScaler()
    
    def update(self, array):
        self._data_arrays.append(array)
    
    def fit(self):
        self.scaler.fit(np.concatenate(self._data_arrays, axis=0))
    
    def scale(self, array):
        return self.scaler.transform(array)

class ClassificationDatasetMemory(InMemoryDataset):
    def __init__(self, data_files, nc, labels_key, data_dir,
            save_dir=None,
            transform=None,
            pre_transform=None,
            pre_filter=None,
            balance='balanced',
            percentage=1.0,
            remove_mask=False,
            unmasked_class=0,
            scale=True,
            scaler=None,
            protein_scaler=None,
            feature_mask=None
        ):
        if(save_dir is None):
            save_dir = data_dir
        self.data_dir = data_dir
        self.save_dir = save_dir
        self.data_files = data_files
        self.labels_key = labels_key
        self.nc = nc
        self.balance = balance
        self.percentage = percentage
        self.remove_mask = remove_mask
        self.unmasked_class = unmasked_class
        self.scale = scale
        self.scaler = scaler
        self.protein_scaler = protein_scaler
        self.transform = transform
        self.pre_filter = pre_filter
        self.pre_transform = pre_transform
        self.feature_mask = feature_mask
        
        super(ClassificationDatasetMemory, self).__init__(save_dir, transform, pre_transform, pre_filter)
        # load data
        self.data, self.slices = torch.load(self.processed_paths[0])
        
        # load scaler
        if self.scale and self.scaler is None:
            self.scaler = load(open(self.processed_paths[1], 'rb'))
        
    @property
    def raw_file_names(self):
        return self.data_files
    
    @property
    def processed_file_names(self):
        m = hashlib.md5()
        args = [
            self.nc,
            self.labels_key,
            self.balance,
            self.percentage,
            self.remove_mask,
            self.unmasked_class,
            self.scale
        ]
        args = "".join([str(_) for _ in args] + list(sorted(self.data_files)))
        m.update(args.encode('utf-8'))
        self.hash_name = m.hexdigest()
        return ['{}.pt'.format(self.hash_name), '{}_scaler.pkl'.format(self.hash_name)]
    
    @property
    def raw_dir(self):
        return self.data_dir

    @property
    def processed_dir(self):
        return self.save_dir
    
    def process(self):
        # get datalist and scaler
        data_list, transforms = _processData(self.raw_paths, self.nc, self.labels_key,
            balance=self.balance, 
            remove_mask=self.remove_mask,
            unmasked_class=self.unmasked_class,
            scaler=self.scaler,
            scale=self.scale,
            pre_filter=self.pre_filter,
            pre_transform=self.pre_transform,
            transform=self.transform,
            feature_mask=self.feature_mask
        )
        
        # save data
        data, slices = self.collate(data_list)
        torch.save((data, slices), self.processed_paths[0])
        
        # save scaler
        scaler = transforms['scaler']
        if scaler is not None:
            self.scaler = scaler
            dump(scaler, open(self.processed_paths[1], "wb"))

def _processData(data_files, nc, labels_key, 
        balance="unmasked",
        remove_mask=False,
        unmasked_class=0,
        scaler=None,
        protein_scaler=None,
        scale=True,
        transform=None,
        pre_filter=None,
        pre_transform=None,
        feature_mask=None,
        label_type="vertex",
        use_seq=True
    ):

    data_list = []
    
    # read and process datafiles
    data_fs = []
    for f in data_files:
        data_arrays = np.load(f, allow_pickle=True)
        
        Y_pwm = data_arrays["Y_pwm"]
        pwm_mask = data_arrays["pwm_mask"]
        Y_hard = data_arrays["Y_hard"]  
        dna_mask = data_arrays["dna_mask"]
        V_prot = data_arrays["V_prot"] 
        X_prot = data_arrays["X_prot"]
        E_prot = data_arrays["E_prot"]
        prot_fnames = data_arrays["prot_feature_names"]
        V_dna = np.array(data_arrays["V_dna"])
        X_dna_point = np.array(data_arrays["X_dna_point"])
        X_dna = data_arrays["X_dna"]
        try:
            prot_vecs = data_arrays["prot_vectors"]
            dna_vecs = data_arrays["dna_vectors"]
        except:
            print(f)
            continue

        #if not use_seq:
        mask = [True]*11 + [True]*6 + [False]*X_dna.shape[1] #11 bead position, 6 chemical groups ADMHPS
        X_dna_point = (X_dna_point.T[mask]).T
        
        try:
            V_dna = np.array(V_dna).astype(float)
        except:
            print("Missing Major/Minor position:", f)
            continue
        dna_fnames = data_arrays["dna_feature_names"]
        
        data = Data(
            x_prot = torch.tensor(X_prot, dtype=torch.float32),
            v_prot = torch.tensor(V_prot, dtype=torch.float32),
            x_dna =  torch.tensor(X_dna, dtype=torch.float32),
            v_dna = torch.tensor(V_dna, dtype=torch.float32),
            x_dna_point = torch.tensor(X_dna_point, dtype=torch.float32),
            y_pwm0 = torch.tensor(Y_pwm[0], dtype=torch.float32),
            y_pwm1 = torch.tensor(Y_pwm[1], dtype=torch.float32),
            y_hard0 = torch.tensor(Y_hard[0], dtype=torch.float32),
            y_hard1 = torch.tensor(Y_hard[1], dtype=torch.float32),
            e_prot = torch.tensor(E_prot, dtype=torch.long),
            prot_vecs = torch.tensor(prot_vecs, dtype=torch.float32),
            dna_vecs = torch.tensor(dna_vecs, dtype=torch.float32),
            dna_mask0 = torch.tensor(dna_mask[0], dtype=torch.bool),
            dna_mask1 = torch.tensor(dna_mask[1], dtype=torch.bool),
            pwm_mask0 = torch.tensor(pwm_mask[0], dtype=torch.bool),
            pwm_mask1 = torch.tensor(pwm_mask[1], dtype=torch.bool)
            )
        data_list.append(data)
        
        data_fs.append(f.split("/")[-1])
    
    
    # filter data
    if pre_filter is not None:
        data_list = [data for data in data_list if pre_filter(data)]
    
    # transform data
    if pre_transform is not None:
        data_list = [pre_transform(data) for data in data_list]
    
    # scale data
    if scale:
        # build a scaler
        if scaler is None:
            scaler = NodeScaler()
            for data in data_list:
                scaler.update(data.x_dna)
            scaler.fit()
            scaler = scaler.scaler
        
        if protein_scaler is None:
            protein_scaler = NodeScaler()
            for data in data_list:
                protein_scaler.update(data.x_prot[:,4:]) 
                # except 4 atom type features on protein, ideally pass non-onehot feature indices here
             
            protein_scaler.fit()
            protein_scaler = protein_scaler.scaler

        # scale node features in each data object
        for data in data_list:
            data.x_dna = torch.tensor(scaler.transform(data.x_dna), dtype=torch.float32)
            data.x_prot[:,4:] = torch.tensor(protein_scaler.transform(data.x_prot[:,4:]), dtype=torch.float32)
    transforms = {
        "scaler": scaler,
        "protein_scaler": protein_scaler,
        "transform": transform,
        "pre_transform": pre_transform,
        "pre_filter": pre_filter
    }
    return data_list, transforms, data_fs

def loadDataset(data_files, nc, labels_key, data_dir, cache_dataset=False, **kwargs):
    if isinstance(data_files, str):
        with open(data_files) as FH:
            data_files = [_.strip() for _ in FH.readlines()]
    
    if cache_dataset:
        dataset = ClassificationDatasetMemory(data_files, nc, labels_key, data_dir, **kwargs)
        transforms = {
            "scaler": dataset.scaler,
            "protein_scaler": dataset.protein_scaler,
            "transform": dataset.transform,
            "pre_transform": dataset.pre_transform,
            "pre_filter": dataset.pre_filter
        }
        info = {
            "num_features": dataset.num_node_features,
            "num_classes": nc,
            "num_instances": len(dataset)
        }
    else:
        data_files = [osp.join(data_dir, f) for f in data_files]
        
        dataset, transforms, data_files = _processData(data_files, nc, labels_key, **kwargs)
        info = {
            "prot_features": int(dataset[0].x_prot.shape[1]),
            "dna_features": int(dataset[0].x_dna.shape[1]),
            "num_classes": nc,
            "num_instances": len(dataset)
        }
    
    return dataset, transforms, info, data_files
