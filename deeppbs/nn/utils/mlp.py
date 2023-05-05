from torch.nn import ReLU, ELU, Identity, Tanh, Dropout, PReLU
from torch.nn import Sequential as Seq, Linear as Lin, BatchNorm1d as BN

def MLP(channels, batch_norm=True, act='relu', bn_kwargs={}, dropout=0.0, dropout_position="right",
        batchnorm_position="right", bias=True):
    
    if batch_norm:
        # bad practice to mix these two
        assert dropout == 0.0
    if dropout > 0.0:
        # bad practice to mix these two
        assert not batch_norm
    
    if isinstance(act, str) or act is None:
        act = [act]*(len(channels)-1)
    
    if isinstance(batch_norm, bool):
        batch_norm = [batch_norm]*(len(channels)-1)
    
    if isinstance(dropout, float):
        dropout = [dropout]*(len(channels)-1)
    
    activation = []
    for a in act:
        if a is None:
            activation.append(Identity)
        elif a == 'relu':
            activation.append(ReLU)
        elif a == 'elu':
            activation.append(ELU)
        elif a == 'tanh':
            activation.append(Tanh)
        elif a == "prelu":
            activation.append(PReLU)
        else:
            raise ValueError("unrecognized keyword: {}".format(act))
    
    layers = []
    for i in range(1, len(channels)):
        if dropout[i-1] > 0 and dropout_position == "left":
            layers.append(Dropout(p=dropout[i-1]))
        
        layers.append(Lin(channels[i-1], channels[i], bias=bias))
        
        if batch_norm[i-1] and batchnorm_position == "left":
            layers.append(BN(channels[i], **bn_kwargs))
        
        layers.append(activation[i-1]())
        
        if batch_norm[i-1] and batchnorm_position == "right":
            layers.append(BN(channels[i], **bn_kwargs))
        
        if dropout[i-1] > 0 and dropout_position == "right":
            layers.append(Dropout(p=dropout[i-1]))
    
    return Seq(*layers)
