import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import sys
from deeppbs.nn import MLP

from torch_geometric.nn import CGConv
from .ppf_conv import PPFConv, point_pair_features
from torch_cluster import radius, radius_graph

def distance(p1, p2):
    return torch.sqrt(torch.sum(torch.square(p1 - p2), dim=1).unsqueeze(1))

def mask_edge_index(edge_index):
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

    p_indices = [0,10]# for phosphate
    s_indices = [1,9]# for phosphate
    M_indices = [2,3,4,5]# for major
    m_indices = [6,7,8]# for minor 
    
    p_edges = []
    s_edges = []
    M_edges = []
    m_edges = []

    for item in edge_index.T:
        item = item.data.cpu().numpy()
        mod = item[1]%11
        if mod in p_indices:
            p_edges.append(item)
        elif mod in s_indices:
            s_edges.append(item)
        elif mod in M_indices:
            M_edges.append(item)
        elif mod in m_indices:
            m_edges.append(item)
    #print(torch.Tensor(m_edges).to(device))
    #sys.exit()
    p_edges = torch.LongTensor(p_edges).to(device)
    s_edges = torch.LongTensor(s_edges).to(device)
    M_edges = torch.LongTensor(M_edges).to(device)
    m_edges = torch.LongTensor(m_edges).to(device)
    return p_edges.T, s_edges.T, M_edges.T, m_edges.T

class BiNet(nn.Module):

    def __init__(self, prot_channels, dna_channels, condition="full", conv="PPFConv", readout="all"):
        super(BiNet, self).__init__()

        self.condition = condition
        self.dna_channels = dna_channels
        self.prot_channels= prot_channels
        self.conv = conv

        self.readout = readout

        if self.conv == "PPFConv":
            local_nn1 = MLP([self.prot_channels + 4 + self.dna_channels, self.dna_channels, self.dna_channels],
                    batch_norm=False)
        
            self.conv1 = PPFConv(local_nn1, global_nn=None, aggr='sum',
                    add_self_loops=False)
        
            local_nn2 = MLP([self.prot_channels + 4 + self.dna_channels, self.dna_channels, self.dna_channels],
                    batch_norm=False)
        
            self.conv2 = PPFConv(local_nn2, global_nn=None, aggr='sum',
                    add_self_loops=False)

            local_nn3 = MLP([self.prot_channels + 4 + self.dna_channels, self.dna_channels, self.dna_channels],
                    batch_norm=False)
            self.conv3 = PPFConv(local_nn3, global_nn=None, aggr='sum',
                    add_self_loops=False)


            local_nn4 = MLP([self.prot_channels + 4 + self.dna_channels, self.dna_channels, self.dna_channels],
                    batch_norm=False)
        
            self.conv4 = PPFConv(local_nn4, global_nn=None, aggr='sum',
                    add_self_loops=False)
        
        elif self.conv == "CGConv":
            self.conv1 = CGConv((self.prot_channels, self.dna_channels), dim=4, aggr='sum')


        
    
    def forward(self, x_dna, v_dna, x_prot, v_prot, prot_vec, dna_vec, add_target_features,
            atom_to_mask):

        if self.condition in ["prot_shape","prot","prot_ag", "prot_shape_ag"]:
            edge_index = radius(v_prot, v_dna, 5).flip(dims=[0])
            if 'str' in str(type(atom_to_mask)):
                np.save(atom_to_mask + "_edge_index.npy",edge_index.data.cpu().numpy())
            elif atom_to_mask is not None:
                edge_index = edge_index[:,edge_index[0, :] != atom_to_mask]
            
            p_edges, s_edges, M_edges, m_edges = mask_edge_index(edge_index)
            
            #edge_weight = 1/distance(v_prot[edge_index[0]], v_dna[edge_index[1]])
            
            #print(v_dna.view(-1,11,3).shape, x_dna.view(-1,11,10).shape, dna_vec.view(-1,11,3).shape)
            ##TODO for different convoloutions on bacbone, major, minor grooves, 
            ## mask edges accordingly
            
            if self.conv == "PPFConv":
                if len(p_edges) > 0 and self.readout in ["all","shape"]:
                    outp, p_conv = self.conv1((x_prot, x_dna), (v_prot, v_dna), (prot_vec, dna_vec),
                        p_edges, add_target_features)
                    outp = F.relu(outp)
                else:
                    outp = torch.zeros_like(x_dna)
                    p_conv = torch.zeros_like(x_dna)

                if len(s_edges) > 0 and self.readout in ["all","shape"]:
                    outs, s_conv = self.conv2((x_prot, x_dna), (v_prot, v_dna), (prot_vec, dna_vec),
                        s_edges, add_target_features)
                    outs = F.relu(outs)
                else:
                    outs= torch.zeros_like(x_dna)
                    s_conv= torch.zeros_like(x_dna)
                
                if len(M_edges) > 0 and self.readout in ["all","base"]:
                    outM, M_conv = self.conv3((x_prot, x_dna), (v_prot, v_dna), (prot_vec, dna_vec),
                        M_edges, add_target_features)
                    outM = F.relu(outM)
                else:
                    outM = torch.zeros_like(x_dna)
                    M_conv = torch.zeros_like(x_dna)
                
                if len(m_edges) > 0 and self.readout in ["all","base"]:
                    outm, m_conv = self.conv4((x_prot, x_dna), (v_prot, v_dna), (prot_vec, dna_vec),
                        m_edges, add_target_features)
                    outm = F.relu(outm)
                else:
                    outm = torch.zeros_like(x_dna)
                    m_conv = torch.zeros_like(x_dna)
                out = outp + outs + outM + outm ## can sum because the target features are 0 and additive
                conv = p_conv + s_conv + M_conv + m_conv

            elif self.conv == "CGConv":
                #edge_weight = 1/distance(v_prot[edge_index[0]], v_dna[edge_index[1]])  
                edge_features = point_pair_features(v_dna[edge_index[1]], v_prot[edge_index[0]],
                    dna_vec[edge_index[1]], prot_vec[edge_index[0]])
                out = F.relu(self.conv1((x_prot, x_dna), edge_index, edge_features))
                conv = out
            #out = F.relu(self.conv1((x_prot, x_dna), edge_index, edge_weight))
            #edge_index = radius(v_prot, v_dna, 5).flip(dims=[0])
            #edge_weight = 1/distance(v_prot[edge_index[0]], v_dna[edge_index[1]])

            #rel_temp = F.relu(self.conv2((x_prot, torch.zeros_like(x_dna)), edge_index, edge_weight))
            #rel_temp = rel_temp.reshape(-1, 11* self.dna_channels).sum(dim=1)
            #rel_temp = 1 + rel_temp#/(torch.max(rel_temp))
            #print(outp.sum(), outM.sum(), outm.sum(), outs.sum())
            #out = torch.cat([outp, outs, outM, outm], dim=1)
            return out.view(-1, 11, self.dna_channels), conv

        elif self.condition in ["shape", "shape_ag","ag"]:
            out = x_dna
            return out.view(-1, 11, self.dna_channels), 1

