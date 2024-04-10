import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import sys
from deeppbs.nn import MLP

from torch_geometric.nn import CGConv, DynamicEdgeConv
from torch_cluster import radius, radius_graph
from deeppbs.nn import ProtEncoder, BiNet

class CNN(nn.Module):
    def __init__(self, dna_channels, hidden_size=8, condition="full"): 
        super(CNN, self).__init__()
        self.condition = condition
        self.conv1 = nn.Conv1d(dna_channels, hidden_size, kernel_size=3, padding='same') 
        self.conv2 = nn.Conv1d(hidden_size, hidden_size, kernel_size = 3, padding='same') 
        self.fc = nn.Conv1d(hidden_size, hidden_size, kernel_size=1)
        self.act = nn.ReLU() 

    def forward(self, x):
        x = self.conv1(x.T)
        x = self.act(x)  
        x = self.conv2(x)  
        x = self.act(x)   
        x = self.fc(x) 
        
        return x.T  

class Model(nn.Module):

    def __init__(self, prot_channels, dna_channels, out_channels=4, condition="full", readout="all", **kwargs):
        super(Model, self).__init__()
        
        self.condition = condition
        assert self.condition in ["prot","prot_ag","prot_shape_ag","shape_ag","shape","prot_shape","ag"]

        if self.condition in ["prot_ag","shape_ag","prot_shape_ag","ag"]:
            self.fn_channels = 6
        else:
            self.fn_channels = 0
        
        self.name = "Model"
        self.hidden_size = 8
        self.dna_channels = dna_channels
        self.binet_reduce_channels = 32
        self.dna_embed_dim = 10
        self.dropout = 0.0
        self.prot_embed_dim = 10
        self.conv = None #for visualizing conv output

        self.embed = MLP([11 + self.fn_channels, self.dna_embed_dim, self.dna_embed_dim],
                dropout=self.dropout)
        
        self.prot_encoder = ProtEncoder(prot_channels, self.prot_embed_dim, condition=self.condition)

        self.binet = BiNet(self.prot_embed_dim, self.dna_embed_dim, condition=self.condition,
        conv="PPFConv", readout=readout) # +dna_channels
        
        self.reduce_nn = MLP([11*self.dna_embed_dim, 2*self.binet_reduce_channels,
            self.binet_reduce_channels], dropout=self.dropout)

        if self.condition in ["prot_ag", "prot", "ag"]:
            self.cnn = CNN(self.binet_reduce_channels, hidden_size=self.hidden_size,
                    condition=self.condition)
        elif self.condition in ["prot_shape","prot_shape_ag","shape_ag"]:
            self.cnn = CNN(self.binet_reduce_channels +  self.dna_channels,
                    hidden_size=self.hidden_size, condition=self.condition)
        #elif self.condition == "shape":
        self.shapecnn = CNN(self.dna_channels, hidden_size=self.hidden_size, condition=self.condition)
        
        self.mlp = MLP([self.hidden_size, self.hidden_size, self.hidden_size, 4],
                dropout=self.dropout)
       
        #self.groove_cnn = CNN(self.dna_embed_dim, hidden_size=4, condition=self.condition,
        #            kernel_size=7)
        #self.bbone_cnn = CNN(self.dna_embed_dim, hidden_size=4, condition=self.condition,
        #        kernel_size=4)


        #self.bbone_shape_cnn = CNN(self.dna_channels + 4, hidden_size=4, condition=self.condition,
        #        kernel_size=3, padding='same')


        #self.bbone_shape_cnn2 = CNN(4, hidden_size=4, condition=self.condition,
        #        kernel_size=3, padding='same')

        self.global_temp = nn.Parameter(torch.randn(1))
        self.device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
        self.outnn = MLP([8,4,4])

    def strandForward(self, e_prot, v_dna, x_dna, x_dna_point, 
            x_prot, v_prot, prot_vecs, dna_vecs, atom_to_mask):
        e_prot = e_prot.T
        v_dna = v_dna.view(-1,3)
        dna_vecs = dna_vecs.view(-1,3)
        x_dna_copy = x_dna
        
        if self.condition != "shape":
            x_dna_point_embed = self.embed(x_dna_point)

            x_prot = self.prot_encoder(x_prot, v_prot, e_prot)
            x_binet, conv = self.binet(x_dna_point_embed, v_dna, x_prot, v_prot, prot_vecs,
                    dna_vecs, add_target_features=True, atom_to_mask=atom_to_mask) 
            
            #x_binet = x_binet.reshape(x_binet.shape[0], -1)
            
            x_binet = x_binet.reshape(x_binet.shape[0], -1)
            x_binet = self.reduce_nn(x_binet)

            #x_groove = x_binet[:,[2,3,4,5,6,7,8],:]
            #x_bbone = x_binet[:,[0,1,9,10],:]
            
            #x_groove = self.groove_cnn(x_groove)
            #x_binet = self.reduce_nn(x_binet)#.reshape(-1, 11)
            
            #x_bbone = self.bbone_cnn(x_bbone)
            if self.condition in ["shape_ag", "prot_shape","prot_shape_ag"]:
                #x_bbone = torch.hstack((x_bbone, x_dna))
                x_binet = torch.hstack((x_binet, x_dna)) 
                #x_binet = torch.hstack((x_binet, torch.zeros_like(x_dna)))
            
            #x_bbone = x_bbone.unsqueeze(0)
            #x_bbone = self.bbone_shape_cnn(x_bbone).T
            #x_bbone = self.bbone_shape_cnn2(x_bbone.unsqueeze(0)).T
            
            
            #x_dnacnn = torch.hstack((x_groove, x_bbone))
            #x_dna = self.outnn(x_dnacnn)

            x_dnacnn = self.cnn(x_binet)
            x_dna = self.mlp(x_dnacnn)
        else:
            x_dnacnn = self.shapecnn(x_dna)
        
            x_dna = self.mlp(x_dnacnn)
        
        #try:
        #    return x_dna*rel_temp[:,None] #torch.cat((x_dna, torch.flip(x_dna, dims=[0,1])), dim=0)/torch.sigmoid(self.global_temp)
        #except:
        return x_dna, conv

    def forward(self, data, atom_to_mask):
        x_dna_point = data.x_dna_point[:,:(11 + self.fn_channels)]
        out1, conv = self.strandForward(data.e_prot, data.v_dna, data.x_dna, x_dna_point, data.x_prot,
                data.v_prot, data.prot_vecs, data.dna_vecs, atom_to_mask)
        
        self.conv = conv

        x_dna_point = x_dna_point.reshape(-1, 11, 11 + self.fn_channels)
        
        switch = torch.LongTensor([10,9,5,4,3,2,8,7,6,1,0]).to(self.device)
        x_dna_point_rc = torch.index_select(x_dna_point, 1, switch)
        v_dna_rc = torch.index_select(data.v_dna, 1, switch)
        dna_vecs_rc = torch.index_select(data.dna_vecs, 1, switch)

        x_dna_point_rc[:,:,:11] = x_dna_point[:,:,:11]
        x_dna_point_rc = torch.flip(x_dna_point_rc, [0]).reshape(-1,11 + self.fn_channels)
        v_dna_rc = torch.flip(v_dna_rc, [0])
        dna_vecs_rc = torch.flip(dna_vecs_rc, [0])
        
        shape_transform = torch.LongTensor([-1, -1, 1, 1, 1, 1, -1, 1, 1, -1, 1, 1, 1,
            1]).to(self.device)

        x_dna_rc = torch.flip(data.x_dna, [0])*shape_transform[None, :]
        tmp = x_dna_rc[0,6:12].clone().detach()
        x_dna_rc[:-1,6:12] = x_dna_rc[1:,6:12]
        x_dna_rc[-1,6:12] = tmp

        out2, _  = self.strandForward(data.e_prot, v_dna_rc, x_dna_rc, x_dna_point_rc, data.x_prot,
                data.v_prot, data.prot_vecs, dna_vecs_rc, atom_to_mask)

        #out = (out1 + out2.flip([0,1]))/2
        
        return torch.cat((out1, out2),
                dim=0)/torch.sigmoid(self.global_temp)
        
        #return torch.cat((out, out.flip([0,1])),
        #        dim=0)/torch.sigmoid(self.global_temp)
