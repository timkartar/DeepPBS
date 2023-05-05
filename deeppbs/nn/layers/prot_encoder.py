import torch
import torch.nn as nn
import torch.nn.functional as F
from deeppbs.nn import MLP

from torch_geometric.nn import CGConv
from torch_cluster import radius_graph



class ProtEncoder(nn.Module):

    def __init__(self, in_channels, hidden_channels, condition="full"):
        super(ProtEncoder, self).__init__()

        self.condition = condition
        self.embed = MLP([in_channels, hidden_channels, hidden_channels])
        self.bond_conv1 = CGConv(hidden_channels, aggr='mean')
        self.bond_conv2 = CGConv(hidden_channels, aggr='mean')

        self.act = nn.ReLU()

        self.radi_conv1 = CGConv(hidden_channels, aggr='mean')
        self.radi_conv2 = CGConv(hidden_channels, aggr='mean')

        self.out_nn = MLP([hidden_channels*2, hidden_channels, hidden_channels])

    def forward(self, x_prot, v_prot, edge_index):
        x_prot = self.embed(x_prot)
        x_prot_bond = self.act(self.bond_conv1(x_prot, edge_index))
        x_prot_bond = self.act(self.bond_conv2(x_prot_bond, edge_index))


        r_edge_index = radius_graph(v_prot, r=4)
        x_prot_radi = self.act(self.radi_conv1(x_prot_bond, r_edge_index))
        x_prot_radi = self.act(self.radi_conv2(x_prot_radi, r_edge_index))

        x_prot = self.act(x_prot_radi)

        return x_prot

