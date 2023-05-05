import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import sys
from deeppbs.nn import MLP

class CNN(nn.Module):

    def __init__(self, dna_channels, hidden_size=8, condition="full", kernel_size=7, padding=0):
        super(CNN, self).__init__()

        self.condition = condition
        # 1D
        self.conv1 = nn.Conv1d(dna_channels, hidden_size, kernel_size=kernel_size, padding=padding)
        #self.conv2 = nn.Conv1d(hidden_size, hidden_size, kernel_size = 3, padding='same')
        #self.fc = nn.Conv1d(hidden_size, hidden_size, kernel_size=1)
        self.act = nn.ReLU()
    def forward(self, x):
        
        x  = x.permute(0,2,1)
        #1D
        x = self.conv1(x)
        x = self.act(x)

        #x = self.conv2(x)
        #x = self.act(x)

        #x = self.fc(x)


        return x.squeeze()

