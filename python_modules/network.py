import torch

import torch.nn as nn

class U_linear_QNet(nn.Module):
    def __init__(self, dim_state, dim_action, size1, size2, size3):
        super(U_linear_QNet, self).__init__()
        
        # Encoder path (up)
        self.up1 = nn.Linear(dim_state, size1)
        self.up2 = nn.Linear(size1, size2)
        self.up3 = nn.Linear(size2, size3)
        
        # Middle layer
        self.fc = nn.Linear(size3, size3)
        
        # Decoder path (down)
        self.down1 = nn.Linear(size3, size2)
        self.down2 = nn.Linear(size2, size1)
        # 修改：将输出层维度改为1，输出每行得分
        self.down3 = nn.Linear(size1, 1)
        
        # Skip connections
        # 修改：输出维度改为1
        self.skip1 = nn.Linear(dim_state, 1)
        self.skip2 = nn.Linear(size1, size1)
        self.skip3 = nn.Linear(size2, size2)
        
    def forward(self, x):
        # Skip connections
        skip1_tensor = self.skip1(x)
        u1 = torch.relu(self.up1(x))
        skip2_tensor = self.skip2(u1)
        u2 = torch.relu(self.up2(u1))
        skip3_tensor = self.skip3(u2)
        u3 = torch.relu(self.up3(u2))
        
        # Middle layer
        d0 = self.fc(u3)
        
        # Decoder path with skip connections
        d1 = torch.relu(self.down1(d0))
        d2 = torch.relu(self.down2(d1 + skip3_tensor))
        d3 = torch.relu(self.down3(d2 + skip2_tensor))
        
        return d3 + skip1_tensor
