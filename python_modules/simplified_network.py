import torch
import torch.nn as nn
import torch.nn.functional as F

class SimpleFeatureQNet(nn.Module):
    """简化的特征网络，直接使用sheet特征而非完整图结构"""
    def __init__(self, feature_dim=10, hidden_dims=[64, 32]):
        super().__init__()
        
        # 构建多层感知机
        layers = []
        input_dim = feature_dim
        for hidden_dim in hidden_dims:
            layers.append(nn.Linear(input_dim, hidden_dim))
            layers.append(nn.ReLU())
            input_dim = hidden_dim
        layers.append(nn.Linear(input_dim, 1))
        
        self.mlp = nn.Sequential(*layers)
        
    def forward(self, x):
        """
        输入:
        x: [batch_size, num_sheets, feature_dim] 包含所有sheet特征的张量
        
        输出:
        q_values: [batch_size, num_sheets] 每个sheet的Q值
        """
        # 将x视为独立样本处理
        batch_size, num_sheets, feature_dim = x.shape
        x_flat = x.view(-1, feature_dim)
        
        # 计算每个sheet的Q值
        q_flat = self.mlp(x_flat)
        q_values = q_flat.view(batch_size, num_sheets)
        
        return q_values
