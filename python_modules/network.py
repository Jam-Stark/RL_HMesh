import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, global_mean_pool, GATConv, GATv2Conv


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
        self.down3 = nn.Linear(size1, dim_action)
        
        # Skip connections
        self.skip1 = nn.Linear(dim_state, dim_action)
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

class HexMeshQNet(nn.Module):
    def __init__(self, node_feat_dim=10, hidden_dim=64, num_gnn_layers=3):
        """
        初始化六面体网格的Q网络
        
        参数:
        node_feat_dim: 节点特征维度，默认为10对应新的State特征维度
                       [energy, boundary_ratio, feature_ratio, endpoints_special, adjacent_features,
                        curvature_metric, normal_variation, dihedral_deviation, min_jacobian, avg_jacobian]
        hidden_dim: GNN隐藏层维度
        num_gnn_layers: GNN层数
        """
        super().__init__()
        # GNN骨干网络
        self.convs = nn.ModuleList()
        self.convs.append(GCNConv(node_feat_dim, hidden_dim))
        for _ in range(num_gnn_layers-1):
            self.convs.append(GCNConv(hidden_dim, hidden_dim))
        
        # 几何特征编码器 - 新增几何特征的编码层
        self.geometry_encoder = nn.Sequential(
            nn.Linear(node_feat_dim, hidden_dim),
            nn.ReLU(inplace=True),
            nn.Linear(hidden_dim, hidden_dim)
        )
        
        # 融合几何特征和图特征的MLP
        self.fusion_layer = nn.Sequential(
            nn.Linear(hidden_dim*2, hidden_dim),
            nn.ReLU(inplace=True)
        )
        
        # MLP将[sheet_emb | global_emb]映射到Q值
        self.q_mlp = nn.Sequential(
            nn.Linear(hidden_dim*2, hidden_dim),
            nn.ReLU(inplace=True),
            nn.Linear(hidden_dim, 1)
        )
    
    def forward(self, x, edge_index, batch, sheet_node_idx, sheet_features=None):
        """
        前向传播
        
        参数:
        x: [N, node_feat_dim] 节点特征
        edge_index: [2, E] 邻接关系
        batch: [N] 批次索引(如果只有1个图，则全为0)
        sheet_node_idx: 列表，每个元素是一个LongTensor，表示一个sheet的节点索引
        sheet_features: [num_sheets, 10] 每个sheet的几何特征向量，包括:
                        [energy, boundary_ratio, feature_ratio, endpoints_special, adjacent_features,
                         curvature_metric, normal_variation, dihedral_deviation, min_jacobian, avg_jacobian]
        """
        # 1) 消息传递
        for conv in self.convs:
            x = F.relu(conv(x, edge_index))
        
        # 2) 全局图嵌入
        #    (每个batch中的图得到一个嵌入，这里batch size=1，所以形状是[1, hidden_dim])
        g_emb = global_mean_pool(x, batch)
        
        # 3) 对于每个sheet，池化其节点嵌入→sheet_emb
        sheet_embs = []
        for idx in sheet_node_idx:
            sheet_embs.append(x[idx].mean(dim=0))
        sheet_embs = torch.stack(sheet_embs, dim=0)  # [num_sheets, hidden_dim]
        
        # 3.5) 如果提供了几何特征，编码并融合
        if sheet_features is not None:
            geo_embs = self.geometry_encoder(sheet_features)  # [num_sheets, hidden_dim]
            sheet_embs = self.fusion_layer(torch.cat([sheet_embs, geo_embs], dim=1))
        
        # 4) 扩展全局嵌入以匹配sheets
        #    (仅适用于batch size=1的情况；否则，需要按图ID对sheets分组)
        global_embs = g_emb.expand(sheet_embs.size(0), -1)  # [num_sheets, hidden_dim]
        
        # 5) 拼接并预测Q值
        h = torch.cat([sheet_embs, global_embs], dim=1)  # [num_sheets, hidden_dim*2]
        q_values = self.q_mlp(h).squeeze(-1)  # [num_sheets]
        return q_values

class AttentionHexMeshQNet(nn.Module):
    def __init__(self, node_feat_dim=10, hidden_dim=64, num_heads=4, num_layers=3):
        """
        使用图注意力网络的六面体网格Q网络
        
        参数:
        node_feat_dim: 节点特征维度
        hidden_dim: 隐藏层维度
        num_heads: 注意力头数量
        num_layers: GAT层数
        """
        super().__init__()
        
        # 初始特征投影层
        self.pre_transform = nn.Linear(node_feat_dim, hidden_dim)
        
        # 多头注意力层
        self.attention_layers = nn.ModuleList()
        self.attention_layers.append(GATv2Conv(hidden_dim, hidden_dim // num_heads, heads=num_heads))
        for _ in range(num_layers-1):
            self.attention_layers.append(GATv2Conv(hidden_dim, hidden_dim // num_heads, heads=num_heads))
        
        # 几何特征编码器
        self.geometry_encoder = nn.Sequential(
            nn.Linear(node_feat_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.ReLU(inplace=True),
            nn.Linear(hidden_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.ReLU(inplace=True)
        )
        
        # 特征融合层 - 使用注意力机制
        self.cross_attention = nn.MultiheadAttention(hidden_dim, num_heads=4, batch_first=True)
        
        # Q值预测MLP
        self.q_mlp = nn.Sequential(
            nn.Linear(hidden_dim*2, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.ReLU(inplace=True),
            nn.Linear(hidden_dim, hidden_dim//2),
            nn.LayerNorm(hidden_dim//2),
            nn.ReLU(inplace=True),
            nn.Linear(hidden_dim//2, 1)
        )
    
    def forward(self, x, edge_index, batch, sheet_node_idx, sheet_features=None):
        # 初始特征投影
        x = self.pre_transform(x)
        
        # 多头注意力消息传递
        for layer in self.attention_layers:
            x = F.elu(layer(x, edge_index))
        
        # 全局图嵌入
        g_emb = global_mean_pool(x, batch)
        
        # Sheet节点池化
        sheet_embs = []
        for idx in sheet_node_idx:
            sheet_embs.append(x[idx].mean(dim=0))
        sheet_embs = torch.stack(sheet_embs, dim=0)
        
        # 如果提供了几何特征，使用交叉注意力融合
        if sheet_features is not None:
            geo_embs = self.geometry_encoder(sheet_features)
            # 使用交叉注意力进行特征融合
            sheet_embs = sheet_embs.unsqueeze(0)  # [1, num_sheets, hidden_dim]
            geo_embs = geo_embs.unsqueeze(0)  # [1, num_sheets, hidden_dim]
            
            # 应用交叉注意力
            fused_embs, _ = self.cross_attention(
                query=sheet_embs, 
                key=geo_embs, 
                value=geo_embs
            )
            sheet_embs = fused_embs.squeeze(0)  # [num_sheets, hidden_dim]
        
        # 扩展全局嵌入
        global_embs = g_emb.expand(sheet_embs.size(0), -1)
        
        # 拼接并预测Q值
        h = torch.cat([sheet_embs, global_embs], dim=1)
        q_values = self.q_mlp(h).squeeze(-1)
        
        return q_values