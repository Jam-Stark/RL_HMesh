import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, global_mean_pool


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
    def __init__(self, node_feat_dim, hidden_dim=64, num_gnn_layers=3):
        super().__init__()
        # GNN backbone
        self.convs = nn.ModuleList()
        self.convs.append(GCNConv(node_feat_dim, hidden_dim))
        for _ in range(num_gnn_layers-1):
            self.convs.append(GCNConv(hidden_dim, hidden_dim))
        # MLP to map [sheet_emb | global_emb] → Q
        self.q_mlp = nn.Sequential(
            nn.Linear(hidden_dim*2, hidden_dim),
            nn.ReLU(inplace=True),
            nn.Linear(hidden_dim, 1)
        )
    
    def forward(self, x, edge_index, batch, sheet_node_idx):
        """
        x:        [N, node_feat_dim]    node features
        edge_index: [2, E]               adjacency
        batch:    [N]                    batch indices (all zeros if 1 graph)
        sheet_node_idx: list of LongTensors,
                        each tensor gives the node‑indices of one sheet
        """
        # 1) Message‐pass
        for conv in self.convs:
            x = F.relu(conv(x, edge_index))
        # 2) Global graph embedding
        #    (one per graph in batch—here batch size is 1, so shape [1, hidden_dim])
        g_emb = global_mean_pool(x, batch)  
        # 3) For each sheet, pool its node embeddings → sheet_emb
        sheet_embs = []
        for idx in sheet_node_idx:
            sheet_embs.append(x[idx].mean(dim=0))
        sheet_embs = torch.stack(sheet_embs, dim=0)      # [num_sheets, hidden_dim]
        # 4) Expand global embedding to match sheets
        #    (only works when batch size=1; otherwise, group sheets by graph id)
        global_embs = g_emb.expand(sheet_embs.size(0), -1)  # [num_sheets, hidden_dim]
        # 5) Concatenate and predict Q
        h = torch.cat([sheet_embs, global_embs], dim=1)  # [num_sheets, hidden_dim*2]
        q_values = self.q_mlp(h).squeeze(-1)             # [num_sheets]
        return q_values