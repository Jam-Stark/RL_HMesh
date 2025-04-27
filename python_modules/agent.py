"""
load nn model
revcieve state from cpp
output action to cpp
update nn model
"""

import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
import numpy as np
import random
import os
from utils import ReplayMemory
import network

class Agent:
    def __init__(self, eta=0.5, gamma=0.99, capacity=100000, batch_size=32, episode=0, use_geometric_features=True):
        self.rewards = []
        self.eta = eta
        self.gamma = gamma
        self.capacity = capacity
        self.batch_size = batch_size
        self.episode = episode
        self.memory = ReplayMemory(self.capacity)
        
        # 是否使用几何特征
        self.use_geometric_features = use_geometric_features
        
        # 使用HexMeshQNet替代之前的U_linear_QNet
        if use_geometric_features:
            # 策略网络(主动更新)
            self.policy_net = network.HexMeshQNet(node_feat_dim=10, hidden_dim=64, num_gnn_layers=3)
            # 目标网络(延迟更新)
            self.target_net = network.HexMeshQNet(node_feat_dim=10, hidden_dim=64, num_gnn_layers=3)
            # 将策略网络参数复制到目标网络
            self.target_net.load_state_dict(self.policy_net.state_dict())
            self.target_net.eval()  # 设置为评估模式
            self.model = self.policy_net  # 兼容性
        else:
            self.policy_net = network.U_linear_QNet(2, 1, 10, 10, 10)
            self.target_net = network.U_linear_QNet(2, 1, 10, 10, 10)
            self.target_net.load_state_dict(self.policy_net.state_dict())
            self.target_net.eval()
            self.model = self.policy_net  # 兼容性
            
        self.optimizer = optim.Adam(self.policy_net.parameters(), lr=0.001)  # 降低学习率以提高稳定性
        self.scheduler = optim.lr_scheduler.StepLR(self.optimizer, step_size=5, gamma=0.5)  # 更平缓的学习率衰减
        
        # 目标网络同步频率
        self.target_update_freq = 10  # 每10次经验回放更新一次目标网络
        self.updates_count = 0

        # 自适应错误过滤机制
        self.error_filter_threshold = 0.05  # 初始非法动作过滤概率（5%）
        self.max_filter_threshold = 0.95   # 最大过滤概率（95%）
        self.filter_increase_rate = 0.01  # 每次更新增加的过滤概率
        self.total_episodes = 0            # 记录总训练回合数
        self.error_encounter_count = 0     # 记录遇到的非法动作次数

    def save_model(self, path):
        """保存模型"""
        torch.save(self.model.state_dict(), path)
        print(f"模型已保存至 {path}")
        
    def save_memory(self, path):
        """保存经验回放记忆"""
        self.memory.save(path)
        
    def save_checkpoint(self, model_path, memory_path):
        """保存模型和经验记忆的检查点"""
        # 创建目录（如果不存在）
        os.makedirs(os.path.dirname(model_path), exist_ok=True)
        os.makedirs(os.path.dirname(memory_path), exist_ok=True)
        
        # 保存模型和记忆
        self.save_model(model_path)
        self.save_memory(memory_path)
        print(f"已完成检查点保存: 模型: {model_path}, 记忆: {memory_path}")

    def load_model(self, path):
        """加载模型"""
        try:
            self.model.load_state_dict(torch.load(path))
            self.model.eval()
            print(f"模型已从 {path} 加载")
            return True
        except FileNotFoundError:
            print(f"模型文件 {path} 不存在")
            return False
        except Exception as e:
            print(f"加载模型时发生错误: {e}")
            return False
    
    def load_memory(self, path):
        """加载经验回放记忆"""
        return self.memory.load(path)
    
    def load_checkpoint(self, model_path, memory_path):
        """加载模型和经验记忆的检查点"""
        model_loaded = self.load_model(model_path)
        memory_loaded = self.load_memory(memory_path)
        
        if model_loaded and memory_loaded:
            print("检查点加载成功")
            return True
        else:
            print("检查点加载不完整或失败")
            return False
            
    def remember(self, state, action, next_state, reward, is_error=False):
        """存储经验，带有错误标记参数"""
        # 如果奖励值非常低（比如 -100000），则认为是非法动作的严厉惩罚
        if reward < -50000:
            is_error = True
            self.error_encounter_count += 1
            print(f"检测到严重错误，增加记忆优先级。错误计数：{self.error_encounter_count}")
        
        # 为错误记忆设置较高的优先级
        priority = 2.0 if is_error else None
        self.memory.push((state, action, next_state, reward, is_error), priority)

    def extract_features(self, state):
        """从完整state中提取几何特征"""
        if not self.use_geometric_features:
            # 仅使用sheet_id和sheet_energy
            return [[item[0], item[1]] for item in state]
        
        # 从C++传来的state列表中提取所有几何特征
        # 假设state的每个元素是一个包含11个元素的列表：
        # [sheet_id, energy, boundary_ratio, feature_ratio, endpoints_special, 
        #  adjacent_features, curvature_metric, normal_variation, dihedral_deviation, 
        #  min_jacobian, avg_jacobian]
        features = []
        ids = []
        for item in state:
            if len(item) >= 11:  # 确保有所有需要的特征
                ids.append(item[0])  # sheet_id
                # 排除sheet_id，收集所有几何特征
                features.append(item[1:])
            else:
                # 如果缺少特征，填充0
                ids.append(item[0])
                features.append([item[1]] + [0.0] * 9)  # 只使用energy，其他填充0
                
        return features, ids

    def should_filter_error(self, episode):
        """根据训练进度决定是否过滤非法动作"""
        # 更新自适应过滤阈值
        self.error_filter_threshold = min(
            self.max_filter_threshold, 
            0.05 + (self.filter_increase_rate * episode)
        )
        
        # 使用随机数决定是否过滤非法动作
        return random.random() < self.error_filter_threshold

    def choose_action(self, state, episode):
        """选择动作，带有自适应错误过滤机制"""
        self.total_episodes = max(self.total_episodes, episode + 1)
        
        sample = random.random()
        eps_threshold = 0.6 / (episode + 1)
        
        if sample < eps_threshold:
            # 随机选择
            action = random.randint(0, len(state) - 1)
            print(f"random action: {action} corresponding sheet id: {state[action][0]} sheet energy: {state[action][1]}")

            # 根据自适应过滤阈值决定是否过滤非法动作
            if self.should_filter_error(episode):
                # 过滤非法动作（energy <= 0）
                attempts = 0
                while state[action][1] <= 0 and attempts < 100:  # 防止无限循环
                    action = random.randint(0, len(state) - 1)
                    print(f"过滤非法动作，重新选择: {action} sheet id: {state[action][0]} energy: {state[action][1]}")
                    attempts += 1
                
                if attempts >= 100:
                    # 如果多次尝试仍找不到合法动作，选择能量最大的sheet
                    print("多次尝试未找到合法动作，选择能量最大的sheet")
                    max_energy = float('-inf')
                    for i, item in enumerate(state):
                        if item[1] > max_energy:
                            max_energy = item[1]
                            action = i
            return action
        else:
            # 贪婪选择
            with torch.no_grad():
                if self.use_geometric_features:
                    # 使用HexMeshQNet
                    features, _ = self.extract_features(state)
                    # 创建模拟的图结构数据（简化版，实际应用中需要真实的图数据）
                    num_nodes = len(state)
                    x = torch.tensor(features, dtype=torch.float32)  # 节点特征
                    
                    # 创建一个简单的全连接图
                    edge_index = []
                    for i in range(num_nodes):
                        for j in range(num_nodes):
                            if i != j:
                                edge_index.append([i, j])
                    edge_index = torch.tensor(edge_index, dtype=torch.long).t()
                    
                    # 单图批处理
                    batch = torch.zeros(num_nodes, dtype=torch.long)
                    
                    # 为每个sheet创建节点索引
                    sheet_node_idx = [torch.tensor([i]) for i in range(num_nodes)]
                    
                    # 预测Q值
                    scores = self.model(x, edge_index, batch, sheet_node_idx, x)
                else:
                    # 使用原有的U_linear_QNet
                    state_tensor = torch.tensor(state, dtype=torch.float32)  # shape: (n, 2)
                    scores = self.model(state_tensor)  # shape: (n, 1)
                
                # 获取Q值
                q_values = scores.view(-1).cpu().numpy()
                
                # 如果启用了错误过滤且训练已达到一定程度
                if self.should_filter_error(episode):
                    # 创建掩码，排除energy <= 0的非法动作
                    valid_mask = np.array([item[1] > 0 for item in state])
                    
                    # 如果所有动作都无效，选择能量最大的
                    if not np.any(valid_mask):
                        print("所有动作均无效，选择能量最大的sheet")
                        action = max(range(len(state)), key=lambda i: state[i][1])
                    else:
                        # 将无效动作的Q值设为负无穷
                        masked_q_values = q_values.copy()
                        masked_q_values[~valid_mask] = float('-inf')
                        action = np.argmax(masked_q_values)
                else:
                    # 不过滤，直接选择Q值最高的动作
                    action = np.argmax(q_values)
                
                print(f"greedy action: {action} corresponding sheet id: {state[action][0]} sheet energy: {state[action][1]}")
                return action

    def get_action(self, state): # 推理使用
        """推理时的动作选择，始终过滤非法动作"""
        with torch.no_grad():
            if self.use_geometric_features:
                # 使用HexMeshQNet
                features, _ = self.extract_features(state)
                # 创建模拟的图结构数据
                num_nodes = len(state)
                x = torch.tensor(features, dtype=torch.float32)
                
                # 创建一个简单的全连接图
                edge_index = []
                for i in range(num_nodes):
                    for j in range(num_nodes):
                        if i != j:
                            edge_index.append([i, j])
                edge_index = torch.tensor(edge_index, dtype=torch.long).t()
                
                # 单图批处理
                batch = torch.zeros(num_nodes, dtype=torch.long)
                
                # 为每个sheet创建节点索引
                sheet_node_idx = [torch.tensor([i]) for i in range(num_nodes)]
                
                # 预测Q值
                scores = self.model(x, edge_index, batch, sheet_node_idx, x)
            else:
                # 使用原有的U_linear_QNet
                state_tensor = torch.tensor(state, dtype=torch.float32)
                scores = self.model(state_tensor)
            
            # 获取Q值并创建掩码，排除energy <= 0的非法动作
            q_values = scores.view(-1).cpu().numpy()
            valid_mask = np.array([item[1] > 0 for item in state])
            
            # 如果所有动作都无效，选择能量最大的
            if not np.any(valid_mask):
                print("所有动作均无效，选择能量最大的sheet")
                action = max(range(len(state)), key=lambda i: state[i][1])
            else:
                # 将无效动作的Q值设为负无穷
                masked_q_values = q_values.copy()
                masked_q_values[~valid_mask] = float('-inf')
                action = np.argmax(masked_q_values)
                
            return action

    def replay(self):
        """经验回放，进行模型训练"""
        if len(self.memory) < self.batch_size:
            return
            
        # 使用优先级采样
        transitions, weights, indices = self.memory.sample(self.batch_size)
        batch = list(zip(*transitions))
        
        # 如果没有返回权重，则创建默认权重
        if weights is None:
            weights = np.ones(self.batch_size)
        weights = torch.tensor(weights, dtype=torch.float32)
        
        if self.use_geometric_features:
            # 处理几何特征的版本
            # 由于图数据的复杂性，我们采用更简单的方法处理每个批次
            losses = []
            errors = []  # 存储TD误差，用于更新优先级
            
            for i in range(self.batch_size):
                state = batch[0][i]
                action = batch[1][i]
                next_state = batch[2][i]
                reward = batch[3][i]
                # 检查是否为错误经验（如果有第5个元素）
                is_error = len(batch) > 4 and batch[4][i]
                
                # 提取当前状态和下一状态的特征
                state_features, _ = self.extract_features(state)
                next_state_features, _ = self.extract_features(next_state)
                
                # 创建图形数据
                num_nodes_state = len(state)
                num_nodes_next = len(next_state)
                
                # 当前状态的图结构
                state_x = torch.tensor(state_features, dtype=torch.float32)
                state_edge_index = []
                for s in range(num_nodes_state):
                    for t in range(num_nodes_state):
                        if s != t:
                            state_edge_index.append([s, t])
                state_edge_index = torch.tensor(state_edge_index, dtype=torch.long).t()
                state_batch = torch.zeros(num_nodes_state, dtype=torch.long)
                state_sheet_node_idx = [torch.tensor([j]) for j in range(num_nodes_state)]
                
                # 下一状态的图结构
                next_state_x = torch.tensor(next_state_features, dtype=torch.float32)
                next_state_edge_index = []
                for s in range(num_nodes_next):
                    for t in range(num_nodes_next):
                        if s != t:
                            next_state_edge_index.append([s, t])
                next_state_edge_index = torch.tensor(next_state_edge_index, dtype=torch.long).t()
                next_state_batch = torch.zeros(num_nodes_next, dtype=torch.long)
                next_state_sheet_node_idx = [torch.tensor([j]) for j in range(num_nodes_next)]
                
                # 使用策略网络计算当前状态的Q值
                state_q_values = self.policy_net(state_x, state_edge_index, state_batch, 
                                              state_sheet_node_idx, state_x)
                current_q_value = state_q_values[action]
                
                # 使用目标网络计算下一状态的最大Q值
                with torch.no_grad():
                    if num_nodes_next > 0:  # 确保下一状态不为空
                        next_state_q_values = self.target_net(next_state_x, next_state_edge_index, 
                                                            next_state_batch, next_state_sheet_node_idx, 
                                                            next_state_x)
                        # 使用Double DQN技术 - 避免Q值过估计
                        # 首先使用策略网络找出最佳动作
                        next_policy_values = self.policy_net(next_state_x, next_state_edge_index,
                                                          next_state_batch, next_state_sheet_node_idx,
                                                          next_state_x)
                        best_action = next_policy_values.argmax().item()
                        # 然后使用目标网络计算该动作的Q值
                        next_q_value = next_state_q_values[best_action]
                    else:
                        next_q_value = torch.tensor(0.0)
                
                # 计算期望Q值，使用截断的n步Q值以提高稳定性
                expected_q_value = torch.tensor(reward) + self.gamma * next_q_value
                
                # 计算TD误差
                td_error = abs((expected_q_value - current_q_value).item())
                errors.append(td_error)
                
                # 使用Huber损失以更好地处理异常值，应用重要性采样权重
                loss = F.smooth_l1_loss(current_q_value.unsqueeze(0), expected_q_value.unsqueeze(0)) * weights[i]
                losses.append(loss)
                
                # 对错误状态给予更多关注
                if is_error:
                    # 增强错误经验的训练强度（增加一次学习）
                    losses.append(loss * 2.0)  # 额外添加一次同样的损失
                
            # 组合所有批次的损失并反向传播
            if losses:
                total_loss = sum(losses) / len(losses)
                self.optimizer.zero_grad()
                total_loss.backward()
                
                # 梯度裁剪以防止梯度爆炸
                torch.nn.utils.clip_grad_norm_(self.policy_net.parameters(), max_norm=1.0)
                
                self.optimizer.step()
                
                # 更新优先级
                self.memory.update_priorities(indices, errors)
                
                # 记录更新次数
                self.updates_count += 1
                
                # 定期同步目标网络
                if self.updates_count % self.target_update_freq == 0:
                    self.target_net.load_state_dict(self.policy_net.state_dict())
                    print(f"目标网络已更新，更新次数: {self.updates_count}")
                
                # 每个epoch结束时更新学习率（不是每个batch）
                if self.updates_count % 100 == 0:
                    self.scheduler.step()
                    print(f"学习率已更新为: {self.optimizer.param_groups[0]['lr']}")
                    print(f"当前错误过滤阈值: {self.error_filter_threshold:.2f}, 总训练回合数: {self.total_episodes}")
        else:
            # 原始版本，无几何特征，也更新为使用目标网络和优先级采样的版本
            state_batch = torch.tensor([item for item in batch[0]], dtype=torch.float32)
            action_batch = torch.tensor([item for item in batch[1]], dtype=torch.long)
            next_state_batch = torch.tensor([item for item in batch[2]], dtype=torch.float32)
            reward_batch = torch.tensor([item for item in batch[3]], dtype=torch.float32)
            
            # 检查是否有错误标志
            error_flags = torch.tensor([len(batch) > 4 and item for item in batch[4]], dtype=torch.bool) if len(batch) > 4 else None
            
            self.optimizer.zero_grad()
            
            # 使用策略网络计算当前Q值
            q_values = self.policy_net(state_batch).squeeze(-1)
            q_value = q_values.gather(1, action_batch.unsqueeze(1)).squeeze(1)
            
            # 使用目标网络计算下一状态最大Q值（Double DQN）
            with torch.no_grad():
                # 使用策略网络选择动作
                next_q_actions = self.policy_net(next_state_batch).squeeze(-1).max(1)[1].unsqueeze(1)
                # 使用目标网络评估动作价值
                next_q_values = self.target_net(next_state_batch).squeeze(-1)
                next_q_value = next_q_values.gather(1, next_q_actions).squeeze(1)
                
            # 计算期望Q值
            expected_q_value = reward_batch + self.gamma * next_q_value
            
            # 计算TD误差（用于更新优先级）
            td_errors = abs(expected_q_value.detach() - q_value.detach()).numpy()
            
            # 计算损失时应用重要性采样权重
            loss = (F.smooth_l1_loss(q_value, expected_q_value, reduction='none') * weights).mean()
            
            # 如果有错误标志，为错误样本增加额外的学习强度
            if error_flags is not None and torch.any(error_flags):
                # 获取错误样本的索引
                error_indices = torch.nonzero(error_flags).squeeze(-1)
                if len(error_indices) > 0:
                    # 为错误样本计算额外的损失
                    error_q_value = q_value[error_indices]
                    error_expected_q = expected_q_value[error_indices]
                    error_weights = weights[error_indices]
                    error_loss = (F.smooth_l1_loss(error_q_value, error_expected_q, reduction='none') * error_weights * 2.0).mean()
                    
                    # 添加到总损失中
                    loss = loss + error_loss
            
            loss.backward()
            torch.nn.utils.clip_grad_norm_(self.policy_net.parameters(), max_norm=1.0)
            self.optimizer.step()
            
            # 更新优先级
            self.memory.update_priorities(indices, td_errors)
            
            # 记录更新并同步目标网络
            self.updates_count += 1
            if self.updates_count % self.target_update_freq == 0:
                self.target_net.load_state_dict(self.policy_net.state_dict())
                print(f"目标网络已更新，更新次数: {self.updates_count}")
                
            # 定期更新学习率
            if self.updates_count % 100 == 0:
                self.scheduler.step()
                print(f"学习率已更新为: {self.optimizer.param_groups[0]['lr']}")
                print(f"当前错误过滤阈值: {self.error_filter_threshold:.2f}, 总训练回合数: {self.total_episodes}")

    def is_finish(self, state):
        #print("is_finish running")
        return not all(row[1] <= 0 for row in state) #if state can not continue, return True, so add "not" here

    def find_sheet_id(self, state, action):
        for i, row in enumerate(state):
            if row[0] == action:
                return i
        return -1

    def print_state(self, state):
        print("state:")
        for row in state:
            print(f"{row[0]} {row[1]}")
