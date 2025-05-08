"""
load nn model
revcieve state from cpp
output action to cpp
update nn model

基于奖励的外部逻辑停止训练 + 部署时监控质量停止
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
    def __init__(self, eta=0.5, gamma=0.99, capacity=100000, batch_size=32, episode=0,
                 use_geometric_features=True,
                 # --- 新增/修改: 使超参数可配置 ---
                 max_steps_per_episode=500,
                 reward_threshold_to_stop=150.0,
                 no_improvement_steps_threshold=50,
                 improvement_tolerance=1.0):

        self.rewards = []
        self.eta = eta
        self.gamma = gamma
        self.capacity = capacity
        self.batch_size = batch_size
        self.episode = episode
        self.memory = ReplayMemory(self.capacity)

        self.target_update_freq = 10
        self.updates_count = 0
        self.error_filter_threshold = 0.05
        self.max_filter_threshold = 0.95
        self.filter_increase_rate = 0.01
        self.total_episodes = 0
        self.error_encounter_count = 0

        # --- 新增/修改: 存储超参数 ---
        self.max_steps_per_episode = max_steps_per_episode
        self.reward_threshold_to_stop = reward_threshold_to_stop
        self.no_improvement_steps_threshold = no_improvement_steps_threshold
        self.improvement_tolerance = improvement_tolerance

        # --- 新增: 每个 episode 需要重置的变量 ---
        self.current_episode_reward = 0.0
        self.current_episode_steps = 0
        self.best_reward_this_episode = -float('inf')
        self.steps_since_last_improvement = 0
        # --- 结束新增/修改 ---
        
        # 检测CUDA设备
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        print(f"Using device: {self.device}")
        
        # 是否使用几何特征
        self.use_geometric_features = use_geometric_features
        
        # 使用HexMeshQNet替代之前的U_linear_QNet
        if use_geometric_features:
            # 策略网络(主动更新)
            self.policy_net = network.HexMeshQNet(node_feat_dim=10, hidden_dim=64, num_gnn_layers=3).to(self.device)
            # 目标网络(延迟更新)
            self.target_net = network.HexMeshQNet(node_feat_dim=10, hidden_dim=64, num_gnn_layers=3).to(self.device)
            # 将策略网络参数复制到目标网络
            self.target_net.load_state_dict(self.policy_net.state_dict())
            self.target_net.eval()  # 设置为评估模式
            self.model = self.policy_net  # 兼容性
        else:
            self.policy_net = network.U_linear_QNet(2, 1, 10, 10, 10).to(self.device)
            self.target_net = network.U_linear_QNet(2, 1, 10, 10, 10).to(self.device)
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
            
    def remember(self, state, action, next_state, single_step_reward, done, is_error=False):
        # --- 检查是否为严重错误 (可以保留或调整) ---
        if single_step_reward < -50000: # 使用传入的单步奖励判断
            is_error = True
            self.error_encounter_count += 1
            print(f"检测到严重错误，增加记忆优先级。错误计数：{self.error_encounter_count}")
        priority = 2.0 if is_error else None # 优先级基于是否错误

        # --- 存储经验到 ReplayMemory ---
        # 注意：这里传递给 memory 的事件元组可能也需要调整，
        # 如果 memory 内部使用了 reward 来定优先级，确保它理解这是单步奖励
        # 或者直接使用这里计算的 priority
        self.memory.push((state, action, next_state, single_step_reward, done, is_error), priority) # 存储单步奖励

        # --- 新增/修改: 更新 episode 状态 ---
        self.current_episode_steps += 1  # 在这里增加步数计数
        self.current_episode_reward += single_step_reward # 累加单步奖励

        # 检查奖励是否有显著改善
        if self.current_episode_reward > self.best_reward_this_episode + self.improvement_tolerance:
            # print(f"  Improvement detected: {self.best_reward_this_episode:.2f} -> {self.current_episode_reward:.2f}") # Debug log
            self.best_reward_this_episode = self.current_episode_reward
            self.steps_since_last_improvement = 0 # 重置计数器
        else:
            self.steps_since_last_improvement += 1 # 没有显著改善，增加计数器
            # print(f"  No significant improvement for {self.steps_since_last_improvement} steps.") # Debug log
        # --- 结束新增/修改 ---

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
    
    def reset_episode_stats(self):
        """在每个episode开始时调用，重置统计信息"""
        self.current_episode_reward = 0.0
        self.current_episode_steps = 0
        self.best_reward_this_episode = -float('inf')
        self.steps_since_last_improvement = 0
        print("Episode stats reset.") # 添加日志确认

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
        # 调整探索率计算：使用指数衰减，并设置最小探索率
        min_epsilon = 0.05  # 最小探索率
        decay_rate = 0.995 # 衰减因子，越接近1衰减越慢
        start_epsilon = 0.9 # 初始探索率
        eps_threshold = max(min_epsilon, start_epsilon * (decay_rate ** episode))
        print(f"Episode: {episode}, Epsilon: {eps_threshold:.4f}") # 打印当前的探索率
        
        if sample < eps_threshold:
            # 随机选择
            action = random.randint(0, len(state) - 1)
            print(f"random action (epsilon={eps_threshold:.4f}): {action} corresponding sheet id: {state[action][0]} sheet energy: {state[action][1]}")

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
            print(f"greedy action (epsilon={eps_threshold:.4f})") # 添加日志区分贪婪选择
            with torch.no_grad():
                if self.use_geometric_features:
                    # 使用HexMeshQNet
                    features, _ = self.extract_features(state)
                    # 创建模拟的图结构数据（简化版，实际应用中需要真实的图数据）
                    num_nodes = len(state)
                    x = torch.tensor(features, dtype=torch.float32).to(self.device)  # 节点特征
                    
                    # 创建一个简单的全连接图
                    edge_index = []
                    for i in range(num_nodes):
                        for j in range(num_nodes):
                            if i != j:
                                edge_index.append([i, j])
                    edge_index = torch.tensor(edge_index, dtype=torch.long).t().to(self.device)
                    
                    # 单图批处理
                    batch = torch.zeros(num_nodes, dtype=torch.long).to(self.device)
                    
                    # 为每个sheet创建节点索引
                    sheet_node_idx = [torch.tensor([i]).to(self.device) for i in range(num_nodes)]
                    
                    # 预测Q值
                    scores = self.model(x, edge_index, batch, sheet_node_idx, x)
                else:
                    # 使用原有的U_linear_QNet
                    state_tensor = torch.tensor(state, dtype=torch.float32).to(self.device)  # shape: (n, 2)
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
                
                print(f"greedy action chosen: {action} corresponding sheet id: {state[action][0]} sheet energy: {state[action][1]}")
                return action

    def get_action(self, state): # 推理使用
        """推理时的动作选择，始终过滤非法动作"""
        with torch.no_grad():
            if self.use_geometric_features:
                # 使用HexMeshQNet
                features, _ = self.extract_features(state)
                # 创建模拟的图结构数据
                num_nodes = len(state)
                x = torch.tensor(features, dtype=torch.float32).to(self.device)
                
                # 创建一个简单的全连接图
                edge_index = []
                for i in range(num_nodes):
                    for j in range(num_nodes):
                        if i != j:
                            edge_index.append([i, j])
                edge_index = torch.tensor(edge_index, dtype=torch.long).t().to(self.device)
                
                # 单图批处理
                batch = torch.zeros(num_nodes, dtype=torch.long).to(self.device)
                
                # 为每个sheet创建节点索引
                sheet_node_idx = [torch.tensor([i]).to(self.device) for i in range(num_nodes)]
                
                # 预测Q值
                scores = self.model(x, edge_index, batch, sheet_node_idx, x)
            else:
                # 使用原有的U_linear_QNet
                state_tensor = torch.tensor(state, dtype=torch.float32).to(self.device)
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
        weights = torch.tensor(weights, dtype=torch.float32).to(self.device)
        
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
                done = batch[4][i]     # 新增：获取 done 状态 (假设 done 是元组中第5个元素，索引为4)
                is_error = batch[5][i] # 修改：假设 is_error 是元组中的第6个元素 (索引为5)
                
                # 提取当前状态和下一状态的特征
                state_features, _ = self.extract_features(state)
                next_state_features, _ = self.extract_features(next_state)
                
                # 创建图形数据
                num_nodes_state = len(state)
                num_nodes_next = len(next_state)
                
                # 当前状态的图结构
                state_x = torch.tensor(state_features, dtype=torch.float32).to(self.device)
                state_edge_index = []
                for s in range(num_nodes_state):
                    for t in range(num_nodes_state):
                        if s != t:
                            state_edge_index.append([s, t])
                state_edge_index = torch.tensor(state_edge_index, dtype=torch.long).t().to(self.device)
                state_batch = torch.zeros(num_nodes_state, dtype=torch.long).to(self.device)
                state_sheet_node_idx = [torch.tensor([j]).to(self.device) for j in range(num_nodes_state)]
                
                # 下一状态的图结构
                next_state_x = torch.tensor(next_state_features, dtype=torch.float32).to(self.device)
                next_state_edge_index = []
                for s in range(num_nodes_next):
                    for t in range(num_nodes_next):
                        if s != t:
                            next_state_edge_index.append([s, t])
                next_state_edge_index = torch.tensor(next_state_edge_index, dtype=torch.long).t().to(self.device)
                next_state_batch = torch.zeros(num_nodes_next, dtype=torch.long).to(self.device)
                next_state_sheet_node_idx = [torch.tensor([j]).to(self.device) for j in range(num_nodes_next)]
                
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
                        next_q_value = torch.tensor(0.0).to(self.device) # Default if no next state
                
                # 计算期望Q值
                if done:
                    expected_q_value = torch.tensor(reward).to(self.device)
                else:
                    expected_q_value = torch.tensor(reward).to(self.device) + self.gamma * next_q_value
                
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
            state_batch = torch.tensor([item for item in batch[0]], dtype=torch.float32).to(self.device)
            action_batch = torch.tensor([item for item in batch[1]], dtype=torch.long).to(self.device)
            next_state_batch = torch.tensor([item for item in batch[2]], dtype=torch.float32).to(self.device)
            reward_batch = torch.tensor([item for item in batch[3]], dtype=torch.float32).to(self.device)
            done_batch = torch.tensor([item for item in batch[4]], dtype=torch.bool).to(self.device) # 新增：获取 done 状态 (假设 done 是元组中第5个元素，索引为4)
            # 假设 is_error 是元组中的第6个元素 (索引5)
            error_flags = torch.tensor([item for item in batch[5]], dtype=torch.bool).to(self.device) # 修改
            
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
            # 如果 done_batch[i] is True, 那么 (self.gamma * next_q_value)[i] 应该为0
            # ~done_batch (布尔型) 在乘法中 True 行为像 1, False 像 0.
            # 所以如果 done is True, ~done_batch is False, 乘积项为0.
            # 如果 done is False, ~done_batch is True, 乘积项为 next_q_value.
            expected_q_value = reward_batch + self.gamma * next_q_value * (~done_batch)
            
            # 计算TD误差（用于更新优先级）
            td_errors = abs(expected_q_value.detach() - q_value.detach()).cpu().numpy()
            
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

    def is_finish(self, state): # state 参数用于能量检查
        """
        检查 episode 是否应该终止。
        返回 False 表示终止，True 表示继续。
        """
        # 1. 检查是否达到最大步数
        if self.current_episode_steps >= self.max_steps_per_episode:
            print(f"Episode finished: Reached max steps ({self.current_episode_steps}/{self.max_steps_per_episode}).")
            return False # 终止

        # 2. 检查累计奖励是否达到目标阈值
        if self.current_episode_reward >= self.reward_threshold_to_stop:
            print(f"Episode finished: Cumulative reward ({self.current_episode_reward:.2f}) reached threshold ({self.reward_threshold_to_stop:.2f}).")
            return False # 终止

        # 3. 检查是否长时间没有显著改善
        if self.steps_since_last_improvement >= self.no_improvement_steps_threshold:
            print(f"Episode finished: No significant reward improvement for {self.steps_since_last_improvement} steps (threshold: {self.no_improvement_steps_threshold}).")
            return False # 终止

        # 4. (可选，但推荐保留) 检查是否已无合法动作（所有 sheet 能量 <= 0）
        if not state or all(row[1] <= 0 for row in state): # state[0] 是 sheet_id, state[1] 是 energy
             print(f"Episode finished: No sheets with positive energy remaining or state is empty.")
             return False # 终止

        # 5. 如果以上条件都不满足，则继续
        return True # 继续

    def find_sheet_id(self, state, action):
        for i, row in enumerate(state):
            if row[0] == action:
                return i
        return -1

    def print_state(self, state):
        print("state:")
        for row in state:
            print(f"{row[0]} {row[1]}")
