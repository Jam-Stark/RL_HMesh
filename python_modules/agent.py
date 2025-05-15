"""
load nn model
revcieve state from cpp
output action to cpp
update nn model

模型主动输出终止信号，辅助cpp判断严重error终止条件
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
                 max_steps_per_episode=500,
                 reward_threshold_to_stop=150.0,
                 no_improvement_steps_threshold=50,
                 improvement_tolerance=1.0,
                 done_loss_weight=0.7):

        # ... (previous __init__ content like self.rewards, self.device, etc.) ...
        # Remove cache variables if they were added for the previous optimization:
        # self._current_state_raw_for_cache = None
        # self._current_step_model_output_cache = None

        self.rewards = []
        self.eta = eta
        self.gamma = gamma
        self.capacity = capacity
        self.batch_size = batch_size
        self.episode = episode # Tracks total episodes processed by agent
        self.memory = ReplayMemory(self.capacity)

        self.target_update_freq = 10
        self.updates_count = 0
        self.error_filter_threshold = 0.05
        self.max_filter_threshold = 0.95
        self.filter_increase_rate = 0.018 
        self.total_episodes = 0 # Can be used to track total episodes agent has seen for epsilon decay
        self.error_encounter_count = 0
        
        # These are now for C++ main loop, but agent can still track its own version if needed for 'remember'
        self.max_steps_per_episode_param = max_steps_per_episode 
        self.reward_threshold_to_stop_param = reward_threshold_to_stop
        self.no_improvement_steps_threshold_param = no_improvement_steps_threshold
        self.improvement_tolerance = improvement_tolerance

        self.current_episode_reward = 0.0
        self.current_episode_steps = 0 # Agent's own step counter for an episode
        self.best_reward_this_episode = -float('inf')
        self.steps_since_last_improvement = 0
        
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        print(f"Using device: {self.device}")
        
        self.use_geometric_features = use_geometric_features
        
        if use_geometric_features:
            self.policy_net = network.AttentionHexMeshQNet(node_feat_dim=10, hidden_dim=64, num_heads=4, num_layers=3).to(self.device)
            self.target_net = network.AttentionHexMeshQNet(node_feat_dim=10, hidden_dim=64, num_heads=4, num_layers=3).to(self.device)
        else:
            # If U_linear_QNet is also modified to output a done signal, reflect that here.
            # For now, assuming it needs to be adapted or this path is handled differently.
            self.policy_net = network.U_linear_QNet(2, 1, 10, 10, 10).to(self.device) 
            self.target_net = network.U_linear_QNet(2, 1, 10, 10, 10).to(self.device)
        
        self.target_net.load_state_dict(self.policy_net.state_dict())
        self.target_net.eval()
        self.model = self.policy_net 

        self.optimizer = optim.Adam(self.policy_net.parameters(), lr=0.001)
        self.scheduler = optim.lr_scheduler.StepLR(self.optimizer, step_size=3, gamma=0.5)
        
        self.done_loss_fn = nn.BCEWithLogitsLoss()
        self.done_loss_weight = done_loss_weight

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
            
# ...existing code...
    def remember(self, state, action, next_state, single_step_reward, done, is_error=False):
        #print("DEBUG: remember: Entered function") # 新增
        # --- 检查是否为严重错误 (可以保留或调整) ---
        if single_step_reward < -50000: # 使用传入的单步奖励判断
            is_error = True
            self.error_encounter_count += 1
            print(f"检测到严重错误，增加记忆优先级。错误计数：{self.error_encounter_count}")
        priority = 2.0 if is_error else None # 优先级基于是否错误
        #print(f"DEBUG: remember: Calculated priority: {priority}") # 新增

        # --- 存储经验到 ReplayMemory ---
        #print("DEBUG: remember: About to call self.memory.push()") # 新增
        self.memory.push((state, action, next_state, single_step_reward, done, is_error), priority) # 存储单步奖励
        #print("DEBUG: remember: Returned from self.memory.push()") # 新增
        #print(f"存储经验: {state}, {action}, {next_state}, {single_step_reward}, {done}, {is_error}") # Debug log
        
        # --- 新增/修改: 更新 episode 状态 ---
        #print("DEBUG: remember: About to update episode stats") # 新增
        self.current_episode_steps += 1  # 在这里增加步数计数
        self.current_episode_reward += single_step_reward # 累加单步奖励
       # print(f"DEBUG: remember: Updated episode_steps: {self.current_episode_steps}, episode_reward: {self.current_episode_reward}") # 新增

        # 检查奖励是否有显著改善
        if self.current_episode_reward > self.best_reward_this_episode + self.improvement_tolerance:
            # print(f"  Improvement detected: {self.best_reward_this_episode:.2f} -> {self.current_episode_reward:.2f}") # Debug log
            self.best_reward_this_episode = self.current_episode_reward
            self.steps_since_last_improvement = 0 # 重置计数器
            #print("DEBUG: remember: Improvement detected.") # 新增
        else:
            self.steps_since_last_improvement += 1 # 没有显著改善，增加计数器
            #print(f"DEBUG: remember: No significant improvement. Steps since last: {self.steps_since_last_improvement}") # 新增
        # --- 结束新增/修改 ---
       # print("DEBUG: remember: Exiting function") # 新增
# ...existing code...

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

    def choose_action(self, state_py_list, episode_num_for_epsilon): # Renamed 'episode' to avoid conflict
        """
        选择动作，并返回网络对当前状态是否应终止的预测。
        返回: (action_idx, state_done_logit)
        """
        # self.total_episodes can be used for epsilon decay based on agent's lifetime episodes
        # or use episode_num_for_epsilon if it's specific to current training run of a file.
        # Using agent's self.episode for overall decay.
        
        q_values_scores = None
        state_done_logit = torch.tensor([0.0]).to(self.device) # Default done logit (e.g. not done)

        # Perform network pass to get Q-values and state_done_logit
        # This pass happens regardless of epsilon-greedy for consistent done_logit output
        if self.use_geometric_features:
            if not state_py_list: # Handle empty state early
                print("Warning: choose_action called with empty state_py_list (geometric).")
                return 0, state_done_logit.item() # Default action, default done logit

            features, _ = self.extract_features(state_py_list)
            num_nodes = len(features)
            if num_nodes == 0:
                # print("Greedy choice: state became empty, no action.")
                return 0, state_done_logit.item() # Default action, default done logit

            x = torch.tensor(features, dtype=torch.float32).to(self.device)
            edge_list = [[i, j] for i in range(num_nodes) for j in range(num_nodes) if i != j]
            edge_index = torch.tensor(edge_list, dtype=torch.long).t().contiguous().to(self.device) if edge_list else torch.empty((2,0), dtype=torch.long).to(self.device)
            batch_tensor = torch.zeros(num_nodes, dtype=torch.long).to(self.device)
            sheet_node_idx_tensor = [torch.tensor([i]).to(self.device) for i in range(num_nodes)]
            sheet_features_tensor = x
            
            with torch.no_grad(): # Ensure no_grad for inference part of action selection
                q_values_scores, state_done_logit = self.model(x, edge_index, batch_tensor, sheet_node_idx_tensor, sheet_features_tensor)
        
        else: # U_linear_QNet path
            if not state_py_list:
                print("Warning: choose_action called with empty state_py_list (non-geometric).")
                return 0, state_done_logit.item()
            
            state_tensor = torch.tensor(state_py_list, dtype=torch.float32).to(self.device)
            if state_tensor.shape[0] == 0:
                 return 0, state_done_logit.item()

            with torch.no_grad():
                network_output = self.model(state_tensor)
                if isinstance(network_output, tuple) and len(network_output) == 2: # Assuming U_linear_QNet also modified
                    q_values_scores, state_done_logit = network_output
                else: # Original U_linear_QNet or unexpected output
                    q_values_scores = network_output
                    # state_done_logit remains default (e.g., 0.0 or not done)
                    # Or handle error if done_logit is expected but not received

        # Epsilon-greedy action selection logic
        sample = random.random()
        # Use episode_num_for_epsilon passed from C++ for current file's training progress
        min_epsilon = 0.05  
        decay_rate = 0.98 
        start_epsilon = 1 
        eps_threshold = max(min_epsilon, start_epsilon * (decay_rate ** episode_num_for_epsilon))

        action_idx = 0
        if not state_py_list: # Should have been caught earlier, but as a safeguard
            return 0, state_done_logit.item() if torch.is_tensor(state_done_logit) else float(state_done_logit)


        if sample < eps_threshold: # Random action
            print(f"random action (epsilon={eps_threshold:.4f})")
            if len(state_py_list) > 0:
                action_idx = random.randint(0, len(state_py_list) - 1)
            
            # Filtering for random action (optional, but was in your original code)
            if len(state_py_list) > 0 and self.should_filter_error(episode_num_for_epsilon): # use episode_num_for_epsilon
                attempts = 0
                # Ensure state_py_list[action_idx] is valid
                while attempts < 100 and (action_idx >= len(state_py_list) or state_py_list[action_idx][1] <= 0):
                    if len(state_py_list) > 0:
                        action_idx = random.randint(0, len(state_py_list) - 1)
                    attempts += 1
                if attempts >= 20 and len(state_py_list) > 0: # Max attempts, pick best available
                    energies = [item[1] for item in state_py_list]
                    if any(e > 0 for e in energies): # if any positive energy
                         action_idx = np.argmax(energies) # Choose max energy
                    # else: action_idx remains as is, or pick based on highest (least negative) energy
        
        else: # Greedy action
            print(f"greedy action (epsilon={eps_threshold:.4f})")
            if q_values_scores is not None and q_values_scores.numel() > 0:
                q_values_numpy = q_values_scores.view(-1).cpu().numpy()
                
                if self.should_filter_error(episode_num_for_epsilon): # use episode_num_for_epsilon
                    valid_mask = np.array([item[1] > 0 for item in state_py_list])
                    if not np.any(valid_mask) and len(state_py_list) > 0:
                        energies = [item[1] for item in state_py_list]
                        action_idx = np.argmax(energies) 
                    elif np.any(valid_mask):
                        masked_q_values = q_values_numpy.copy()
                        masked_q_values[~valid_mask] = float('-inf')
                        action_idx = np.argmax(masked_q_values)
                    # else: if state_py_list is empty, action_idx remains 0
                elif len(state_py_list) > 0:
                    action_idx = np.argmax(q_values_numpy)
                # else: if state_py_list is empty, action_idx remains 0
            elif len(state_py_list) > 0: # Q-scores are None or empty, but state is not
                action_idx = random.randint(0, len(state_py_list) - 1) # Fallback to random

        # Ensure action_idx is valid for state_py_list before printing
        if len(state_py_list) > 0 and action_idx < len(state_py_list):
            print(f"Chosen action: {action_idx} (Sheet ID: {state_py_list[action_idx][0]}), Done Logit: {state_done_logit.item():.4f}")
        elif len(state_py_list) == 0:
            print(f"Chosen action: {action_idx} (State empty), Done Logit: {state_done_logit.item():.4f}")
        else: # action_idx out of bounds or other issue
            print(f"Chosen action: {action_idx} (Error with state or action index), Done Logit: {state_done_logit.item():.4f}")

        # Return the chosen action index and the state_done_logit (as a Python float)
        final_logit = state_done_logit.item() if torch.is_tensor(state_done_logit) else float(state_done_logit)
        return action_idx, final_logit


    def get_action(self, state): # 推理使用
        """推理时的动作选择，始终过滤非法动作"""
        with torch.no_grad():
            if self.use_geometric_features:
                features, _ = self.extract_features(state)
                num_nodes = len(state) # Changed from len(features) to len(state) for consistency with choose_action
                if num_nodes == 0: # Handle empty state
                    print("State is empty in get_action. No action possible.")
                    # 根据您的逻辑返回一个默认动作或者抛出错误
                    # 例如，如果环境中没有可选动作，可能返回一个特定的“无操作”代码
                    return 0 # 或者其他合适的值，例如 None 或 -1，取决于您的C++代码如何处理

                x = torch.tensor(features, dtype=torch.float32).to(self.device)
                
                edge_index_list = [] # Renamed for clarity
                for i in range(num_nodes):
                    for j in range(num_nodes):
                        if i != j:
                            edge_index_list.append([i, j])
                
                # Ensure edge_index is correctly shaped even if edge_index_list is empty
                if not edge_index_list:
                    edge_index = torch.empty((2, 0), dtype=torch.long).to(self.device)
                else:
                    edge_index = torch.tensor(edge_index_list, dtype=torch.long).t().contiguous().to(self.device)
                
                batch_tensor = torch.zeros(num_nodes, dtype=torch.long).to(self.device)
                sheet_node_idx_tensor = [torch.tensor([i]).to(self.device) for i in range(num_nodes)]
                # Assuming sheet_features are the same as node features 'x' for AttentionHexMeshQNet
                sheet_features_tensor = x 

                # --- MODIFIED SECTION START ---
                # self.model (policy_net) 现在返回 q_values 和 state_done_logit
                # 我们主要需要 q_values_scores 来选择动作。
                # state_done_logit 在这里被接收但没有直接用于选择动作的索引。
                # 如果推理循环需要知道网络是否认为当前状态是终止状态，可以考虑让 get_action 也返回这个信息。
                q_values_scores, _state_done_logit = self.model(x, edge_index, batch_tensor, sheet_node_idx_tensor, sheet_features_tensor)
                # _state_done_logit 被接收以匹配网络输出，但如果 get_action 仅用于选择动作，则可能不直接使用它。
                # --- MODIFIED SECTION END ---
            else: # U_linear_QNet path
                state_tensor = torch.tensor(state, dtype=torch.float32).to(self.device)
                if state_tensor.shape[0] == 0 : # check for empty state in non-geo path
                    print("State is empty in get_action (non-geometric). No action possible.")
                    return 0 # 或者其他合适的值
                
                # --- MODIFIED SECTION START ---
                # 假设 U_linear_QNet 也被修改以返回 (q_values, state_done_logit)
                # 如果 U_linear_QNet 没有修改，这里的逻辑需要调整回只接收一个输出
                network_output = self.model(state_tensor)
                if isinstance(network_output, tuple) and len(network_output) == 2:
                    q_values_scores, _state_done_logit = network_output
                else: # 假设原始行为或意外的输出
                    q_values_scores = network_output
                    # _state_done_logit 在这种情况下未定义，如果后续逻辑依赖它，需要处理
                # --- MODIFIED SECTION END ---

            q_values = q_values_scores.view(-1).cpu().numpy()
            
            # 如果q_values为空 (例如，如果num_nodes为0导致网络没有有效输入/输出)
            if q_values.size == 0:
                if len(state) > 0: # 状态不为空，但Q值为空，可能是网络问题或num_nodes=0后处理问题
                    print("Warning: Q-values are empty in get_action, but state is not. Choosing random valid action.")
                    #尝试选择一个随机的合法动作作为后备
                    valid_actions = [idx for idx, item in enumerate(state) if item[1] > 0]
                    if valid_actions:
                        return random.choice(valid_actions)
                    elif len(state) > 0: # 没有合法的正能量动作，但有动作可选
                        return random.randint(0, len(state) -1)
                    else: # 状态也为空
                        return 0 # 或者其他合适的值
                else: # 状态和Q值都为空
                     print("State and Q-values are empty in get_action.")
                     return 0 # 或者其他合适的值


            valid_mask = np.array([item[1] > 0 for item in state])
            
            action = 0 # Default action
            if len(state) > 0 : # 确保状态列表不为空
                if not np.any(valid_mask):
                    print("所有动作均无效(无正能量)，选择(能量)最大的sheet (可能是负能量或0)")
                    # 即使所有都是非正能量，也选择一个，通常是能量“最高”（最不差）的那个
                    action = np.argmax([item[1] for item in state]) if state else 0
                else:
                    masked_q_values = q_values.copy()
                    masked_q_values[~valid_mask] = float('-inf') # 将无效动作的Q值设为负无穷
                    action = np.argmax(masked_q_values)
            else:
                print("State is empty at the end of get_action, returning default action 0.")

            return action

    def replay(self):
        """经验回放，进行模型训练"""
        print("DEBUG: replay: Entered function") 
        if len(self.memory) < self.batch_size:
            return
            
        transitions, weights, indices = self.memory.sample(self.batch_size)
        batch = list(zip(*transitions))
        
        if weights is None:
            weights = np.ones(self.batch_size)
        weights_tensor = torch.tensor(weights, dtype=torch.float32).to(self.device) # Renamed to avoid conflict
        
        if self.use_geometric_features: # Assuming AttentionHexMeshQNet is used here
            q_losses = []
            done_losses = [] # Store done losses separately
            td_errors_list = [] 
            
            for i in range(self.batch_size):
                state, action, next_state, reward, done_flag, is_error = batch[0][i], batch[1][i], batch[2][i], batch[3][i], batch[4][i], batch[5][i]

                # --- Prepare current state tensors ---
                state_features, _ = self.extract_features(state)
                num_nodes_state = len(state_features)
                if num_nodes_state == 0: continue # Skip if state is empty
                state_x = torch.tensor(state_features, dtype=torch.float32).to(self.device)
                state_edge_list = [[s, t] for s in range(num_nodes_state) for t in range(num_nodes_state) if s != t]
                state_edge_index = torch.tensor(state_edge_list, dtype=torch.long).t().contiguous().to(self.device) if state_edge_list else torch.empty((2,0), dtype=torch.long).to(self.device)
                state_batch_tensor = torch.zeros(num_nodes_state, dtype=torch.long).to(self.device)
                state_sheet_node_idx = [torch.tensor([j]).to(self.device) for j in range(num_nodes_state)]
                
                # --- Prepare next state tensors ---
                next_state_features, _ = self.extract_features(next_state)
                num_nodes_next = len(next_state_features)
                next_state_x = torch.tensor(next_state_features, dtype=torch.float32).to(self.device) if num_nodes_next > 0 else torch.empty((0, state_x.shape[1])).to(self.device)
                next_state_edge_list = [[s, t] for s in range(num_nodes_next) for t in range(num_nodes_next) if s != t]
                next_state_edge_index = torch.tensor(next_state_edge_list, dtype=torch.long).t().contiguous().to(self.device) if num_nodes_next > 0 and next_state_edge_list else torch.empty((2,0), dtype=torch.long).to(self.device)
                next_state_batch_tensor = torch.zeros(num_nodes_next, dtype=torch.long).to(self.device) if num_nodes_next > 0 else torch.empty((0), dtype=torch.long).to(self.device)
                next_state_sheet_node_idx = [torch.tensor([j]).to(self.device) for j in range(num_nodes_next)] if num_nodes_next > 0 else []

                # --- Q-value for current state and predicted done for current state ---
                # policy_net output: q_action_values, state_level_done_logit
                current_q_action_values, current_state_done_logit = self.policy_net(state_x, state_edge_index, state_batch_tensor, state_sheet_node_idx, state_x)
                current_q_value = current_q_action_values[action]

                # --- Done loss for current state ---
                done_target = torch.tensor([1.0 if done_flag else 0.0], dtype=torch.float32).to(self.device)
                # Ensure current_state_done_logit is correctly shaped for BCEWithLogitsLoss ([N] or [N,1])
                # state_done_predictor outputs [1] for single instance, so it should be fine.
                # If g_emb was [batch, hidden_dim] -> done_logit [batch], so for single item, it's [1]
                loss_d = self.done_loss_fn(current_state_done_logit, done_target) 
                done_losses.append(loss_d * weights_tensor[i])


                # --- Next Q-value using Double DQN ---
                next_q_value_max = torch.tensor(0.0).to(self.device) # Default if next state is terminal or empty
                if not done_flag and num_nodes_next > 0 :
                    with torch.no_grad():
                        # Target net for Q-value evaluation of next state
                        next_q_target_action_values, _ = self.target_net(next_state_x, next_state_edge_index, next_state_batch_tensor, next_state_sheet_node_idx, next_state_x)
                        # Policy net for best action selection in next state
                        next_q_policy_action_values, _ = self.policy_net(next_state_x, next_state_edge_index, next_state_batch_tensor, next_state_sheet_node_idx, next_state_x)
                        
                        if next_q_policy_action_values.numel() > 0: # Check if there are any q-values (i.e. actions)
                            best_next_action = next_q_policy_action_values.argmax().item()
                            next_q_value_max = next_q_target_action_values[best_next_action]
                            # If next_q_policy_action_values is empty, next_q_value_max remains 0.0
                
                expected_q_value = torch.tensor(reward, dtype=torch.float32).to(self.device)
                if not done_flag:
                    expected_q_value += self.gamma * next_q_value_max
                
                loss_q = F.smooth_l1_loss(current_q_value.unsqueeze(0), expected_q_value.unsqueeze(0)) * weights_tensor[i]
                q_losses.append(loss_q)
                
                td_errors_list.append(abs((expected_q_value - current_q_value).item()))
                
                if is_error: # Augment loss for error states
                    q_losses.append(loss_q * 2.0) # Add Q-loss again
                    if done_losses: done_losses.append(loss_d * weights_tensor[i] * 2.0) # Add done_loss again

            if q_losses: # Check if there are any losses to process
                total_q_loss = sum(q_losses) / len(q_losses) if q_losses else torch.tensor(0.0).to(self.device)
                total_done_loss = sum(done_losses) / len(done_losses) if done_losses else torch.tensor(0.0).to(self.device)
                
                total_loss = total_q_loss + self.done_loss_weight * total_done_loss

                self.optimizer.zero_grad()
                total_loss.backward()
                torch.nn.utils.clip_grad_norm_(self.policy_net.parameters(), max_norm=1.0)
                self.optimizer.step()
                
                self.memory.update_priorities(indices, td_errors_list) # Use td_errors_list
        
        else: # Non-geometric features path (U_linear_QNet)
            # This part assumes U_linear_QNet is NOT modified to output a done signal.
            # If it IS modified, this section needs similar changes to the geometric path.
            # For brevity, I'll assume it's not changed for now. If it is, let me know.
            state_b = torch.tensor([item for item in batch[0]], dtype=torch.float32).to(self.device)
            action_b = torch.tensor([item for item in batch[1]], dtype=torch.long).to(self.device)
            next_state_b = torch.tensor([item for item in batch[2]], dtype=torch.float32).to(self.device)
            reward_b = torch.tensor([item for item in batch[3]], dtype=torch.float32).to(self.device)
            done_b = torch.tensor([item for item in batch[4]], dtype=torch.bool).to(self.device)
            error_flags = torch.tensor([item for item in batch[5]], dtype=torch.bool).to(self.device)

            self.optimizer.zero_grad()
            
            # --- Current Q-values ---
            # Assuming U_linear_QNet only outputs Q-values
            q_values_policy_all = self.policy_net(state_b).squeeze(-1) # Shape: [B, num_actions] or [B] if num_actions=1
            # If U_linear_QNet returns (q_values, done_logit_state), unpack here:
            # q_values_policy_all, done_logits_policy_state = self.policy_net(state_b)
            # q_values_policy_all = q_values_policy_all.squeeze(-1)
            # And calculate done_loss:
            # done_targets_b = done_b.float() # if done_logits_policy_state is [B]
            # done_loss_b = self.done_loss_fn(done_logits_policy_state, done_targets_b)

            # Adjust gather for q_values_policy_all shape
            if q_values_policy_all.ndim == 1: # if model output is [B] (e.g. single action value per state)
                q_value = q_values_policy_all 
            else: # if model output is [B, num_actions]
                q_value = q_values_policy_all.gather(1, action_b.unsqueeze(1)).squeeze(1)


            # --- Next Q-values (Double DQN) ---
            with torch.no_grad():
                # Assuming U_linear_QNet only outputs Q-values
                next_q_actions_policy_all = self.policy_net(next_state_b).squeeze(-1)
                # If U_linear_QNet returns (q_values, _), unpack here:
                # next_q_actions_policy_all, _ = self.policy_net(next_state_b) 
                # next_q_actions_policy_all = next_q_actions_policy_all.squeeze(-1)

                # Assuming U_linear_QNet only outputs Q-values
                next_q_values_target_all = self.target_net(next_state_b).squeeze(-1)
                # If U_linear_QNet returns (q_values, _), unpack here:
                # next_q_values_target_all, _ = self.target_net(next_state_b)
                # next_q_values_target_all = next_q_values_target_all.squeeze(-1)


                if next_q_actions_policy_all.ndim == 1: # if model output is [B]
                    # This case implies a fixed action or single Q value per state from U_linear_QNet, which is unusual for DQN with multiple actions.
                    # Max over non-existent dimension or direct use depends on U_linear_QNet's exact output for a state.
                    # For simplicity, if it's [B], we assume it's the value of the best action already or a single action.
                    next_q_value_target = next_q_values_target_all # Or handle appropriately if actions need to be chosen
                else: # if model output is [B, num_actions]
                    best_next_actions = next_q_actions_policy_all.max(1)[1].unsqueeze(1)
                    next_q_value_target = next_q_values_target_all.gather(1, best_next_actions).squeeze(1)


            expected_q_value = reward_b + self.gamma * next_q_value_target * (~done_b)
            
            td_errors_numpy = abs(expected_q_value.detach() - q_value.detach()).cpu().numpy()
            
            loss_q_b = (F.smooth_l1_loss(q_value, expected_q_value, reduction='none') * weights_tensor).mean()
            
            # --- Combine with done loss if U_linear_QNet was modified ---
            # total_loss_b = loss_q_b + self.done_loss_weight * done_loss_b 
            # For now, assuming no done_loss_b for U_linear_QNet
            total_loss_b = loss_q_b

            if error_flags is not None and torch.any(error_flags):
                error_indices = torch.nonzero(error_flags).squeeze(-1)
                if len(error_indices) > 0:
                    error_q_value = q_value[error_indices]
                    error_expected_q = expected_q_value[error_indices]
                    error_weights = weights_tensor[error_indices] # use weights_tensor
                    # Recalculate error loss with reduction='none' before multiplying by weights
                    error_loss_elements = F.smooth_l1_loss(error_q_value, error_expected_q, reduction='none')
                    error_loss = (error_loss_elements * error_weights * 2.0).mean()
                    total_loss_b = total_loss_b + error_loss
            
            total_loss_b.backward()
            torch.nn.utils.clip_grad_norm_(self.policy_net.parameters(), max_norm=1.0)
            self.optimizer.step()
            
            self.memory.update_priorities(indices, td_errors_numpy)

        # --- Common update logic ---
        self.updates_count += 1
        if self.updates_count % self.target_update_freq == 0:
            self.target_net.load_state_dict(self.policy_net.state_dict())
            print(f"目标网络已更新，更新次数: {self.updates_count}")
            
        if self.updates_count % 10 == 0:
            self.scheduler.step()
            print(f"学习率已更新为: {self.optimizer.param_groups[0]['lr']}")
            print(f"当前错误过滤阈值: {self.error_filter_threshold:.2f}, 总训练回合数: {self.total_episodes}")

    def is_finish(self, current_raw_state): # state_py_list from C++
        """
        检查 episode 是否应该终止。
        现在会包含来自网络的终止信号预测。
        返回 False 表示终止，True 表示继续。
        """
        # 1. 检查是否达到最大步数
        if self.current_episode_steps >= self.max_steps_per_episode:
            print(f"Episode finished: Reached max steps ({self.current_episode_steps}/{self.max_steps_per_episode}).")
            return False 

        # 2. 检查累计奖励是否达到目标阈值
        if self.current_episode_reward >= self.reward_threshold_to_stop:
            print(f"Episode finished: Cumulative reward ({self.current_episode_reward:.2f}) reached threshold ({self.reward_threshold_to_stop:.2f}).")
            return False

        # 3. 检查是否长时间没有显著改善
        if self.steps_since_last_improvement >= self.no_improvement_steps_threshold:
            print(f"Episode finished: No significant reward improvement for {self.steps_since_last_improvement} steps (threshold: {self.no_improvement_steps_threshold}).")
            return False

        # 4. 检查是否已无合法动作（所有 sheet 能量 <= 0）
        if not current_raw_state or all(row[1] <= 0 for row in current_raw_state): 
             print(f"Episode finished: No sheets with positive energy remaining or state is empty.")
             return False
        
        # --- MODIFIED SECTION START ---
        # 5. Query the network for its 'done' prediction for the current state
        if self.use_geometric_features : # and isinstance(self.policy_net, network.AttentionHexMeshQNet)
            # Prepare state for network input
            features, _ = self.extract_features(current_raw_state)
            num_nodes = len(features)
            if num_nodes == 0: # If state becomes empty after extraction, consider it done.
                print(f"Episode finished: State became empty after feature extraction.")
                return False

            x = torch.tensor(features, dtype=torch.float32).to(self.device)
            edge_list = [[i, j] for i in range(num_nodes) for j in range(num_nodes) if i != j]
            edge_index = torch.tensor(edge_list, dtype=torch.long).t().contiguous().to(self.device) if edge_list else torch.empty((2,0), dtype=torch.long).to(self.device)
            batch_tensor = torch.zeros(num_nodes, dtype=torch.long).to(self.device)
            sheet_node_idx_tensor = [torch.tensor([i]).to(self.device) for i in range(num_nodes)]
            sheet_features_tensor = x # Assuming sheet_features are same as x for AttentionHexMeshQNet

            with torch.no_grad():
                _, state_done_logit = self.policy_net(x, edge_index, batch_tensor, sheet_node_idx_tensor, sheet_features_tensor)
            
            # The state_done_logit is for the entire state. If batch_size in training is >1, 
            # g_emb is [B, H] and state_done_logit is [B]. Here, for a single state eval, it should be [1] or scalar.
            done_prob = torch.sigmoid(state_done_logit.squeeze()).item() # Ensure it's a scalar probability

            # Define a threshold for the network's done signal
            network_done_threshold = 0.5 # This can be tuned
            if done_prob > network_done_threshold:
                print(f"Episode finished: Network predicted done (prob: {done_prob:.4f} > {network_done_threshold}).")
                return False # Terminate if network says so
        # else: non-geometric path - if U_linear_QNet is also modified, add similar logic here.
        
        # --- MODIFIED SECTION END ---
        return True

    def find_sheet_id(self, state, action):
        for i, row in enumerate(state):
            if row[0] == action:
                return i
        return -1

    def print_state(self, state):
        print("state:")
        for row in state:
            print(f"{row[0]} {row[1]}")
