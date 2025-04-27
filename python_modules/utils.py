import random
import pickle
import numpy as np

#memory buffer
class ReplayMemory:
    def __init__(self, capacity):
        self.capacity = capacity
        self.memory = []
        self.priorities = []
        self.use_priorities = True  # 是否启用优先级采样
        self.alpha = 0.6  # 优先级权重，控制采样偏差程度
        self.beta = 0.4  # 初始重要性采样权重，会随训练进度增加
        self.beta_increment = 0.001  # beta增长率
        self.epsilon = 0.01  # 为避免零优先级
        
    def push(self, event, priority=None):
        """添加一个事件到记忆中，可选指定优先级"""
        # 事件结构: (state, action, next_state, reward, is_error)
        # is_error表示是否为非法动作错误
        
        if priority is None:
            # 默认优先级基于事件的第4个元素(reward)和第5个元素(is_error)
            if len(event) >= 5 and event[4]:  # 如果是错误选择
                priority = 2.0  # 给予高优先级
            elif len(event) >= 4:
                # 对于负奖励，给予较高优先级
                priority = 1.0 + abs(min(event[3], 0)) * 0.1
            else:
                priority = 1.0
        
        if len(self.memory) >= self.capacity:
            # 如果到达容量，替换最低优先级的条目
            if self.use_priorities:
                idx = np.argmin(self.priorities)
                self.memory[idx] = event
                self.priorities[idx] = priority
            else:
                # 如果不使用优先级，则遵循原有的FIFO策略
                del self.memory[0]
                self.memory.append(event)
                if self.use_priorities:
                    del self.priorities[0]
                    self.priorities.append(priority)
        else:
            # 容量未满，直接添加
            self.memory.append(event)
            if self.use_priorities:
                self.priorities.append(priority)

    def sample(self, batch_size):
        """根据优先级采样记忆"""
        if self.use_priorities and len(self.priorities) > 0:
            # 增加beta随训练进度
            self.beta = min(1.0, self.beta + self.beta_increment)
            
            # 计算采样概率
            probs = np.array(self.priorities) ** self.alpha
            probs /= probs.sum()
            
            # 加权随机采样
            indices = np.random.choice(len(self.memory), batch_size, p=probs)
            samples = [self.memory[idx] for idx in indices]
            
            # 计算重要性采样权重
            weights = (len(self.memory) * probs[indices]) ** (-self.beta)
            weights /= weights.max()  # 归一化权重
            
            return samples, weights, indices
        else:
            # 如果不使用优先级或没有优先级信息，使用原始随机采样
            return random.sample(self.memory, batch_size), None, None

    def update_priorities(self, indices, errors):
        """更新采样记忆的优先级"""
        if self.use_priorities and indices is not None:
            for idx, error in zip(indices, errors):
                self.priorities[idx] = error + self.epsilon

    def __len__(self):
        return len(self.memory)
    
    def save(self, path):
        """保存记忆到文件"""
        with open(path, 'wb') as f:
            pickle.dump({
                'memory': self.memory,
                'priorities': self.priorities if self.use_priorities else []
            }, f)
        print(f"记忆已保存至 {path}")
    
    def load(self, path):
        """从文件加载记忆"""
        try:
            with open(path, 'rb') as f:
                data = pickle.load(f)
                
                # 支持新旧格式
                if isinstance(data, dict):
                    self.memory = data.get('memory', [])
                    self.priorities = data.get('priorities', [])
                    # 如果加载的优先级不完整，重新初始化
                    if len(self.priorities) != len(self.memory):
                        self.priorities = [1.0] * len(self.memory)
                else:
                    # 旧格式，仅有记忆列表
                    self.memory = data
                    self.priorities = [1.0] * len(self.memory)
                    
            print(f"从 {path} 加载了 {len(self.memory)} 条记忆")
            return True
        except FileNotFoundError:
            print(f"文件 {path} 不存在")
            return False
        except Exception as e:
            print(f"加载记忆时发生错误: {e}")
            return False