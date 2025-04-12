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
from utils import ReplayMemory
import network

class Agent:
    def __init__(self, eta=0.5, gamma=0.99, capacity=100000, batch_size=32, episode=0):
        self.rewards = []
        self.eta = eta
        self.gamma = gamma
        self.capacity = capacity
        self.batch_size = batch_size
        self.episode = episode
        self.memory = ReplayMemory(self.capacity)
        self.model =  network.U_linear_QNet(2, 1, 10, 10, 10)
        self.optimizer = optim.Adam(self.model.parameters(), lr=0.1)
        self.scheduler = optim.lr_scheduler.StepLR(self.optimizer, step_size=1, gamma=0.1)

    def save_model(self, path):
        torch.save(self.model.state_dict(), path)

    def remember(self, state, action, next_state, reward):
        self.memory.push((state, action, next_state, reward))

    def choose_action(self, state, episode):
        sample = random.random()
        sample = 0
        eps_threshold = 0.6 / (episode + 1)
        if sample < eps_threshold:
            action = random.randint(0, len(state) - 1)
            print(f"random action: {action} corresponding sheet id: {state[action][0]} sheet energy: {state[action][1]}")

            #choose a action that has energy >0 to test
            while state[action][1] <= 0:
                action = random.randint(0, len(state) - 1)
                print(f"random action: {action} corresponding sheet id: {state[action][0]} sheet energy: {state[action][1]}")
            return action
        else:
            with torch.no_grad():
                state_tensor = torch.tensor(state, dtype=torch.float32)  # shape: (n, 2)
                scores = self.model(state_tensor)  # shape: (n, 1)
                action = torch.argmax(scores.view(-1)).item()
                print(f"greedy action: {action} corresponding sheet id: {state[action][0]} sheet energy: {state[action][1]}")
                return action

    def get_action(self, state): # 推理使用
        with torch.no_grad():
            state_tensor = torch.tensor(state, dtype=torch.float32)
            scores = self.model(state_tensor)
            action = torch.argmax(scores.view(-1)).item()
            return action

    def replay(self):
        if len(self.memory) < self.batch_size:
            return
        transitions = self.memory.sample(self.batch_size)
        batch = list(zip(*transitions))
        state_batch = torch.tensor(batch[0], dtype=torch.float32)
        action_batch = torch.tensor(batch[1], dtype=torch.long)
        next_state_batch = torch.tensor(batch[2], dtype=torch.float32)
        reward_batch = torch.tensor(batch[3], dtype=torch.float32)
        print(f"state_batch: {state_batch}")
        print(f"action_batch: {action_batch}")
        print(f"next_state_batch: {next_state_batch}")
        print(f"reward_batch: {reward_batch}")
        self.optimizer.zero_grad()
        # 修改：去掉最后一维，形状变为 (batch, n)
        q_values = self.model(state_batch).squeeze(-1)
        next_q_values = self.model(next_state_batch).squeeze(-1)
        q_value = q_values.gather(1, action_batch.unsqueeze(1)).squeeze(1)
        next_q_value = next_q_values.max(1)[0].detach()
        expected_q_value = reward_batch + self.gamma * next_q_value
        loss = F.smooth_l1_loss(q_value, expected_q_value)
        loss.backward()
        self.optimizer.step()
        self.scheduler.step()

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
