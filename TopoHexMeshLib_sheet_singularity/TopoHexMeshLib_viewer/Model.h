#pragma once
#include <torch/torch.h>
#include <torch/script.h>
#include <iostream>
#include <vector>

#include "settings.h"

// Define the equivalent of Linear_QNet in C++
class Linear_QNet :public torch::nn::Module {
    torch::nn::Linear fc1 = nullptr, fc2 = nullptr, fc3 = nullptr, fc4 = nullptr;


public:
    Linear_QNet() {
        //torch::Device device(torch::kCUDA);
        fc1 = register_module("fc1", torch::nn::Linear(state_scale[1], size1));
        fc2 = register_module("fc2", torch::nn::Linear(size1, size1));
        fc3 = register_module("fc3", torch::nn::Linear(size1, dim_action));
        fc4 = register_module("fc4", torch::nn::Linear(state_scale[1], dim_action));
        // Move all model parameters to GPU if available
        if (torch::cuda::is_available()) {
            this->to(torch::kCUDA);
            std::cout << "Model moved to GPU" << std::endl;
        }
    }

    torch::Tensor forward(torch::Tensor x) {
        //torch::Device device(torch::kCUDA);
        //torch::Tensor new_x=x.clone();
        //x = x.to(device);
        //}//将tensor转移到GPU上
        torch::Tensor identity = fc4->forward(x).to(torch::kCUDA);
        x = torch::relu(fc1->forward(x));
        x = torch::relu(fc2->forward(x));
        x = fc3->forward(x);

        return x + identity;
    }

    void save(std::string file_name = "../data/model/model.pt") {
        torch::serialize::OutputArchive output_archive;
        output_archive.save_to(file_name);
    }
};

class U_linear_QNet :public torch::nn::Module {
    torch::nn::Linear up1 = nullptr, up2 = nullptr, up3 = nullptr, fc = nullptr, down1 = nullptr,
        down2 = nullptr, down3 = nullptr, skip1 = nullptr, skip2 = nullptr, skip3 = nullptr;

public:
    U_linear_QNet()
    {
        up1 = register_module("up1", torch::nn::Linear(dim_state, size1));
        up2 = register_module("up2", torch::nn::Linear(size1, size2));
        up3 = register_module("up3", torch::nn::Linear(size2, size3));

        fc = register_module("fc", torch::nn::Linear(size3, size3));

        down1 = register_module("down1", torch::nn::Linear(size3, size2));
        down2 = register_module("down2", torch::nn::Linear(size2, size1));
        down3 = register_module("down3", torch::nn::Linear(size1, dim_action));

        skip1 = register_module("skip1", torch::nn::Linear(dim_state, dim_action));
        skip2 = register_module("skip2", torch::nn::Linear(size1, size1));
        skip3 = register_module("skip3", torch::nn::Linear(size2, size2));
    }

    torch::Tensor forward(torch::Tensor x) {
        //if (torch::cuda::is_available()) {
        //    torch::Device device(torch::kCUDA);
        //    x = x.to(device);
        //}//将tensor转移到GPU上
        torch::Tensor skip1_tensor = skip1->forward(x);
        torch::Tensor u1 = torch::relu(up1->forward(x));
        torch::Tensor skip2_tensor = skip2->forward(u1);
        torch::Tensor u2 = torch::relu(up2->forward(u1));
        torch::Tensor skip3_tensor = skip3->forward(u2);
        torch::Tensor u3 = torch::relu(up3->forward(u2));

        torch::Tensor d0 = fc->forward(u3);

        torch::Tensor d1 = torch::relu(down1->forward(d0));
        torch::Tensor d2 = torch::relu(down2->forward(d1 + skip3_tensor));
        torch::Tensor d3 = torch::relu(down3->forward(d2 + skip2_tensor));


        return d3 + skip1_tensor;
    }
};

// 设置优化器&计算梯度（废弃）
class QTrainer {
    torch::optim::Adam optimizer;
    torch::nn::MSELoss criterion;

public:
    QTrainer(Linear_QNet& model, float lr = 0.1, float gamma = 0.99)
        : optimizer(model.parameters(), torch::optim::AdamOptions(lr)),
        criterion(torch::nn::MSELossOptions()) {}

    QTrainer& operator=(const QTrainer& other) = default;

    void train_step_short(std::vector<double[2]> state, std::vector<int> action, double reward,
        std::vector<double[2]> next_state, bool done, Linear_QNet& model, float gamma = 0.9)
    {
        //turn to tensor
        torch::Tensor state_tensor = torch::from_blob(state.data(), { static_cast<long long>(state.size()),2 }, torch::kDouble).clone();//kDouble和kFloat64是同种数据类型
        torch::Tensor next_state_tensor = torch::from_blob(next_state.data(), { static_cast<long long>(next_state.size()),2 }, torch::kDouble).clone();
        torch::Tensor action_tensor = torch::tensor(action, torch::kInt).clone();
        torch::Tensor reward_tensor = torch::tensor(reward, torch::kDouble).clone();

        auto pred = model.forward(state_tensor);
        auto target = pred.clone();

        auto Q_new = reward_tensor.item<double>();
        if (!done)
            Q_new = reward_tensor.item<double>() + gamma * torch::max(model.forward(next_state_tensor)).item<double>();
        target[action_tensor.item<int>()] = Q_new;

        optimizer.zero_grad();
        auto loss = criterion(target, pred);
        loss.backward();
        optimizer.step();






    }
    void train_step_long(std::vector<std::vector<double[2]>> states, std::vector<std::vector<int>> actions, std::vector<double> rewards,
        std::vector<std::vector<double[2]>> next_states, std::vector<bool> dones, Linear_QNet& model, float gamma = 0.9)
    {
        // Convert states to tensor
        torch::Tensor states_tensor = torch::from_blob(states.data(), { static_cast<long long>(states.size()),static_cast<long long>(states[0].size()), 2 }, torch::kDouble);
        torch::Tensor actions_tensor = torch::from_blob(actions.data(), { static_cast<long long>(actions.size()), static_cast<long long>(actions[0].size()) }, torch::kInt);
        torch::Tensor rewards_tensor = torch::from_blob(rewards.data(), { static_cast<long long>(rewards.size()) }, torch::kDouble);
        torch::Tensor next_states_tensor = torch::from_blob(next_states.data(), { static_cast<long long>(next_states.size()),static_cast<long long>(next_states[0].size()), 2 }, torch::kDouble);

        auto pred = model.forward(states_tensor);
        auto target = pred.clone();

        for (int idx = 0; idx < dones.size(); ++idx) {
            auto Q_new = rewards_tensor[idx];
            if (!dones[idx]) {
                Q_new = rewards_tensor[idx] + gamma * torch::max(model.forward(next_states_tensor[idx])).item<float>();
            }
            target[idx][torch::argmax(actions_tensor[idx]).item<int>()] = Q_new;
        }

        optimizer.zero_grad();
        auto loss = criterion(target, pred);
        //auto loss = criterion.forward(target, pred);
        loss.backward();
        optimizer.step();
    }
};

