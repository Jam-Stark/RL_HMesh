#pragma once
#include <torch/torch.h>
#include <vector>
#include <deque>
#include <random>

constexpr int action_scale = 100;
constexpr int state_scale[] = { 100,2 };

constexpr int dim_state = 2;
constexpr int dim_action = 1;

constexpr int CAPACITY = 100000;

constexpr int size1 = 128;
constexpr int size2 = 256;
constexpr int size3 = 512;

constexpr int BATCH_SIZE = 32;

constexpr int MODEL_TYPE[] = { 1,0 };
extern std::vector<int> STEP_INERVAL;



// 定义经验结构体
struct Transition {
    std::vector<std::vector<double>> state;
    int action;
    std::vector<std::vector<double>> next_state;
    double reward;

    Transition(const std::vector<std::vector<double>>& s, int a, const std::vector<std::vector<double>>& ns, double r)
        : state(s), action(a), next_state(ns), reward(r) {}
    bool is_next_state_valid() const {
        // 这里假设状态是一个 std::vector<double>
        return !state.empty(); // 如果状态非空，则返回 true；否则返回 false
    }
    bool is_next_state_valid2() const {
        //当前state中的所有第二列元素都>0，返回true
        for (int i = 0; i < next_state.size(); i++)
        {
            if (next_state[i][1] > 0) {
                std::cout << "is_next_state_valid2: valid" << std::endl;
                return true;
            }
        }
        std::cout << "is_next_state_valid2: invalid" << std::endl;
        return false;
    }
};

// 定义经验回放缓冲区类
class ReplayMemory {
public:
    ReplayMemory(int capacity) : capacity(capacity) {}
    ReplayMemory() :capacity(0) {}
    // 从缓冲区中随机采样一个批次的经验
    std::vector<Transition> sample(int batch_size) {
        std::vector<Transition> batch;
        std::sample(memory.begin(), memory.end(), std::back_inserter(batch), batch_size, std::mt19937{std::random_device{}()});
        return batch;
    }

    // 将经验添加到缓冲区中
    void push(const Transition& transition) {
        memory.push_back(transition);
        if (memory.size() > capacity) {
            memory.pop_front();
        }
    }

    int size() const {
        return memory.size();
    }

private:
    int capacity;
    std::deque<Transition> memory;
};

//继承torch::optim::LRScheduler类，实现自定义的学习率调度器
class MyLRScheduler : public torch::optim::LRScheduler {
public:
    MyLRScheduler(torch::optim::Optimizer& optimizer, std::vector<int> _step_interval, double _gamma) : torch::optim::LRScheduler(optimizer), step_interval(_step_interval), gamma(_gamma) {}

protected:
    virtual std::vector<double> get_lrs() override {
        // 获取当前的迭代次数
        unsigned current_step = step_count_;

        // 根据当前的迭代次数动态计算学习率
        double lr = this->get_current_lrs()[0];  // 初始学习率,假设只有一个参数

        for (int i = 0; i < step_interval.size(); i++) {
            if (current_step < step_interval[i]) {
                break;
            }
            lr *= gamma;
        }

        return { lr };
    }

private:
    std::vector<int> step_interval;
    double gamma;
};

class GradCalcSwtitch //（废弃）
{
public:
    GradCalcSwtitch() {
        torch::NoGradGuard();
    };
    ~GradCalcSwtitch() {};
};

