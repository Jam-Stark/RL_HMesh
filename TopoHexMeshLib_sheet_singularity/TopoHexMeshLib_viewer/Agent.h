#pragma once
#include <iostream>
#include <fstream>
#include <random>
#include <torch/torch.h>
#include <torch/script.h>
#include "settings.h"

#include "Model.h"


class Agent
{
	/*定义智能体，包括其policy,action,reward function*/
public:
	/*
	设定超参数，网络参数，memory
	*/
	std::vector<double> rewards;
	double eta, gamma;
	int capacity, batch_size, episode, num;
	//memory包含 state action next_state reward
	ReplayMemory memory;
	Linear_QNet* model;

	torch::optim::Adam* optimizer;
	//torch::optim::LRScheduler scheduler; LRScheduler是一个抽象类
	MyLRScheduler* scheduler;


	Agent(double _eta = 0.5, double _gamma = 0.99, int _capacity = 100000, int _batch_size = BATCH_SIZE, int _episode = 0);
	~Agent();

	void save_model(std::string path)
	{
		//torch::save(model, path);
		model->save(path);
	}

	void remember(std::vector < std::vector<double> > state, int action, std::vector<std::vector<double>> next_state, double reward);

	int choose_action(std::vector < std::vector<double> > state, int _episode);

	int get_action(std::vector<std::vector<double>> state);

	//void train_long_memory();
 //   void train_short_memory();

	void replay();

	/*--tool--*/
	bool is_finish(std::vector<std::vector<double>> state)
	{
		//当前state中的所有第二列元素都<=0，返回true
		for (int i = 0; i < state.size(); i++)
		{
			if (state[i][1] > 0)
				return false;
		}
		return true;
	}
	int find_sheet_id(std::vector<std::vector<double>> state, int action)
	{
		for (int i = 0; i < state.size(); i++)
		{
			if (state[i][0] == action)
			{
				return i;
			}
		}
		return -1;
	}
	void print_state(std::vector<std::vector<double>> state)
	{
		std::cout << "state: " << std::endl;
		for (std::vector<std::vector<double>>::iterator i = state.begin(); i != state.end(); i++)
		{
			std::cout << i->data()[0] << " " << i->data()[1] << std::endl;
		}
	}
	//将vector状态转换为张量
	torch::Tensor convertToTensor2(const std::vector<std::vector<double>>& state) {
		int rows = state.size(); // 外部向量的大小，即行数
		int cols = state[0].size(); // 假设每个内部向量的长度相同，即列数

		// 动态分配一个连续的数组来模拟二维数组
		auto* state_array = new float[rows * cols]; // 创建一个足够大的数组

		// 填充数组
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				state_array[i * cols + j] = static_cast<float>(state[i][j]);
			}
		}

		// 使用 from_blob 创建 Tensor，注意指定尺寸
		torch::Tensor tensor = torch::from_blob(state_array, { rows, cols }, torch::kFloat32).to(torch::kCUDA);

		// 因为 from_blob 不拷贝数据，所以需要 clone 来创建数据副本
		torch::Tensor tensor_clone = tensor.clone().to(torch::kCUDA);

		// 删除原始数组以释放内存
		delete[] state_array;

		return tensor_clone; // 返回克隆后的张量，现在它拥有自己的数据副本
	}
	//vector转tensor
	template<typename T>
	torch::Tensor convertToTensor(const std::vector<T>& vec) {
		int size = vec.size(); // 外部向量的大小

		// 动态分配一个连续的数组来模拟二维数组
		T* _array = new T[size]; // 创建一个足够大的数组

		// 填充数组
		for (int i = 0; i < size; i++) {
				_array[i ] = vec[i];
			
		}

		if(T==double)
		torch::Tensor tensor = torch::from_blob(_array, { size }, torch::kFloat32);
		else
		torch::Tensor tensor = torch::from_blob(_array, { size }, torch::kInt64);
		// 因为 from_blob 不拷贝数据，所以需要 clone 来创建数据副本
		torch::Tensor tensor_clone = tensor.clone();

		// 删除原始数组以释放内存
		delete[] _array;

		return tensor_clone; // 返回克隆后的张量，现在它拥有自己的数据副本
	}

};