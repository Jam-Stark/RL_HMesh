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
	/*���������壬������policy,action,reward function*/
public:
	/*
	�趨�����������������memory
	*/
	std::vector<double> rewards;
	double eta, gamma;
	int capacity, batch_size, episode, num;
	//memory���� state action next_state reward
	ReplayMemory memory;
	Linear_QNet* model;

	torch::optim::Adam* optimizer;
	//torch::optim::LRScheduler scheduler; LRScheduler��һ��������
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
		//��ǰstate�е����еڶ���Ԫ�ض�<=0������true
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
	//��vector״̬ת��Ϊ����
	torch::Tensor convertToTensor2(const std::vector<std::vector<double>>& state) {
		int rows = state.size(); // �ⲿ�����Ĵ�С��������
		int cols = state[0].size(); // ����ÿ���ڲ������ĳ�����ͬ��������

		// ��̬����һ��������������ģ���ά����
		auto* state_array = new float[rows * cols]; // ����һ���㹻�������

		// �������
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				state_array[i * cols + j] = static_cast<float>(state[i][j]);
			}
		}

		// ʹ�� from_blob ���� Tensor��ע��ָ���ߴ�
		torch::Tensor tensor = torch::from_blob(state_array, { rows, cols }, torch::kFloat32).to(torch::kCUDA);

		// ��Ϊ from_blob ���������ݣ�������Ҫ clone ���������ݸ���
		torch::Tensor tensor_clone = tensor.clone().to(torch::kCUDA);

		// ɾ��ԭʼ�������ͷ��ڴ�
		delete[] state_array;

		return tensor_clone; // ���ؿ�¡���������������ӵ���Լ������ݸ���
	}
	//vectorתtensor
	template<typename T>
	torch::Tensor convertToTensor(const std::vector<T>& vec) {
		int size = vec.size(); // �ⲿ�����Ĵ�С

		// ��̬����һ��������������ģ���ά����
		T* _array = new T[size]; // ����һ���㹻�������

		// �������
		for (int i = 0; i < size; i++) {
				_array[i ] = vec[i];
			
		}

		if(T==double)
		torch::Tensor tensor = torch::from_blob(_array, { size }, torch::kFloat32);
		else
		torch::Tensor tensor = torch::from_blob(_array, { size }, torch::kInt64);
		// ��Ϊ from_blob ���������ݣ�������Ҫ clone ���������ݸ���
		torch::Tensor tensor_clone = tensor.clone();

		// ɾ��ԭʼ�������ͷ��ڴ�
		delete[] _array;

		return tensor_clone; // ���ؿ�¡���������������ӵ���Լ������ݸ���
	}

};