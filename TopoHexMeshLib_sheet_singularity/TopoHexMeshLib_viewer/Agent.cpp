#include "Agent.h"

Agent::Agent(double _eta, double _gamma, int _capacity, int _batch_size, int _episode)
{
	this->eta = _eta;
	this->gamma = _gamma;
	this->capacity = _capacity;

	this->batch_size = _batch_size;
	this->episode = _episode;

	this->model = new Linear_QNet();
	optimizer = new  torch::optim::Adam(model->parameters(), torch::optim::AdamOptions(0.1));
	scheduler = new MyLRScheduler(*optimizer, STEP_INERVAL, _gamma);

	this->memory = ReplayMemory(this->capacity);
	this->num = 0;

}

Agent::~Agent()
{
}

void Agent::remember(std::vector<std::vector<double>> state, int action, std::vector<std::vector<double>> next_state, double reward)
{
	auto t = Transition(state, action, next_state, reward);
	this->memory.push(t);
}

int Agent::choose_action(std::vector < std::vector<double> > state, int _episode)
{
	//print_state(state);
	if (is_finish(state))
		return -1;
	//std::cout << state << std::endl;
	this->episode = _episode;
	double eps = 0.6 * 1 / (1 + _episode);
	unsigned seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::default_random_engine RAND_GENERATER(seed);
	std::uniform_real_distribution<double> DISTRIBUSION(0.0, 1.0);
	double rand = DISTRIBUSION(RAND_GENERATER);
	std::cout << "rand:" << rand << " eps: " << eps << std::endl;
	if (rand < eps) {//Random action
		std::uniform_int_distribution<int> distributionInScale(0, state.size() - 1);
		int key = distributionInScale(RAND_GENERATER);
		//int key_id = find_sheet_id(state,key);
		while (state[key][1] <= 0) {
			key = distributionInScale(RAND_GENERATER);
			//key_id = find_sheet_id(state, key);
		}
		std::cout << "Random action: " << key << std::endl;
		return key;
	}
	else {
		//Greedy action
		this->model->eval();
		//����������ʵ��with�﷨�������ݶȼ���
		torch::NoGradGuard no_grad;
		//torch::Tensor state_tensor = torch::from_blob(state.data(), { static_cast<long long>(state.size()), static_cast<long long>(state[0].size()) }, torch::kFloat);
		torch::Tensor state_tensor = convertToTensor2(state);
		//state_tensor = state_tensor.view({ 1, -1 });
		//std::cout << "state_tensor:" << state_tensor << std::endl;
		torch::Tensor new_action = this->model->forward(state_tensor);
		//new_action = torch::softmax(new_action, 1);
		//new_action=new_action.view({ 1, -1 });
		//std::cout << "new_action: " << new_action << std::endl;
		for (int i = 0; i < state.size(); i++) {
			if (state[i][1] <= 0)
				new_action[i] = -1000000;//���Ƿ���action��ֵ����Ϊ�����ʹ�䲻��ѡ��
		}
		//std::cout << "masked_new_action: " << new_action << std::endl;
		new_action = torch::softmax(new_action, 0);
		//std::cout<<"softmaxed_new_action: "<<new_action<<std::endl;

		int action = new_action.argmax().item<int>();//ѡȡ���ʷֲ�������ֵ
		std::cout << "Greedy action:" << action << std::endl;
		return action;
	}

}

int Agent::get_action(std::vector < std::vector<double> > state)
{
		//Greedy action
		this->model->eval();
		//����������ʵ��with�﷨�������ݶȼ���
		torch::NoGradGuard no_grad;
		//torch::Tensor state_tensor = torch::from_blob(state.data(), { static_cast<long long>(state.size()), static_cast<long long>(state[0].size()) }, torch::kFloat);
		torch::Tensor state_tensor = convertToTensor2(state);
		//state_tensor = state_tensor.view({ 1, -1 });
		//std::cout << "state_tensor:" << state_tensor << std::endl;
		torch::Tensor new_action = this->model->forward(state_tensor);
		//new_action = torch::softmax(new_action, 1);
		//new_action=new_action.view({ 1, -1 });
		//std::cout << "new_action: " << new_action << std::endl;
		for (int i = 0; i < state.size(); i++) {
			if (state[i][1] <= 0)
				new_action[i] = -1000000;//���Ƿ���action��ֵ����Ϊ�����ʹ�䲻��ѡ��
		}
		//std::cout << "masked_new_action: " << new_action << std::endl;
		new_action = torch::softmax(new_action, 0);
		//std::cout<<"softmaxed_new_action: "<<new_action<<std::endl;

		int action = new_action.argmax().item<int>();//ѡȡ���ʷֲ�������ֵ
		std::cout << "Greedy action:" << action << std::endl;
		return action;

}

void Agent::replay()
{
	if (this->memory.size() < this->batch_size)
		return;
	else {
		std::vector<Transition> sample = this->memory.sample(this->batch_size);
		// ��ʼ����������
		std::vector<torch::Tensor> state_tensors;
		std::vector<int> action_tensors;
		std::vector<double> reward_tensors;
		std::vector<torch::Tensor> next_state_tensors;

		// ���������е�ÿ������
		int sample_i = 0;
		for (auto& transition : sample) {
			// �������е�����ת��Ϊ��������ӵ���Ӧ������������
			state_tensors.push_back(convertToTensor2(transition.state));
			//state_tensors.push_back(torch::from_blob(transition.state.data(), { static_cast<long long>(transition.state.size()), static_cast<long long>(transition.state[0].size()) }, torch::kFloat).clone());
			action_tensors.push_back(transition.action);
			reward_tensors.push_back(transition.reward);
			if (!is_finish(transition.next_state))
				next_state_tensors.push_back(convertToTensor2(transition.next_state));

			//std::cout<< "state "<< sample_i << " : " << std::endl;
			sample_i++;
			//print_state(transition.next_state);
		}
		
		// ����������ת��Ϊ����
		//std::cout << "state_tensors: " << state_tensors << std::endl;
		torch::Tensor state_batch = torch::stack(state_tensors).to(torch::kCUDA);
		//std::cout << "state_batch: " << state_batch << std::endl;
		
		//torch::Tensor action_batch = torch::tensor(action_tensors, torch::kInt64).unsqueeze(1).unsqueeze(1);
		//std::cout << "action_tensors: " << action_tensors << std::endl;
		torch::Tensor action_batch = torch::tensor(action_tensors, torch::kInt64).unsqueeze(1).unsqueeze(2).to(torch::kCUDA);
		//std::cout<< "action_batch: " << action_batch << std::endl;

		//torch::Tensor reward_batch = torch::tensor(reward_tensors, torch::kFloat).unsqueeze(1);
		//std::cout << "reward_tensors: " << reward_tensors << std::endl;
		torch::Tensor reward_batch = torch::tensor(reward_tensors, torch::kFloat).unsqueeze(1).to(torch::kCUDA);
		//std::cout << "reward_batch: " << reward_batch << std::endl;

		torch::Tensor non_final_next_states_batch = torch::stack(next_state_tensors).to(torch::kCUDA);
		//std::cout << "non_final_next_states_batch: " << non_final_next_states_batch << std::endl;

		//����train model
		this->model->train();
		//pred
		auto outputs = this->model->forward(state_batch);
		//std::cout << "outputs: " << outputs << std::endl;

		auto state_action_values = torch::gather(outputs, 1, action_batch);
		//std::cout << "state_action_values: " << state_action_values  << std::endl;

		//gt
		std::vector<int> _non_final_mask(sample.size());//std::vector<bool> ��һ��������������Բ���ֵ���н��յĴ洢��ͨ��ʹ��λѹ������ʾÿ������ֵ��������˵��std::vector<bool> �е�ÿ������ֵ����ֱ�Ӵ洢��һ���ֽ���,����ʹ��int���
		std::transform(sample.begin(), sample.end(), _non_final_mask.begin(), [](const Transition& t) { return t.is_next_state_valid2(); });//��sample�е�ÿ��Ԫ��ת��Ϊboolֵ
		//for (auto& transition : sample) {
		//	if (transition.is_next_state_valid2())
		//		_non_final_mask.push_back(1);
		//	else
		//		_non_final_mask.push_back(0);
		//}
		//std::cout << "_non_final_mask: " << _non_final_mask << std::endl;
		torch::Tensor non_final_mask = torch::from_blob(_non_final_mask.data(), { static_cast<long long>(_non_final_mask.size()) }, torch::kInt).to(torch::kCUDA);
		//std::cout << "non_final_mask: " << non_final_mask << std::endl;//[ CPUFloatType{BATCH_SIZE} ]
		torch::Tensor next_state_values = torch::zeros({ this->batch_size }).to(torch::kCUDA);
		//std::cout << "next_state_values: " << next_state_values << std::endl;//[ CPUFloatType{BATCH_SIZE} ]

		//next_state_tensors[non_final_mask] = this->model->forward(non_final_next_states_batch).max(1)[0].detach();
		// ����ģ�͵��������ȡ���ֵ
		torch::Tensor model_output = model->forward(non_final_next_states_batch).to(torch::kCUDA);
		//std::cout << "model_output: " << model_output << std::endl;

		// ����ÿ��������Ԫ�ص����ֵ
		auto max_results = model_output.max(1);
		//std::cout << "max_results: " << max_results << std::endl;

		// ��ȡ���ֵ��������״Ϊ{BATCH_SIZE, 1, 1}
		torch::Tensor max_values = std::get<0>(max_results);

		// ȥ�����һ��ά�ȣ���״��Ϊ{32, 1}
		//max_values = max_values.squeeze(-1); // ʹ��squeezeȥ����СΪ1��ά��
		//std::cout << "max_values: " << max_values << std::endl;

		//max_values �����{32,1}����non_final_mask��{32}��������Ҫ��max_values��ά�ȵ���Ϊ{32}
		//max_values = max_values.view({ -1 });

		// ���� non_final_mask �е�ֵ������ next_state_values
		//non_final_mask = non_final_mask.view({ -1 });
		next_state_values = torch::zeros_like(non_final_mask, torch::kFloat).unsqueeze(1); // ��ʼ��Ϊ��
		//std::cout << "next_state_values: " << next_state_values << std::endl;   //[ CPUFloatType{32} ]

		// ʹ����Ԫ�������������mask��ģ��Ԥ�⶯����Qֵ���� next_state_values
		next_state_values.index_put_({ non_final_mask.to(torch::kBool) }, max_values.detach());
		//std::cout << "calced next_state_values: " << next_state_values << std::endl;
		
		torch::Tensor expected_state_action_values = (reward_batch + next_state_values * this->gamma).to(torch::kCUDA);
		//std::cout << "expected_state_action_values: " << expected_state_action_values << std::endl;
		this->num++;
		state_action_values = state_action_values.squeeze(-1);
		std::cout << "new state_action_values: " << state_action_values << std::endl;
		// ������ʧ
		auto loss = torch::nn::functional::smooth_l1_loss(state_action_values, expected_state_action_values);

		this->optimizer->zero_grad();
		loss.backward();
		this->optimizer->step();
		this->scheduler->step();
	}
}
