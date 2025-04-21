#ifndef STATE_FUNCTIONS_H
#define STATE_FUNCTIONS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <vector>
#include <string>
#include <ctime>
#include "tools.h"
#include "meshQuality.h"
#include "sheet_operation.h"
#include "singularity_number.h"
#include "topoMesh.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

using namespace HMeshLib;

// State类定义
class State {
	public:
		std::vector<int> sheet_id;
		std::vector<double> sheet_energy;
        int sheet_num;
        // 默认构造函数
        State() {
            sheet_num = 0;
        }
		
		// 清空 state 数据
		void clear() {
			sheet_id.clear();
			sheet_energy.clear();
            sheet_num = 0;
		}
		// 返回状态大小
		size_t size() const {
			return sheet_id.size();
		}
		// 添加一条记录
		void add(int id, double energy) {
			sheet_id.push_back(id);
			sheet_energy.push_back(energy);
            sheet_num++;
		}
		// 重载赋值操作符
		State& operator=(const State& other) {
			if (this != &other) {
				sheet_id = other.sheet_id;
				sheet_energy = other.sheet_energy;
			}
			return *this;
		}
		
		//根据sheet_id查找sheet_energy
		double get_energy(int id) {
			for (size_t i = 0; i < sheet_id.size(); i++) {
				if (sheet_id[i] == id) {
					return sheet_energy[i];
				}
			}
			return 0.0;
		}

        // 检查state中除了当前action外是否还有其他sheet_energy > 0的情况, 有则返回true,此时为异常
        bool check_error_choice(int action) {
            for (size_t i = 0; i < sheet_id.size(); i++) {
                if (sheet_id[i] != action && sheet_energy[i] <= 0) {
                    return true;
                }
            }
            return false;
        }

        // check if alraedy have the sheet_id
        bool check_sheet_id(int id) {
            for (size_t i = 0; i < sheet_id.size(); i++) {
                if (sheet_id[i] == id) {
                    return true;
                }
            }
            return false;
        }

        void print() {
            if(sheet_id.size() == 0) {
                std::cout << "state is empty" << std::endl;
            }
            else {
                std::cout << "state: " << std::endl;
                for (size_t i = 0; i < sheet_id.size(); i++) {
                    std::cout << "sheet id: " << sheet_id[i] << " energy: " << sheet_energy[i] << std::endl;
                }
            }
        }
		
	};

// ────────────────────
std::string format_matrix(const State& state);
void log(int episode, const State& state, int action, const State& next_state, float reward);
void fill_state(State& state, int origin_size);
int find_sheet_id(State state, int action);
bool check_in_state(State state, int action);
std::vector<TE*> get_sheet_byId(TMesh* mesh, int sheet_id, sheet_operation<TMesh>& sheet_op);
void calc_state(TMesh* mesh, State& state, sheet_operation<TMesh>& sheet_op);
int play_action(int action, int done, State& state, TMesh& tmesh, sheet_operation<TMesh>& sheet_op, get_singularity_number<TMesh>& get_singularity_num_op);
float calc_reward(int current_singularity_num, get_singularity_number<TMesh> get_singularity_num_op);
void print_state(const State& state);
void compare_diff(const State& state1, const State& state2);
py::list state_to_list(const State& state);
// ────────────────────

std::string format_matrix(const State& state) {
    std::stringstream ss;
    for (size_t i = 0; i < state.size(); i++) {
        if (state.sheet_energy[i] > 0) {
            ss << state.sheet_id[i] << " " << state.sheet_energy[i] << "\n";
        }
    }
    return ss.str();
}

// 修改训练日志保存位置
void log(int episode, const State& state, int action, const State& next_state, float reward) {
    // 从环境变量获取会话ID
    char* session_id_env = nullptr;
    size_t len = 0;
    _dupenv_s(&session_id_env, &len, "RL_HMESH_SESSION_ID");
    std::string session_id = session_id_env ? std::string(session_id_env) : "unknown";
    if (session_id_env) free(session_id_env);
    
    // 创建当前会话的日志目录
    std::string log_dir = "f:\\RL_HMesh\\logs\\" + session_id;
    std::string log_filename = log_dir + "\\training.log";
    
    std::ofstream log_file(log_filename, std::ios::app);
    
    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);
    std::tm now_tm;
    localtime_s(&now_tm, &now_time);
    
    std::stringstream ss;
    ss << std::put_time(&now_tm, "%Y-%m-%d %H:%M:%S");
    std::string time_str = ss.str();
    
    if (log_file) {
        log_file << "Training Time: " << time_str << "\n"
                 << "Episode: " << episode << "\n"
                 << "State:\n" << format_matrix(state)
                 << "Action: [" << state.sheet_id[action] << ", " << state.sheet_energy[action] << "]\n"
                 << "Next State:\n" << format_matrix(next_state)
                 << "Reward: " << reward << "\n"
                 << "---------------------------------\n";
    } else {
        std::cerr << "Error opening log file: " << log_filename << "\n";
    }
    
    log_file.close();
}

void fill_state(State& state, int origin_size) {
    if (state.size() < static_cast<size_t>(origin_size)) {
        for (int i = 1; i <= origin_size; i++) {
            bool key = true;
            for (size_t j = 0; j < state.size(); j++) {
                if (state.sheet_id[j] == i) {
                    key = false;
                    break;
                }
            }
            if (key) {
                state.add(i, 0);
            }
        }
    } else {
        std::cout << "state is full" << std::endl;
    }
}

int find_sheet_id(State state, int action) {
    for (size_t i = 0; i < state.size(); i++) {
        if (state.sheet_id[i] == action) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

bool check_in_state(State state, int action) {
    for (size_t i = 0; i < state.size(); i++) {
        if (state.sheet_id[i] == action) {
            return true;
        }
    }
    return false;
}

std::vector<TE*> get_sheet_byId(TMesh* mesh, int sheet_id, sheet_operation<TMesh>& sheet_op) {
    std::cout << "sheet_id: " << sheet_id << std::endl;
    for (TMesh::MEIterator eite(mesh); !eite.end(); eite++) {
        TE* be = *eite;
        if (be->sheet() == sheet_id) {
            std::cout << "find the sheet" << std::endl;
            std::vector<TE*> sheet = sheet_op.get_one_sheet(be);
            std::cout << "find sheet done " << std::endl;
            return sheet;
        }
    }
    std::cout << "can't find the sheet" << std::endl;
    return {};
}

void calc_state(TMesh* mesh, State& state, sheet_operation<TMesh>& sheet_op) {
    // 获取最大的sheet ID
    int max_sheet_id = 0;
    for (TMesh::MEIterator eite(mesh); !eite.end(); eite++) {
        TE* be = *eite;
        int sheet_id = be->sheet();
        if (sheet_id > max_sheet_id) {
            max_sheet_id = sheet_id;
        }
    }
    std::cout << "max_sheet_id: " << max_sheet_id << std::endl;

    state.print();
    // 按顺序处理每个sheet ID，确保所有sheet都被添加到state中
    for (int sheet_id = 1; sheet_id <= max_sheet_id; sheet_id++) {
        // 如果这个sheet ID已经在state中，跳过
        if (state.check_sheet_id(sheet_id)){
            std::cout << "sheet_id: " << sheet_id << " already in state" << std::endl;
            continue;
        }
            
        // 找到具有这个sheet ID的边
        TE* edge_with_id = nullptr;
        for (TMesh::MEIterator eite(mesh); !eite.end(); eite++) {
            TE* be = *eite;
            if (be->sheet() == sheet_id) {
                edge_with_id = be;
                std::cout << "find edge with sheet_id: " << sheet_id << std::endl;
                break;
            }
        }
        
        // 如果找不到具有这个sheet ID的边，跳过
        if (edge_with_id == nullptr){
            //std::cout << "can't find edge with sheet_id: " << sheet_id << std::endl;
            continue;
        }
            
        // 获取sheet并计算能量
        std::vector<TE*> sheet = sheet_op.get_one_sheet(edge_with_id);
        double sheet_energy = sheet_op.predict_sheet_collapse_energy(sheet);
        std::cout << "sheet_id: " << sheet_id << " sheet_energy: " << sheet_energy << std::endl;
        state.add(sheet_id, sheet_energy);
    }
    std::cout << "state size: " << state.sheet_num << std::endl;
}

// 修改play_action函数，整合日志
int play_action(int action, int done, State& state, TMesh& tmesh, sheet_operation<TMesh>& sheet_op, get_singularity_number<TMesh>& get_singularity_num_op) {
    // 不再使用固定的state size
    std::cout << "action: " << action << " state_energy: " << state.sheet_energy[action] << std::endl;
    
    // 从环境变量获取会话ID
    char* session_id_env = nullptr;
    size_t len = 0;
    _dupenv_s(&session_id_env, &len, "RL_HMESH_SESSION_ID");
    std::string session_id = session_id_env ? std::string(session_id_env) : "unknown";
    if (session_id_env) free(session_id_env);
    
    // 创建当前会话的日志目录
    std::string log_dir = "f:\\RL_HMesh\\logs\\" + session_id;
    
    // 创建特定于此操作的日志文件
    std::string action_log_filename = log_dir + "\\action_" + 
                                    std::to_string(action) + "_sheet_" + 
                                    std::to_string(state.sheet_id[action]) + ".log";
    
    std::ofstream action_log(action_log_filename, std::ios::app);
    
    if (action_log.is_open()) {
        auto now = std::chrono::system_clock::now();
        std::time_t now_time = std::chrono::system_clock::to_time_t(now);
        std::tm now_tm;
        localtime_s(&now_tm, &now_time);
        
        char timestamp[64];
        std::strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", &now_tm);
        
        action_log << "==== Action Log: " << timestamp << " ====" << std::endl;
        action_log << "Session ID: " << session_id << std::endl;
        action_log << "Action: " << action << " Sheet ID: " << state.sheet_id[action] << " Energy: " << state.sheet_energy[action] << std::endl;
    }
    
    if (state.sheet_energy[action] <= 0 && state.check_error_choice(action)) {
        std::cout << "sheet_energy is 0, can't collapse" << std::endl;
        done = 0; // 0表示 因为模型选择错误导致的异常，应赋予极大惩罚
        if (action_log.is_open()) {
            action_log << "Done status: " << done << std::endl;
            action_log << "========== End of Action Log ==========" << std::endl << std::endl;
            action_log.close();
        }
        return done;
    } 
    else if (state.sheet_energy[action] <= 0) {
        done = 1; // 1表示 mesh优化结束，重置开始下一个episode
        if (action_log.is_open()) {
            action_log << "Done status: " << done << std::endl;
            action_log << "========== End of Action Log ==========" << std::endl << std::endl;
            action_log.close();
        }
        return done;
    }
    else {
        done = 2; // 2表示优化未结束，正常选择下一个sheet
        sheet_op.collapse_one_sheet2(get_sheet_byId(&tmesh, state.sheet_id[action], sheet_op));
        std::cout << "collapse sheet: " << state.sheet_id[action] << std::endl;
        sheet_op.get_mesh_sheet_number();
        std::cout << "get_mesh_sheet_number" << std::endl;
        get_singularity_num_op.generate_singularity_number(&tmesh);
        std::cout << "generate_singularity_number" << std::endl;
        tmesh.write_Qhex("data/test.Qhex");
        std::cout << "write_Qhex" << std::endl;
        state.clear();
        std::cout << "clear state" << std::endl;
        calc_state(&tmesh, state, sheet_op);
        std::cout << "calc_state" << std::endl;
        // 不再调用fill_state，让state包含所有的sheet
        if (action_log.is_open()) {
            action_log << "Done status: " << done << std::endl;
            action_log << "========== End of Action Log ==========" << std::endl << std::endl;
            action_log.close();
        }
        return done;
    }
}

float calc_reward(int current_singularity_num, get_singularity_number<TMesh> get_singularity_num_op) {
    std::cout << "new singularity_num: " << get_singularity_num_op.singualarity_id 
              << " current_singularity_num: " << current_singularity_num << std::endl;
    if (get_singularity_num_op.singualarity_id < current_singularity_num) 
        return 10.0;
    else if (get_singularity_num_op.singualarity_id > current_singularity_num) 
        return -10.0;
    else
        return -1.0;
}

void print_state(const State& state) {
    std::cout << "state: " << std::endl;
    for (size_t i = 0; i < state.size(); i++) {
        std::cout << state.sheet_id[i] << " " << state.sheet_energy[i] << std::endl;
    }
}

void compare_diff(const State& state1, const State& state2) {
    for (size_t i = 0; i < state1.size(); i++) {
        if (state1.sheet_energy[i] != state2.sheet_energy[i]) {
            std::cout << "id " << state1.sheet_id[i] << " state1: " << state1.sheet_energy[i] 
                      << " state2: " << state2.sheet_energy[i] << std::endl;
        }
    }
}

py::list state_to_list(const State& state) {
    py::list result;
    for (size_t i = 0; i < state.size(); i++) {
        py::list row;
        row.append(state.sheet_id[i]);
        row.append(state.sheet_energy[i]);
        result.append(row);
    }
    return result;
}

#endif // STATE_FUNCTIONS_H
