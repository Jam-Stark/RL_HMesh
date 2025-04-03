#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <algorithm>
#include <random>

#include "tools.h"
#include "meshQuality.h"
#include "sheet_operation.h"
#include "singularity_number.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#define PYBIND11_DETAILED_ERROR_MESSAGES
#include <pybind11/embed.h>
#include <pybind11/stl.h>
namespace py = pybind11;

std::string path;//path of input file

#include "state_functions.h"

int main(int argc, char* argv[])
{
    try {
        // 创建基础日志目录
        std::string log_base_dir = "F:/RL_HMesh/logs";
        system(("mkdir \"" + log_base_dir + "\" 2>NUL").c_str());
        
        // 创建时间戳
        auto now = std::chrono::system_clock::now();
        std::time_t now_time = std::chrono::system_clock::to_time_t(now);
        std::tm now_tm;
        localtime_s(&now_tm, &now_time);
        
        char timestamp[32];
        std::strftime(timestamp, sizeof(timestamp), "%Y%m%d_%H%M%S", &now_tm);
        std::string session_id = std::string(timestamp);
        
        // 创建当前会话的日志子目录
        std::string log_dir = log_base_dir + "/" + session_id;
        system(("mkdir \"" + log_dir + "\" 2>NUL").c_str());
        
        // 设置全局会话ID供其他组件使用
        // 使用环境变量设置会话ID
        _putenv_s("RL_HMESH_SESSION_ID", session_id.c_str());
        
        std::string main_log_filename = log_dir + "/main.log";
        std::ofstream main_log(main_log_filename);
        
        // 重定向cout到日志文件
        std::streambuf* cout_buf = std::cout.rdbuf();
        std::cout.rdbuf(main_log.rdbuf());
        
        std::cout << "Log started at " << timestamp << std::endl;
        std::cout << "Session ID: " << session_id << std::endl;
        std::cout << "================== Main Program Log ==================" << std::endl;
        
        py::scoped_interpreter guard{};
        py::module_ sys = py::module_::import("sys");
        sys.attr("path").attr("append")("F:/RL_HMesh/python_modules");
        py::module_ agent_file = py::module_::import("agent");
		py::object agent = agent_file.attr("Agent")();
        try {
            py::module_::import("agent");
        } catch (py::error_already_set &e) {
            std::cerr << "1.2 Failed to import agent module: " << e.what() << std::endl;
        }
        std::cout << "Python agent module found." << std::endl;

        if (argc != 2)
        {
            std::cout << "2 Usage: input.Qhex" << std::endl;
            return -1;
        }
        std::string mesh_name(argv[1]);
        {
            std::ifstream fin(mesh_name);
            if (!fin.good()){
                std::cerr << "2.1 File not found or cannot be opened: " << mesh_name << std::endl;
                return -1;
            }
        }
        int pos = mesh_name.rfind("/", mesh_name.length());
        path = mesh_name.substr(0, pos+1);
        std::cout << mesh_name << std::endl;

        for(int episode = 0; episode < 50; episode++){
            TMesh tmesh;
            sheet_operation<TMesh> sheet_op(&tmesh);
            get_singularity_number<TMesh> get_singularity_num_op(&tmesh);
            if (strutil::endsWith(mesh_name, ".Qhex"))
            {
                tmesh.load_Qhex(mesh_name.c_str());
                //std::cout << "3.2 finish read" << std::endl;
            }
            else
            {
                std::cout << "the wrong input file " << std::endl;
            }
            std::cout << "mesh file finish read" << std::endl;

            int current_singularity_num = 0;
            float total_reward = 0.0;
            State state;
            state.print();
            // 首先重置所有边的sheet属性
		    sheet_op.reset_sheet_attributes();
            sheet_op.get_mesh_sheet_number();
            sheet_op.compute_edge_energy();
            calc_state(&tmesh, state, sheet_op);

            std::cin.get();
            get_singularity_num_op.generate_singularity_number(&tmesh);
            current_singularity_num = get_singularity_num_op.singualarity_id;

            std::cout << "3.4 episode: " << episode << std::endl;
            std::cout << "current_singularity_num: " << current_singularity_num << std::endl;

            while(agent.attr("is_finish")(state_to_list(state)).cast<bool>()){
                State state_snapshot = state;  // 对 state 做快照
                int action = agent.attr("choose_action")(state_to_list(state_snapshot),episode).cast<int>();
                std::cout << "action :" << action << std::endl;
                int done = play_action(action, 1, state, tmesh, sheet_op, get_singularity_num_op);
                std::cout << "done: " << done << std::endl;
                if(done==0){
                    std::cout << "Critical error: sheet energy <= 0. Ending episode with penalty." << std::endl;
                    total_reward += -1000000.0f; // 极大负奖励
                    log(episode, state_snapshot, action, state, total_reward); // 使用快照
                    agent.attr("remember")(state_to_list(state_snapshot), action, state_to_list(state), total_reward);
                    break; // 立即结束当前回合训练
                }
                else if(done==1){
                    std::cout << "Optimization finished. Ending episode." << std::endl;
                    total_reward += calc_reward(current_singularity_num, get_singularity_num_op);
                    std::cout << "total_reward: " << total_reward << std::endl;
                    log(episode, state_snapshot, action, state, total_reward);
                    agent.attr("remember")(state_to_list(state_snapshot), action, state_to_list(state), total_reward);
                    break; // 立即结束当前回合训练
                }
                else if(done==2){
                    std::cout << "Optimization not finished. Continuing." << std::endl;
                }
                total_reward += calc_reward(current_singularity_num, get_singularity_num_op);
                std::cout << "total_reward: " << total_reward << std::endl;
                log(episode, state_snapshot, action, state, total_reward);
                std::cout << "log finished" << std::endl;
                agent.attr("remember")(state_to_list(state_snapshot), action, state_to_list(state), total_reward);
                std::cout << "remember finished" << std::endl;
                // 下次循环会重新快照 state
            }
            agent.attr("replay")();
            agent.attr("save_model")("python_modules/models/model.pt");
        }
        std::cout << "3.4 Training finished." << std::endl;
        
        // 训练结束后恢复cout
        std::cout.rdbuf(cout_buf);
        std::cout << "Training completed. Logs saved to " << log_dir << std::endl;
        
        return 0;
    }
    catch(py::error_already_set &e) {
        std::cerr << "Pybind11 error: " << e.what() << std::endl;
        PyErr_Print();
        return -1;
    }
    catch(const std::exception &ex) {
        std::cerr << "Exception: " << ex.what() << std::endl;
        return -1;
    }
    catch(...) {
        std::cerr << "Unknown exception occurred." << std::endl;
        return -1;
    }
    return 0;
}
