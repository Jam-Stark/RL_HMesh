#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <algorithm>
#include <random>

#include "tools.h"
#include "meshQuality.h"
#include "TopoHexMeshLib_sheet_singularity/src/sheet_operation.h"
#include "singularity_number.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <filesystem> // 添加 filesystem 头文件

#define PYBIND11_DETAILED_ERROR_MESSAGES
#include <pybind11/embed.h>
#include <pybind11/stl.h>
namespace py = pybind11;

std::string path;//path of input file

#include "TopoHexMeshLib_sheet_singularity/src/state_functions.h"

// 检测当前构建模式
inline std::string GetBuildMode() {
#if defined(_DEBUG) && !defined(NDEBUG)
    return "Debug";
#else
    return "Release";
#endif
}

// 用法说明函数
void print_usage() {
    std::cout << "用法: RL_Mesh_demo <mesh_directory_path> [--load-model] [--model-path=path] [--memory-path=path]" << std::endl;
    std::cout << "参数:" << std::endl;
    std::cout << "  <mesh_directory_path> 包含 .Qhex 文件的文件夹路径" << std::endl;
    std::cout << "  --load-model          加载之前保存的模型和记忆(可选，默认不加载)" << std::endl;
    std::cout << "  --model-path=path     指定模型文件路径(可选，默认: python_modules/models/model.pt)" << std::endl;
    std::cout << "  --memory-path=path    指定记忆文件路径(可选，默认: python_modules/models/memory.pkl)" << std::endl;
}

int main(int argc, char* argv[])
{
    try {
        // 获取当前构建模式
        std::string build_mode = GetBuildMode();
        std::cout << "Current build mode: " << build_mode << std::endl;
        
        // Create base log directory
        std::string log_base_dir = "F://RL_HMesh//logs";
        system(("mkdir \"" + log_base_dir + "\" 2>NUL").c_str());
        
        // 根据构建模式创建子目录
        std::string mode_log_dir = log_base_dir + "/" + build_mode;
        system(("mkdir \"" + mode_log_dir + "\" 2>NUL").c_str());
        
        // Create timestamp
        auto now = std::chrono::system_clock::now();
        std::time_t now_time = std::chrono::system_clock::to_time_t(now);
        std::tm now_tm;
        localtime_s(&now_tm, &now_time);
        
        char timestamp[32];
        std::strftime(timestamp, sizeof(timestamp), "%Y%m%d_%H%M%S", &now_tm);
        std::string session_id = std::string(timestamp);
        
        // Create current session log subdirectory (在构建模式子目录下)
        std::string log_dir = mode_log_dir + "/" + session_id;
        system(("mkdir \"" + log_dir + "\" 2>NUL").c_str());
        
        // Set global session ID for other components
        _putenv_s("RL_HMESH_SESSION_ID", session_id.c_str());
        _putenv_s("RL_HMESH_BUILD_MODE", build_mode.c_str());
        
        // Main log file for standard output and errors
        std::string main_log_filename = log_dir + "/main.log";
        std::ofstream main_log(main_log_filename);
        if (!main_log.is_open()) {
            std::cerr << "Failed to open main log file: " << main_log_filename << std::endl;
            return -1;
        }
        

        
        // Create debug log file for detailed debug info
        std::string debug_log_filename = log_dir + "/debug.log";
        std::ofstream debug_log(debug_log_filename);
        if (!debug_log.is_open()) {
            std::cout << "Failed to open debug log file: " << debug_log_filename << std::endl;
            return -1;
        }


        debug_log << "Initializing Python interpreter..." << std::endl;
        py::scoped_interpreter guard{};
        debug_log << "Python interpreter initialized successfully" << std::endl;
        
        debug_log << "Importing sys module..." << std::endl;
        py::module_ sys = py::module_::import("sys");
        debug_log << "sys module imported successfully" << std::endl;
        
        debug_log << "Adding Python module path..." << std::endl;
        sys.attr("path").attr("append")("F:/RL_HMesh/python_modules");
        debug_log << "Python path added successfully: " << py::str(sys.attr("path")).cast<std::string>() << std::endl;
        
        debug_log << "Attempting to import agent module..." << std::endl;
        py::module_ agent_module;
        py::object agent;
        
        try {
            agent_module = py::module_::import("agent");
            debug_log << "Agent module imported successfully" << std::endl;
            
            debug_log << "Creating Agent instance..." << std::endl;
            agent = agent_module.attr("Agent")();
            debug_log << "Agent instance created successfully" << std::endl;
            
            // Check agent module functionality
            debug_log << "Checking agent module functionality..." << std::endl;
            if (py::hasattr(agent, "is_finish")) {
                debug_log << "Agent has is_finish method" << std::endl;
            } else {
                debug_log << "ERROR: Agent doesn't have is_finish method!" << std::endl;
            }
            
            if (py::hasattr(agent, "choose_action")) {
                debug_log << "Agent has choose_action method" << std::endl;
            } else {
                debug_log << "ERROR: Agent doesn't have choose_action method!" << std::endl;
            }
            
            if (py::hasattr(agent, "remember")) {
                debug_log << "Agent has remember method" << std::endl;
            } else {
                debug_log << "ERROR: Agent doesn't have remember method!" << std::endl;
            }
            
            if (py::hasattr(agent, "replay")) {
                debug_log << "Agent has replay method" << std::endl;
            } else {
                debug_log << "ERROR: Agent doesn't have replay method!" << std::endl;
            }
        } catch (py::error_already_set &e) {
            debug_log << "Failed to import agent module or create Agent instance: " << e.what() << std::endl;
            std::cerr << "Python error: " << e.what() << std::endl;
            PyErr_Print();
            return -1;
        }
        
        // 解析命令行参数
        if (argc < 2)
        {
            print_usage();
            return -1;
        }

        std::filesystem::path mesh_directory_path(argv[1]); // 获取文件夹路径
        if (!std::filesystem::is_directory(mesh_directory_path)) {
            std::cerr << "错误: 提供的路径不是一个有效的文件夹: " << argv[1] << std::endl;
            print_usage();
            return -1;
        }

        bool load_model = false;
        std::string model_path = "python_modules/models/model.pt";
        std::string memory_path = "python_modules/models/memory.pkl";
        
        // 解析其他命令行参数
        for(int i = 2; i < argc; i++)
        {
            std::string arg(argv[i]);
            if(arg == "--load-model") {
                load_model = true;
            }
            else if(arg.find("--model-path=") == 0) {
                model_path = arg.substr(13);  // 截取等号后面的部分
            }
            else if(arg.find("--memory-path=") == 0) {
                memory_path = arg.substr(14);  // 截取等号后面的部分
            }
            else {
                std::cout << "unkown args: " << arg << std::endl;
                print_usage();
                return -1;
            }
        }
        
        // 如果指定了加载模型选项，尝试加载模型和记忆
        if(load_model) {
            std::cout << "load exsited model&memory" << std::endl;
            bool load_success = agent.attr("load_checkpoint")(model_path, memory_path).cast<bool>();
            if(load_success) {
                std::cout << "Successfully loaded" << std::endl;
            } else {
                std::cout << "Warning: failed to load, trying to train " << std::endl;
            }
        }

        // 遍历文件夹中的文件
        for (const auto& entry : std::filesystem::directory_iterator(mesh_directory_path)) {
            if (entry.is_regular_file() && entry.path().extension() == ".Qhex") {
                std::string mesh_name = entry.path().string(); // 获取当前处理的文件完整路径
                std::cout << "\n=====================================================" << std::endl;
                std::cout << "Start to process mesh: " << mesh_name << std::endl;
                std::cout << "=====================================================" << std::endl;

                // 检查文件是否存在且可读 (虽然directory_iterator通常能保证，但多一层检查更安全)
                {
                    std::ifstream fin(mesh_name);
                    if (!fin.good()){
                        std::cerr << "错误: 文件无法找到或打开: " << mesh_name << std::endl;
                        continue; // 跳过这个文件，处理下一个
                    }
                }
                // 更新全局路径变量 (如果其他地方需要)
                // int pos = mesh_name.rfind("/", mesh_name.length());
                // if (pos != std::string::npos) {
                //     path = mesh_name.substr(0, pos + 1);
                // } else {
                //     path = "./"; // 或者其他默认值
                // }


        for(int episode = 0; episode < 100; episode++){
            std::cout << "Episode: " << episode << " for file: " << mesh_name << std::endl;

        /*---------- log setting ----------*/
            TMesh tmesh;
            sheet_operation<TMesh> sheet_op(&tmesh);
            get_singularity_number<TMesh> get_singularity_num_op(&tmesh);

            // 加载 .Qhex 文件
            std::cout << "Loading mesh file: " << mesh_name << std::endl;
            tmesh.load_Qhex(mesh_name.c_str());
            //std::cout << "Mesh loaded, checking mesh data..." << std::endl;

            //std::cout<< "boundary count by face: " << tmesh.count_boundary_byF() << std::endl;
            //std::cout<< "boundary count by edge: " << tmesh.count_boundary_byE() << std::endl;
            //std::cout<< "sharp edge count: " << tmesh.count_sharp_edges() << std::endl;
            //std::cout<< "corner count: " << tmesh.count_corner() << std::endl;

            // 使用正确的方法访问网格数据
            std::cout << "Vertex count: " << tmesh.vs.size() << std::endl;
            std::cout << "Edge count: " << tmesh.es.size() << std::endl;
            std::cout << "Face count: " << tmesh.fs.size() << std::endl;
            std::cout << "Cell count: " << tmesh.hs.size() << std::endl;
            std::cout << "Mesh file successfully loaded" << std::endl;

            // Redirect stdout and stderr to main log
            std::streambuf* cout_buf = std::cout.rdbuf();
            std::cout.rdbuf(main_log.rdbuf());
            std::streambuf* cerr_buf = std::cerr.rdbuf();
            std::cerr.rdbuf(main_log.rdbuf());
            std::cout << "Log started at " << timestamp << std::endl;
            std::cout << "Session ID: " << session_id << std::endl;
            std::cout << "Build Mode: " << build_mode << std::endl;
            std::cout << "================== Main Program Log for " << mesh_name << " ==================" << std::endl; // 添加文件名到日志

            int current_singularity_num = 0;
            float total_reward = 0.0;
            State state;
            //std::cout << "Initial state:" << std::endl;
            state.print();
            
            // First reset all edge sheet attributes
            //std::cout << "Resetting edge sheet attributes..." << std::endl;
            //sheet_op.reset_sheet_attributes();
            std::cout << "Calculating mesh sheet number..." << std::endl;
            sheet_op.get_mesh_sheet_number();
            std::cout << "Computing edge energy..." << std::endl;
            sheet_op.compute_edge_energy();
            
            std::cout << "Calculating state..." << std::endl;
            calc_state(&tmesh, state, sheet_op);
            std::cout << "State calculation complete, current state:" << std::endl;
            state.print();

            //std::cin.get(); // 等待用户输入，确保状态打印完成后再继续

            std::cout << "Generating singularity numbers..." << std::endl;
            get_singularity_num_op.generate_singularity_number(&tmesh);
            current_singularity_num = get_singularity_num_op.singualarity_id;

                std::cout << "Episode: " << episode << " for file: " << mesh_name << std::endl; // 明确指出是哪个文件的episode
                std::cout << "Current singularity number: " << current_singularity_num << std::endl;

                // 检查状态是否有效：如果所有边的能量都小于等于0，则跳过此 episode
                bool all_energy_non_positive = true; // Assume true initially
                if (state.size() > 0) { // 确保状态不为空
                    // 使用索引循环访问 sheet_energy
                    for (size_t i = 0; i < state.size(); ++i) {
                        // Access energy using the correct member vector
                        if (state.sheet_energy[i] > 0) { // Check if any energy is positive
                            all_energy_non_positive = false; // Found a positive energy, so not all are non-positive
                            break; // No need to check further
                        }
                    }
                } else {
                    // 如果状态本身为空，也视为无效状态 (或者说，没有正能量边)
                    all_energy_non_positive = true;
                }

                // If all energies are non-positive (<= 0), skip this episode
                if (all_energy_non_positive) {
                    std::cout << "INFO: All edge energies are non-positive or state is empty. No further optimization possible for this episode for file " << mesh_name << std::endl;
                    debug_log << "All edge energies are non-positive or state is empty after calculation for file " << mesh_name << ", aborting episode" << std::endl;
                    break; // Skip this file and try next
                }

                while(agent.attr("is_finish")(state_to_list(state)).cast<bool>()){
                    State state_snapshot = state;  // Take a snapshot of the state
                    int action = agent.attr("choose_action")(state_to_list(state_snapshot), episode).cast<int>();
                    std::cout << "Action: " << action << std::endl;
                    
                    debug_log << "State size before action: " << state.size() << std::endl;
                    int done = play_action(action, 1, state, tmesh, sheet_op, get_singularity_num_op);
                    debug_log << "State size after action: " << state.size() << std::endl;
                    std::cout << "Action result: " << done << std::endl;

                    if(done==0){
                        std::cout << "Critical error: sheet energy <= 0. Ending episode with penalty." << std::endl;
                        total_reward += -100000.0f; // Large negative reward
                        //log(episode, state_snapshot, action, state, total_reward); // Use snapshot
                        agent.attr("remember")(state_to_list(state_snapshot), action, state_to_list(state), total_reward);
                        break; // End current training episode immediately
                    }
                    else if(done==1){
                        std::cout << "Optimization finished. Ending episode." << std::endl;
                        total_reward += calc_reward(current_singularity_num, get_singularity_num_op, state_snapshot, action, state);
                        std::cout << "Total reward: " << total_reward << std::endl;
                        //log(episode, state_snapshot, action, state, total_reward);  
                        agent.attr("remember")(state_to_list(state_snapshot), action, state_to_list(state), total_reward);
                        break; // End current training episode immediately
                    }
                    else if(done==2){
                        std::cout << "Optimization not finished. Continuing." << std::endl;
                    }
                    total_reward += calc_reward(current_singularity_num, get_singularity_num_op, state_snapshot, action, state);
                    std::cout << "Total reward: " << total_reward << std::endl;
                    //log(episode, state_snapshot, action, state, total_reward);
                    std::cout << "Log completed" << std::endl;
                    agent.attr("remember")(state_to_list(state_snapshot), action, state_to_list(state), total_reward);
                    std::cout << "Memory update completed" << std::endl;

                }
                
                if (state.size() > 0) {
                    std::cout << "Training normally ended, saving model" << std::endl;
                    agent.attr("replay")();

                } else {
                    std::cout << "Training abnormally ended, state is empty" << std::endl;
                }

            //     // save to data\results
            //     std::string result_dir = "F://RL_HMesh//data//results//";
            //     system(("mkdir \"" + result_dir + "\" 2>NUL").c_str());

            // // 从mesh_name提取文件名部分 (不含扩展名)
            // std::string base_filename = entry.path().stem().string();

            // std::string result_file = result_dir + base_filename + "_processed_" + std::to_string(episode) + ".Qhex";
            // std::cout << "Saving collapsed mesh to: " << result_file << std::endl;
            // tmesh.write_Qhex(result_file.c_str());


                // Restore cout and cerr at the end of training for this file
                std::cout.rdbuf(cout_buf);
                std::cerr.rdbuf(cerr_buf);
            } // end episode loop for one file

            std::cout << "process down: " << mesh_name << std::endl;
                    
            // 确保目录存在
            system("mkdir \"python_modules\\models\" 2>NUL");
            
            // 使用save_checkpoint函数同时保存模型和记忆
            std::cout << "Saving model and replay memory..." << std::endl;
            agent.attr("save_checkpoint")(model_path, memory_path);
            std::cout << "Model and replay memory saved successfully" << std::endl;
            } // end if .Qhex file
        } // end directory iteration loop

        std::cout << "所有 .Qhex 文件处理完成。日志保存在 " << log_dir << std::endl;
        // Restore original cout/cerr buffers if they were redirected outside the loop (though in this code they are redirected inside)
        // std::cout.rdbuf(original_cout_buf); // Need to save original buffers at the start if needed
        // std::cerr.rdbuf(original_cerr_buf);
        return 0; // Ensure return 0 is outside the main loop
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
