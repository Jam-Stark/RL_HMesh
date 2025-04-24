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

int main(int argc, char* argv[])
{
    try {
        /*---------- log setting ----------*/
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
        
        // Redirect stdout and stderr to main log
        std::streambuf* cout_buf = std::cout.rdbuf();
        std::cout.rdbuf(main_log.rdbuf());
        std::streambuf* cerr_buf = std::cerr.rdbuf();
        std::cerr.rdbuf(main_log.rdbuf());
        
        // Create debug log file for detailed debug info
        std::string debug_log_filename = log_dir + "/debug.log";
        std::ofstream debug_log(debug_log_filename);
        if (!debug_log.is_open()) {
            std::cout << "Failed to open debug log file: " << debug_log_filename << std::endl;
            return -1;
        }
        
        std::cout << "Log started at " << timestamp << std::endl;
        std::cout << "Session ID: " << session_id << std::endl;
        std::cout << "Build Mode: " << build_mode << std::endl;
        std::cout << "================== Main Program Log ==================" << std::endl;
        
        /*---------- log setting ----------*/

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
        
        if (argc != 2)
        {
            std::cout << "Usage: RL_Mesh_demo <mesh_file_path>" << std::endl;
            return -1;
        }
        std::string mesh_name(argv[1]);
        {
            std::ifstream fin(mesh_name);
            if (!fin.good()){
                std::cerr << "File not found or cannot be opened: " << mesh_name << std::endl;
                return -1;
            }
        }
        int pos = mesh_name.rfind("/", mesh_name.length());
        path = mesh_name.substr(0, pos+1);

        for(int episode = 0; episode < 1; episode++){
            TMesh tmesh;
            sheet_operation<TMesh> sheet_op(&tmesh);
            get_singularity_number<TMesh> get_singularity_num_op(&tmesh);
            if (strutil::endsWith(mesh_name, ".Qhex"))
            {
                std::cout << "Loading mesh file: " << mesh_name << std::endl;
                tmesh.load_Qhex(mesh_name.c_str());
                //std::cout << "Mesh loaded, checking mesh data..." << std::endl;
                
                // 使用正确的方法访问网格数据
                std::cout << "Vertex count: " << tmesh.vs.size() << std::endl;
                std::cout << "Edge count: " << tmesh.es.size() << std::endl;
                std::cout << "Face count: " << tmesh.fs.size() << std::endl;
                std::cout << "Cell count: " << tmesh.hs.size() << std::endl;
            }
            else
            {
                std::cout << "Invalid input file format: " << mesh_name << std::endl;
                return -1;
            }
            std::cout << "Mesh file successfully loaded" << std::endl;

            int current_singularity_num = 0;
            float total_reward = 0.0;
            State state;
            //std::cout << "Initial state:" << std::endl;
            state.print();
            
            // First reset all edge sheet attributes
            std::cout << "Resetting edge sheet attributes..." << std::endl;
            sheet_op.reset_sheet_attributes();
            std::cout << "Calculating mesh sheet number..." << std::endl;
            sheet_op.get_mesh_sheet_number();
            std::cout << "Computing edge energy..." << std::endl;
            sheet_op.compute_edge_energy();
            
            std::cout << "Calculating state..." << std::endl;
            calc_state(&tmesh, state, sheet_op);
            std::cout << "State calculation complete, current state:" << std::endl;
            state.print();

            // std::cout << "Vertex count: " << tmesh.vs.size() << std::endl;
            // std::cout << "Edge count: " << tmesh.es.size() << std::endl;
            // std::cout << "Face count: " << tmesh.fs.size() << std::endl;
            // std::cout << "Cell count: " << tmesh.hs.size() << std::endl;

            std::cout << "Generating singularity numbers..." << std::endl;
            get_singularity_num_op.generate_singularity_number(&tmesh);
            current_singularity_num = get_singularity_num_op.singualarity_id;

            std::cout << "Episode: " << episode << std::endl;
            std::cout << "Current singularity number: " << current_singularity_num << std::endl;

            if (state.size() == 0) {
                std::cout << "ERROR: State is empty, cannot continue with reinforcement learning algorithm" << std::endl;
                debug_log << "State is empty after calculation, aborting episode" << std::endl;
                continue; // Skip this episode and try next one
            }
            
            while(agent.attr("is_finish")(state_to_list(state)).cast<bool>()){
                State state_snapshot = state;  // Take a snapshot of the state
                int action = agent.attr("choose_action")(state_to_list(state_snapshot), episode).cast<int>();
                std::cout << "Action: " << action << std::endl;
                
                debug_log << "State size before action: " << state.size() << std::endl;
                int done = play_action(action, 1, state, tmesh, sheet_op, get_singularity_num_op);
                debug_log << "State size after action: " << state.size() << std::endl;
                std::cout << "Action result: " << done << std::endl;

                // std::cout << "Vertex count: " << tmesh.vs.size() << std::endl;
                // std::cout << "Edge count: " << tmesh.es.size() << std::endl;
                // std::cout << "Face count: " << tmesh.fs.size() << std::endl;
                // std::cout << "Cell count: " << tmesh.hs.size() << std::endl;

                if(done==0){
                    std::cout << "Critical error: sheet energy <= 0. Ending episode with penalty." << std::endl;
                    total_reward += -1000000.0f; // Large negative reward
                    log(episode, state_snapshot, action, state, total_reward); // Use snapshot
                    agent.attr("remember")(state_to_list(state_snapshot), action, state_to_list(state), total_reward);
                    break; // End current training episode immediately
                }
                else if(done==1){
                    std::cout << "Optimization finished. Ending episode." << std::endl;
                    total_reward += calc_reward(current_singularity_num, get_singularity_num_op);
                    std::cout << "Total reward: " << total_reward << std::endl;
                    log(episode, state_snapshot, action, state, total_reward);  
                    agent.attr("remember")(state_to_list(state_snapshot), action, state_to_list(state), total_reward);
                    break; // End current training episode immediately
                }
                else if(done==2){
                    std::cout << "Optimization not finished. Continuing." << std::endl;
                }
                total_reward += calc_reward(current_singularity_num, get_singularity_num_op);
                std::cout << "Total reward: " << total_reward << std::endl;
                log(episode, state_snapshot, action, state, total_reward);
                std::cout << "Log completed" << std::endl;
                agent.attr("remember")(state_to_list(state_snapshot), action, state_to_list(state), total_reward);
                std::cout << "Memory update completed" << std::endl;

            }
            
            if (state.size() > 0) {
                std::cout << "Training normally ended, saving model" << std::endl;
                agent.attr("replay")();
                agent.attr("save_model")("python_modules/models/model.pt");
            } else {
                std::cout << "Training abnormally ended, state is empty" << std::endl;
            }

        // save to data\results
        std::string result_dir = "F://RL_HMesh//data//results//";
        system(("mkdir \"" + result_dir + "\" 2>NUL").c_str());
        
        // 从mesh_name提取文件名部分
        std::string base_filename = mesh_name.substr(mesh_name.find_last_of("/\\") + 1);
        base_filename = base_filename.substr(0, base_filename.find_last_of("."));
        
        std::string result_file = result_dir + base_filename + "_" + std::to_string(episode) + ".Qhex";
        std::cout << "Saving collapsed  mesh to: " << result_file << std::endl;
        tmesh.write_Qhex(result_file.c_str());
        }
        std::cout << "Training completed." << std::endl;

        // Restore cout and cerr at the end of training
        std::cout.rdbuf(cout_buf);
        std::cerr.rdbuf(cerr_buf);
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
