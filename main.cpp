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
#include <filesystem>
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
    std::cout << "  --inference 是否使用现有模型推理验证算法（默认进行训练）" << std::endl;
    std::cout << "  --GA 是否使用贪心算法计算一遍进行对比（默认不使用）" << std::endl;
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
        std::ofstream test_log(log_dir + "/test.log");

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
        
        agent_module = py::module_::import("agent");
        debug_log << "Agent module imported successfully" << std::endl;
        
        debug_log << "Creating Agent instance..." << std::endl;
        agent = agent_module.attr("Agent")(
            // --- 传递新的超参数 ---
            py::arg("max_steps_per_episode") = 10, // 示例值
            py::arg("reward_threshold_to_stop") = 100.0, // 示例值
            py::arg("no_improvement_steps_threshold") = 5, // 示例值
            py::arg("improvement_tolerance") = 1.0, // 示例值
            // --- 你可能还需要传递其他原有的构造函数参数 ---
            py::arg("eta") = 0.5, // 示例
            py::arg("gamma") = 0.99, // 示例
            py::arg("use_geometric_features") = true // 示例
            // ... 其他参数 ...
        );
        
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

        bool train_or_inference = false; // false for train, true for inference
        bool use_GA = false; // false for train, true for inference
        
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
            else if(arg.find("--inference"))
                train_or_inference = true;
            else if(arg.find("--GA"))
                use_GA = true;
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

        // 计算文件夹中的文件总数
        int total_files = 0;
        for (const auto& _ : std::filesystem::directory_iterator(mesh_directory_path)) {
            if (_.is_regular_file() && _.path().extension() == ".Qhex") {
                total_files++;
            }
        }
        
        std::cout << "Total .Qhex files found: " << total_files << std::endl;
        int total_episodes;
        int file_count = 0;
        // 遍历文件夹中的文件
        for (const auto& entry : std::filesystem::directory_iterator(mesh_directory_path)) {
            // 更新已处理文件数量
            if (entry.is_regular_file() && entry.path().extension() == ".Qhex") {
                file_count++;
                std::cout << "Processing file " << file_count << " of " << total_files << " (" 
                          << (file_count * 100 / total_files) << "%)" << std::endl;
            }

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

            if(train_or_inference==false){

                for(int episode = 0; episode < 50; episode++){
                    std::cout << "Episode: " << episode << " for file: " << mesh_name << std::endl;
                    agent.attr("reset_episode_stats")(); // <-- 在 episode 循环开始处调用
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

                    int original_hex_count = tmesh.hs.size();
                    int original_edge_count = tmesh.es.size();
                    int original_vertex_count = tmesh.vs.size();
                    int original_face_count = tmesh.fs.size();

                    // 使用正确的方法访问网格数据
                    std::cout << "Vertex count: " << original_vertex_count << std::endl;
                    std::cout << "Edge count: " << original_edge_count << std::endl;
                    std::cout << "Face count: " << original_face_count << std::endl;
                    std::cout << "Cell count: " <<  original_hex_count << std::endl;
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
                    //state.print();
                    
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
                    //state.print();

                    std::cout << "Generating singularity numbers..." << std::endl;
                    get_singularity_num_op.generate_singularity_number(&tmesh);
                    current_singularity_num = get_singularity_num_op.singualarity_id;
                    int original_singularity_num = current_singularity_num;

                    double original_jacobian = calculate_average_min_jacobian(&tmesh);

                    std::cout << "Original minimum Jacobian: " << original_jacobian << std::endl;

                   // std::cin.get(); // 等待用户输入，防止程序自动继续


                    std::cout << "Episode: " << episode << " for file: " << mesh_name << std::endl; // 明确指出是哪个文件的episode
                    std::cout << "Original singularity number: " <<original_singularity_num << std::endl;

                    // 检查状态是否有效：如果所有边的能量都小于等于0，则跳过此 episode
                    int energy_non_positive_count = 0; // Assume true initially
                    if (state.size() > 0) { // 确保状态不为空
                        // 使用索引循环访问 sheet_energy
                        for (size_t i = 0; i < state.size(); ++i) {
                            // Access energy using the correct member vector
                            if (state.sheet_energy[i] > 0) { // Check if any energy is positive
                                energy_non_positive_count ++; // Found a positive energy, so not all are non-positive 
                            }
                        }
                    }
                    // If all energies are non-positive (<= 0), skip this episode
                    if (energy_non_positive_count<=10) {
                        std::cout << "INFO: All edge energies are non-positive or state is empty. No further optimization possible for this episode for file " << mesh_name << std::endl;
                        debug_log << "All edge energies are non-positive or state is empty after calculation for file " << mesh_name << ", aborting episode" << std::endl;
                        break; // Skip this file and try next
                    }

                    float total_reward_for_logging = 0.0; // 可以保留一个用于日志的累计
                    int current_episode_steps_cpp = 0;
                    int prev_singularity_num_for_calc_reward = current_singularity_num;
                    //int original_singularity_num_at_episode_start = current_singularity_num;
                   while(true) {
                        State state_snapshot_St = state; // 当前状态 S_t
                        py::list py_St = state_to_list(state_snapshot_St);

                        py::tuple action_result = agent.attr("choose_action")(py_St, episode).cast<py::tuple>();
                        int action_At = action_result[0].cast<int>();
                        float done_logit_for_St = action_result[1].cast<float>();
                        // 使用 sigmoid 将 logit 转换为概率，然后与阈值比较更标准
                        bool network_thinks_St_is_done = (1.0f / (1.0f + std::exp(-done_logit_for_St))) > 0.5f; 
                        //bool network_thinks_St_is_done = (done_logit_for_St > 0.0f); // logit > 0 等价于 sigmoid(logit) > 0.5

                        std::cout << "Action: " << action_At << ", Done logit: " << done_logit_for_St 
                                << ", NetworkThinksDone: " << network_thinks_St_is_done << std::endl;

                        bool force_continue_override = false;
                        // **在这里加入您判断是否要覆盖网络“提前终止”的逻辑**
                        // 例如：如果是 训练 的first sevaral steps，并且初始状态并非真的不可操作，则强制继续
                        if (network_thinks_St_is_done == true && current_episode_steps_cpp <3 /* && !is_genuinely_terminal(state_snapshot_St) */) {
                            std::cout << "Network suggested DONE at step 0. Overriding to CONTINUE." << std::endl;
                            force_continue_override = true;
                        }

                        float reward_for_remember = 0.0f;
                        bool done_flag_for_remember = false;
                        bool is_error_for_remember = false;
                        py::list py_next_state_for_remember = py_St; // 默认情况下，如果原地终止，next_state 就是 current_state

                        if (network_thinks_St_is_done ==true && force_continue_override==false) {
                            // 情况1：网络认为 S_t 终止，C++不覆盖 (例如，不是第一步，或者S_t确实可能是一个合理的终止点)
                            // play_action 不执行。计算基于 S_t 的终局奖励。
                            get_singularity_num_op.generate_singularity_number(&tmesh); // 确保tmesh是S_t状态
                            int current_sing_at_St = get_singularity_num_op.singualarity_id;
                            double current_jac_at_St = calculate_average_min_jacobian(&tmesh);

                            // 调用 evaluate_overall_mesh_quality 时，请确保参数顺序正确
                            reward_for_remember = evaluate_overall_mesh_quality(
                                original_singularity_num, current_sing_at_St,
                                original_jacobian, current_jac_at_St);
                            
                            done_flag_for_remember = true; // 因为 训练 在此终止
                            // py_next_state_for_remember 已经是 py_St (或一个表示终止的特殊状态)
                            std::cout << "Terminating based on network prediction for S_t. Terminal Reward: " << reward_for_remember << std::endl;

                        } else {
                            // 情况2：网络认为 S_t 不终止，或者网络认为终止但被C++覆盖了(force_continue_override = true)
                            // 因此，实际执行动作 A_t
                            debug_log << "State size before action: " << state_snapshot_St.size() << std::endl; // 用 state_snapshot_St
                            // play_action会修改 state，使其从 S_t 变为 S_{t+1}
                            int done_code_from_env = play_action(action_At, 1, state, tmesh, sheet_op, get_singularity_num_op, original_hex_count);
                            debug_log << "State size after action: " << state.size() << std::endl; // state 现在是 S_{t+1}
                            std::cout << "play_action result code: " << done_code_from_env << std::endl;

                            py_next_state_for_remember = state_to_list(state); // S_{t+1}

                            // 计算单步奖励 R_step (S_t, A_t) -> S_{t+1}
                            // 注意 calc_reward 的第一个参数是执行动作前的奇异线数量
                            reward_for_remember = calc_reward(prev_singularity_num_for_calc_reward, get_singularity_num_op, state_snapshot_St, action_At, state);
                            std::cout << "Step reward from calc_reward: " << reward_for_remember << std::endl;

                            if (done_code_from_env == 0) { // C++ 环境错误
                                std::cout << "Terminating due to environment error (play_action returned 0)." << std::endl;
                                reward_for_remember = -100000.0f; // 重写奖励为巨大惩罚
                                done_flag_for_remember = true;
                                is_error_for_remember = true;
                            } else if (done_code_from_env == 1) { // C++ 环境判定成功终止 (例如，网格元素数量达到阈值)
                                std::cout << "Terminating due to C++ environment success (play_action returned 1)." << std::endl;
                                // 此时，S_{t+1} 是终态。我们可以用 evaluate_overall_mesh_quality 计算其终局奖励
                                get_singularity_num_op.generate_singularity_number(&tmesh); // 确保tmesh是S_{t+1}状态
                                int current_sing_at_St_plus_1 = get_singularity_num_op.singualarity_id;
                                double current_jac_at_St_plus_1 = calculate_average_min_jacobian(&tmesh);
                                
                                reward_for_remember = evaluate_overall_mesh_quality(
                                    original_singularity_num, current_sing_at_St_plus_1,
                                    original_jacobian, current_jac_at_St_plus_1);
                                std::cout << "Terminal reward for S_{t+1} (env. success): " << reward_for_remember << std::endl;
                                done_flag_for_remember = true;
                            } else if(done_code_from_env==2){ // done_code_from_env == 2 (正常进行)
                                // 如果是被 force_continue_override 的情况，即使网络之前说S_t要done，
                                // 但因为我们强制继续了，所以实际的 done_flag_for_remember 应该是 False。
                                done_flag_for_remember = false; 
                                // 这里还可以加入其他C++侧的终止条件判断，例如最大步数
                                // if (current_episode_steps_cpp >= MAX_STEPS_PER_EPISODE) {
                                //    done_flag_for_remember = true;
                                //    // 如果因为最大步数终止，也可以考虑调用 evaluate_overall_mesh_quality 计算最终状态 S_{t+1} 的质量作为奖励
                                // }
                                std::cout << "Continuing to next step (play_action returned 2)." << std::endl;
                            } else if(done_code_from_env == 3) { //判断正常实施action后的next_state是否应该为终止state,如果是则 done_flag_for_remember = true,用evaluate_overall_mesh_quality(S_{t+1})
                                std::cout << "Terminating based on action returned 3 for S_{t+1} is end state." << std::endl;
                                done_flag_for_remember = true;

                                // 关键: 重新计算奖励为基于 S_{t+1} 的终局奖励
                                get_singularity_num_op.generate_singularity_number(&tmesh); // 确保是 S_{t+1} 的奇异线数量
                                int current_sing_at_St_plus_1 = get_singularity_num_op.singualarity_id;
                                double current_jac_at_St_plus_1 = calculate_average_min_jacobian(&tmesh); // S_{t+1} 的雅可比值

                                reward_for_remember = evaluate_overall_mesh_quality(
                                    original_singularity_num, // 初始的奇异线数量
                                    current_sing_at_St_plus_1,  // S_{t+1} 的奇异线数量
                                    original_jacobian,         // 初始的雅可比值
                                    current_jac_at_St_plus_1   // S_{t+1} 的雅可比值
                                );

                                reward_for_remember += 100.0f; // 网络确实抵达了终止状态再结束的额外奖励
                                if(network_thinks_St_is_done == true)
                                    reward_for_remember += 100.0f; // 如果网络也认为是终止状态，给额外奖励
                                std::cout << "Terminal reward for S_{t+1} (agent.is_finish() triggered): " << reward_for_remember << std::endl;
                            }
                        }

                        // 记录经验
                        agent.attr("remember")(py_St, action_At, py_next_state_for_remember, reward_for_remember, done_flag_for_remember, is_error_for_remember);
                        total_reward_for_logging += reward_for_remember;

                        // 学习
                        agent.attr("replay")();

                        // 更新用于下一次 calc_reward 的奇异线计数 (应该是 S_{t+1} 的计数)
                        // 只有在 训练 未实际终止时才需要更新，为下一次迭代做准备
                        if (done_flag_for_remember == false) { // 如果 训练 还没有因为这一步而终止
                            std::cout<<"before get_singularity_num_op.generate_singularity_number(&tmesh): "<< prev_singularity_num_for_calc_reward<<std::endl;
                            get_singularity_num_op.generate_singularity_number(&tmesh); // 获取 S_{t+1} 的奇异线数
                            prev_singularity_num_for_calc_reward = get_singularity_num_op.singualarity_id;
                            current_episode_steps_cpp++;
                            std::cout << "Current episode steps (C++): " << current_episode_steps_cpp << std::endl;
                            std::cout << "Current singularity number: " << prev_singularity_num_for_calc_reward << std::endl;
                        }

                        if (done_flag_for_remember == true) { // 如果任何原因导致了这一步是终止步
                            std::cout << "EPISODE: " << episode << " terminated. Total logged reward: " << total_reward_for_logging 
                                    << ". Final reward for remember: " << reward_for_remember 
                                    << ". Done flag for remember: " << done_flag_for_remember << std::endl;
                            break; // 退出 while(true) 循环
                        }
                    } // 结束 while(true)
                
                    if (state.size() > 0) {
                        std::cout << "Training normally ended, saving model" << std::endl;
                        agent.attr("replay")();

                    } else {
                        std::cout << "Training abnormally ended, state is empty" << std::endl;
                    }

                    // Restore cout and cerr at the end of training for this file
                    std::cout.rdbuf(cout_buf);
                    std::cerr.rdbuf(cerr_buf);
                    } // end episode loop for one file
                }
                else if(train_or_inference==true){ //使用模型去实际运行，测试训练结果

                        // Redirect stdout and stderr to main log
                        std::streambuf* cout_buf = std::cout.rdbuf();
                        std::cout.rdbuf(test_log.rdbuf());
                        std::streambuf* cerr_buf = std::cerr.rdbuf();
                        std::cerr.rdbuf(test_log.rdbuf());
                        std::cout << "Log started at " << timestamp << std::endl;
                        std::cout << "================== Test Program Log for " << mesh_name << " ==================" << std::endl; // 添加文件名到日志

                        TMesh tmesh;
                        sheet_operation<TMesh> sheet_op(&tmesh);
                        get_singularity_number<TMesh> get_singularity_num_op(&tmesh);

                        // 加载 .Qhex 文件
                        std::cout << "Loading mesh file: " << mesh_name << std::endl;
                        tmesh.load_Qhex(mesh_name.c_str());

                        int original_hex_count = tmesh.hs.size();
                        int original_edge_count = tmesh.es.size();
                        int original_vertex_count = tmesh.vs.size();
                        int original_face_count = tmesh.fs.size();

                        // 使用正确的方法访问网格数据
                        std::cout << "Vertex count: " << original_vertex_count << std::endl;
                        std::cout << "Edge count: " << original_edge_count << std::endl;
                        std::cout << "Face count: " << original_face_count << std::endl;
                        std::cout << "Cell count: " <<  original_hex_count << std::endl;
                        std::cout << "Mesh file successfully loaded" << std::endl;

                        int current_singularity_num = 0;
                        State state;
                        
                        std::cout << "Calculating mesh sheet number..." << std::endl;
                        sheet_op.get_mesh_sheet_number();
                        std::cout << "Computing edge energy..." << std::endl;
                        sheet_op.compute_edge_energy();
                        
                        std::cout << "Calculating state..." << std::endl;
                        calc_state(&tmesh, state, sheet_op);
                        std::cout << "State calculation complete, current state:" << std::endl;
                        state.print();

                        std::cout << "Generating singularity numbers..." << std::endl;
                        get_singularity_num_op.generate_singularity_number(&tmesh);
                        current_singularity_num = get_singularity_num_op.singualarity_id;
                        int original_singularity_num = current_singularity_num;
                        std::cout << "Initial singularity number: " << current_singularity_num << std::endl;

                        double original_jacobian = calculate_average_min_jacobian(&tmesh);

                        int energy_non_positive_count = 0; 
                        if (state.size() > 0) { 
                            for (size_t i = 0; i < state.sheet_energy.size(); ++i) { 
                                if (state.sheet_energy[i] > 0) { 
                                    energy_non_positive_count++; 
                                }
                            }
                        }
                        
                        if (energy_non_positive_count <= 0 && state.size() > 0) { 
                            std::cout << "INFO: All edge energies are non-positive or state is empty. No further optimization possible for this file " << mesh_name << std::endl;
                            debug_log << "All edge energies are non-positive or state is empty after initial calculation for file " << mesh_name << ", skipping inference loop." << std::endl;
                        } else if (state.size() == 0) {
                             std::cout << "INFO: Initial state is empty. No optimization possible for this file " << mesh_name << std::endl;
                             debug_log << "Initial state is empty for file " << mesh_name << ", skipping inference loop." << std::endl;
                        }
                        else
                        {
                            int max_inference_steps = 200; 
                            int inference_steps_count = 0;

                            std::cout << "Starting inference process for mesh: " << mesh_name << std::endl;

                            while (true) { 
                                if (inference_steps_count >= max_inference_steps) {
                                    std::cout << "Inference: Reached maximum inference steps (" << max_inference_steps << "). Stopping." << std::endl;
                                    break;
                                }

                                if (state.size() == 0) {
                                    std::cout << "Inference: State is empty. Stopping." << std::endl;
                                    break;
                                }

                                bool positive_energy_exists = false;
                                for (double energy : state.sheet_energy) {
                                    if (energy > 0) {
                                        positive_energy_exists = true;
                                        break;
                                    }
                                }
                                if (!positive_energy_exists) {
                                    std::cout << "Inference: No sheets with positive energy remaining. Stopping." << std::endl;
                                    break;
                                }

                                std::cout << "Inference Step " << inference_steps_count + 1 << "/" << max_inference_steps << std::endl;

                                py::list current_state_py = state_to_list(state);
                                if (current_state_py.empty()) {
                                    std::cout << "Inference: state_to_list returned empty list. Stopping." << std::endl;
                                    break;
                                }
                                
                                // 使用choose_action代替get_action，获取动作和模型对终止的预测
                                py::tuple action_result = agent.attr("choose_action")(current_state_py, 0).cast<py::tuple>();
                                int action_idx = action_result[0].cast<int>();
                                float done_logit = action_result[1].cast<float>();
                                
                                // 打印终止预测结果
                                std::cout << "Action: " << action_idx << ", Done logit: " << done_logit << std::endl;
                                
                                // 检查模型是否预测应该停止
                                if (done_logit > 0.5) {
                                    std::cout << "Inference: Model predicted episode should end (done_logit = " 
                                            << done_logit << "). Stopping." << std::endl;
                                    debug_log << "Inference: Model predicted episode should end (done_logit = " 
                                            << done_logit << ") for file " << mesh_name << ". Stopping." << std::endl;
                                    break;
                                }

                                if (action_idx < 0 || static_cast<size_t>(action_idx) >= state.size()) {
                                    std::cerr << "Inference Error: Agent returned invalid action index: " << action_idx << " for state size " << state.size() << std::endl;
                                    debug_log << "Inference Error: Agent returned invalid action index: " << action_idx << " for state size " << state.size() << std::endl;
                                    break; 
                                }
                                
                                std::cout << "  Chosen action (index in state): " << action_idx
                                        << ", Sheet ID: " << state.sheet_id[action_idx]
                                        << ", Energy: " << state.sheet_energy[action_idx]
                                        << ""<< std::endl;
                                debug_log << "Inference Step " << inference_steps_count + 1 << ": Action index " << action_idx << ", Sheet ID " << state.sheet_id[action_idx] << " for file " << mesh_name << std::endl;

                                int play_status = play_action(action_idx, 2, state, tmesh, sheet_op, get_singularity_num_op, original_hex_count);

                                std::cout << "  Action play_status: " << play_status << std::endl;
                                debug_log << "  Play_action status: " << play_status << " for file " << mesh_name << std::endl;

                                if (play_status == 0) {
                                    std::cout << "Inference: Action resulted in an error (status 0 from play_action). Stopping." << std::endl;
                                    debug_log << "Inference: play_action returned 0 (error) for file " << mesh_name << ". Stopping." << std::endl;
                                    break;
                                } else if (play_status == 1) {
                                    std::cout << "Inference: Mesh collapsed to a point where it's too small (status 1 from play_action). Stopping." << std::endl;
                                    debug_log << "Inference: play_action returned 1 (mesh too small) for file " << mesh_name << ". Stopping." << std::endl;
                                    break;
                                }

                                get_singularity_num_op.generate_singularity_number(&tmesh);
                                current_singularity_num = get_singularity_num_op.singualarity_id;
                                std::cout << "  Singularity number after action: " << current_singularity_num << std::endl;
                                debug_log << "  Singularity number after action: " << current_singularity_num << " for file " << mesh_name << std::endl;

                                inference_steps_count++;
                            }
                             if (inference_steps_count >= max_inference_steps && max_inference_steps > 0) { 
                                std::cout << "Inference: Reached maximum inference steps (" << inference_steps_count << "/" << max_inference_steps << ") for file " << mesh_name << "." << std::endl;
                                debug_log << "Inference: Reached maximum inference steps for file " << mesh_name << "." << std::endl;
                            }
                        }
                        // 打印算法优化后的网格各项指标对比
                        std::cout << "Final mesh statistics after inference:" << std::endl;
                        std::cout << "Original hex count: " << original_hex_count << std::endl;
                        std::cout << "Original face count: " << original_face_count << std::endl;
                        std::cout << "Original Vertex count: " << original_vertex_count << std::endl;
                        std::cout << "Original Edge count: " << original_edge_count << std::endl;
                        std::cout << "Original singularity number: " << original_singularity_num << std::endl;
                        std::cout<< "Original jacobian metric: " << original_jacobian << std::endl;

                        std::cout << "Final Vertex count: " << tmesh.vs.size() << std::endl;
                        std::cout << "Final Edge count: " << tmesh.es.size() << std::endl;
                        std::cout << "Face count: " << tmesh.fs.size() << std::endl;
                        std::cout << "Final Cell count: " << tmesh.hs.size() << std::endl;
                        std::cout << "Final Current singularity number: " << current_singularity_num << std::endl;
                        //jacobian quality
                        std::cout<< "Final jacobian metric: " <<calculate_average_min_jacobian(&tmesh) << std::endl;
                        
                        //save to data\results
                        std::string result_dir = "F://RL_HMesh//data//results//";
                        system(("mkdir \"" + result_dir + "\" 2>NUL").c_str());
                        std::string base_filename = entry.path().stem().string();
                        std::string result_file = result_dir + base_filename + "_processed.Qhex";
                        std::cout << "Saving collapsed mesh to: " << result_file << std::endl;
                        tmesh.write_Qhex(result_file.c_str());

                        if(use_GA){
                            std::cout << "\n--- Starting Greedy Algorithm (GA) for comparison ---" << std::endl;
                            // 1. 为 GA 创建独立的网格和操作对象
                            TMesh ga_tmesh;
                            sheet_operation<TMesh> ga_sheet_op(&ga_tmesh);
                            get_singularity_number<TMesh> ga_get_singularity_num_op(&ga_tmesh);

                            std::cout << "GA: Loading original mesh file: " << mesh_name << std::endl;
                            ga_tmesh.load_Qhex(mesh_name.c_str()); // 加载网格副本
                            int ga_original_hex_count = ga_tmesh.hs.size();
                            std::cout << "GA: Original hex count: " << ga_original_hex_count << std::endl;

                            // 2. 计算 GA 的初始状态
                            State ga_state;
                            std::cout << "GA: Calculating initial mesh sheet number..." << std::endl;
                            ga_sheet_op.get_mesh_sheet_number();
                            std::cout << "GA: Computing initial edge energy..." << std::endl;
                            ga_sheet_op.compute_edge_energy();
                            std::cout << "GA: Calculating initial state..." << std::endl;
                            calc_state(&ga_tmesh, ga_state, ga_sheet_op);
                            std::cout << "GA: Initial state calculation complete. State size: " << ga_state.size() << std::endl;
                            // ga_state.print(); // 可选: 打印 GA 初始状态

                            ga_get_singularity_num_op.generate_singularity_number(&ga_tmesh);
                            int ga_current_singularity_num = ga_get_singularity_num_op.singualarity_id;
                            std::cout << "GA: Initial singularity number: " << ga_current_singularity_num << std::endl;

                            int ga_steps_count = 0;
                            int max_ga_steps = 300; // GA 的安全步数上限

                            // 3. GA 循环
                            while(ga_steps_count < max_ga_steps) {
                                if (ga_state.size() == 0) {
                                    std::cout << "GA: State is empty. Stopping." << std::endl;
                                    break;
                                }

                                int action_idx_ga = -1;
                                double max_energy_ga = 0.0; // 只考虑能量为正的 sheet

                                for (size_t i = 0; i < ga_state.sheet_energy.size(); ++i) {
                                    if (ga_state.sheet_energy[i] > max_energy_ga) {
                                        max_energy_ga = ga_state.sheet_energy[i];
                                        action_idx_ga = static_cast<int>(i);
                                    }
                                }

                                if (action_idx_ga == -1 || max_energy_ga <= 1e-9) { // 没有能量为正的 sheet
                                    std::cout << "GA: No sheets with positive energy remaining or no valid action. Stopping." << std::endl;
                                    break;
                                }

                                std::cout << "GA Step " << ga_steps_count + 1 << ": Choosing action (index in state): " << action_idx_ga
                                          << ", Sheet ID: " << ga_state.sheet_id[action_idx_ga]
                                          << ", Energy: " << ga_state.sheet_energy[action_idx_ga] << std::endl;
                                debug_log << "GA Step " << ga_steps_count + 1 << ": Action index " << action_idx_ga << ", Sheet ID " << ga_state.sheet_id[action_idx_ga] << " for file " << mesh_name << std::endl;
                                
                                // play_action 会就地修改 ga_state
                                int ga_play_status = play_action(action_idx_ga, 2, ga_state, ga_tmesh, ga_sheet_op, ga_get_singularity_num_op, ga_original_hex_count);

                                std::cout << "GA: Action play_status: " << ga_play_status << std::endl;
                                debug_log << "GA: Play_action status: " << ga_play_status << " for file " << mesh_name << std::endl;

                                if (ga_play_status == 0) { // play_action 发生错误
                                    std::cout << "GA: Action resulted in an error (status 0). Stopping." << std::endl;
                                    debug_log << "GA: play_action returned 0 (error) for file " << mesh_name << ". Stopping." << std::endl;
                                    break;
                                } else if (ga_play_status == 1) { // 网格过小
                                    std::cout << "GA: Mesh collapsed significantly (status 1). Stopping." << std::endl;
                                    debug_log << "GA: play_action returned 1 (mesh too small) for file " << mesh_name << ". Stopping." << std::endl;
                                    break;
                                }
                                // 若 ga_play_status == 2, 操作成功, ga_state 已更新

                                ga_get_singularity_num_op.generate_singularity_number(&ga_tmesh);
                                ga_current_singularity_num = ga_get_singularity_num_op.singualarity_id;
                                std::cout << "GA: Singularity number after action: " << ga_current_singularity_num << std::endl;
                                debug_log << "GA: Singularity number after action: " << ga_current_singularity_num << " for file " << mesh_name << std::endl;

                                ga_steps_count++;
                            }
                            if (ga_steps_count >= max_ga_steps) {
                                std::cout << "GA: Reached maximum GA steps (" << max_ga_steps << ")." << std::endl;
                                debug_log << "GA: Reached maximum GA steps for file " << mesh_name << "." << std::endl;
                            }

                            // 4. 保存 GA 处理后的网格
                            std::string ga_result_file = result_dir + base_filename + "_GA_processed.Qhex";
                            std::cout << "GA: Saving collapsed mesh to: " << ga_result_file << std::endl;
                            ga_tmesh.write_Qhex(ga_result_file.c_str());
                            std::cout << "--- Greedy Algorithm (GA) processing finished ---" << std::endl;
                        }
                }

            std::cout << "process completed: " << mesh_name << std::endl;
                    
            // 确保目录存在
            system("mkdir \"python_modules\\models\" 2>NUL");
            
            // 使用save_checkpoint函数同时保存模型和记忆
            std::cout << "Saving model and replay memory..." << std::endl;
            agent.attr("save_checkpoint")(model_path, memory_path);
            std::cout << "Model and replay memory saved successfully" << std::endl;
            } // end if .Qhex file
        } // end directory iteration loop

        std::cout << "All .Qhex processed。logs saved at " << log_dir << std::endl;
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
