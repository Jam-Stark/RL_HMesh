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
        
        /* 新增state 参数 Direct Feature/Boundary Involvement (Highest Importance Maybe)
        1 sheet_on_boundary_ratio:
        Calculation: Percentage of edges in the sheet that are boundary edges (TE::m_boundary == true) or percentage of unique vertices in the sheet that are boundary vertices (TV::m_boundary == true).
        Interpretation: A high value means collapsing this sheet will directly alter the boundary shape. This is a very strong signal against collapse unless specifically desired.

        2 sheet_on_feature_ratio:
        Calculation: Percentage of edges in the sheet that are marked as sharp/feature edges (TE::m_sharp == true or potentially other project-specific flags like TE::m_feature_curve).
        Interpretation: A high value indicates the sheet is the feature. Collapsing it will likely destroy or severely damage that feature. This is another strong signal against collapse.

        3 endpoints_are_corners_or_boundary:
        Calculation: Check the two extreme vertices of the sheet path. Are they marked as corners (TV::m_corner == true, TV::m_feature_vertex == true) or are they boundary vertices (TV::m_boundary == true)? You could use boolean flags for each end or a combined count/flag.
        Interpretation: Collapsing a sheet ending at a geometric corner or on the boundary is highly likely to distort the shape at that critical location.

        4 adjacent_feature_edges:
        Calculation: For each edge in the sheet, count how many of its neighboring edges (within adjacent faces/hexes) are marked as sharp/feature edges. Calculate an average or maximum count for the sheet.
        Interpretation: Indicates if the sheet runs parallel and close to a feature. Collapsing it might cause the adjacent feature to distort or snap undesirably to the new vertices.*/

        std::vector<double> sheet_on_boundary_ratio;
        std::vector<double> sheet_on_feature_ratio;
        std::vector<bool> sheet_endpoints_are_corners_or_boundary;
        std::vector<int> sheet_adjacent_feature_edges_count;

        /* Geometric Configuration (Proxy for Shape/Curvature):
        1 Sheet Curvature Metric:
        Calculation:
        Calculate the angle between consecutive edge segments along the sheet path. Find the average or maximum deviation from 180 degrees (straight).
        Alternatively, calculate the ratio: (Euclidean distance between sheet endpoints) / (Total geometric length of sheet edges). A value close to 1 means straight, smaller values mean more curved.
        Interpretation: Highly curved sheets often trace curved features of the original geometry. Collapsing them replaces the curve with potentially fewer, straighter segments, thus altering the shape.

        2 Average Normal Variation Along Sheet:
        Calculation: For each edge in the sheet, find the angle between the m_normal vectors of its two endpoint vertices (TV). Average these angles over the sheet.
        Interpretation: High average angle variation indicates the sheet lies on a region of high surface curvature. Collapsing the sheet might flatten this region.

        3 Dihedral Angle Deviation (m_total_angle analysis):
        Calculation: Analyze the TE::m_total_angle for edges in the sheet (pre-calculated or calculate on the fly).
        For internal edges, how much does the total_angle deviate from 360 degrees?
        For boundary edges, how much does the total_angle deviate from 180 degrees (or the specific angle if it's a sharp feature)? Calculate average or maximum deviation for the sheet.
        Interpretation: Large deviations often correspond to sharp features (intended) or areas of high curvature (also part of the shape). Collapsing edges with extreme dihedral angles is risky for shape preservation.
        The ideal_degree calculation in sheet_operation already uses this implicitly, but making the raw deviation explicit might help the network.*/

        std::vector<double> sheet_curvature_metric;

        // 存储法线变化的向量
        std::vector<double> sheet_normal_variation;
        
        // 存储二面角偏差的向量
        std::vector<double> sheet_dihedral_angle_deviation;

        /* Neighborhood Quality (Indicates Geometric Stress):
        Min/Average Scaled Jacobian of Adjacent Hexes:
        Calculation: Find all hexahedra (TH*) adjacent to the edges (TE*) of the sheet. Calculate the minimum and average Scaled Jacobian (CMeshQuality) for these hexes.
        Interpretation: If the surrounding elements are already of low quality or highly distorted,
        a collapse (which is a significant topological and geometric change) is more likely to lead to further degradation or invalid elements, potentially flattening or warping the local shape.*/

        // 存储缩放雅可比值的向量（最小值和平均值）
        std::vector<double> sheet_min_scaled_jacobian;
        std::vector<double> sheet_avg_scaled_jacobian;

        int sheet_num;
        // 默认构造函数
        State() {
            sheet_num = 0;
        }
		
		// 清空 state 数据
		void clear() {
			sheet_id.clear();
			sheet_energy.clear();
            sheet_on_boundary_ratio.clear();
            sheet_on_feature_ratio.clear();
            sheet_endpoints_are_corners_or_boundary.clear();
            sheet_adjacent_feature_edges_count.clear();
            sheet_curvature_metric.clear();
            sheet_normal_variation.clear();
            sheet_dihedral_angle_deviation.clear();
            sheet_min_scaled_jacobian.clear();
            sheet_avg_scaled_jacobian.clear();
            sheet_num = 0;
		}
		// 返回状态大小
		size_t size() const {
			return sheet_id.size();
		}
		// 添加一条记录
		void add(int id, double energy,double boundary_ratio,
            double feature_ratio,
            bool endpoints_special,
            int adjacent_features,
            double curvature_metric, 
            double normal_variation,
            double dihedral_deviation,
            double min_jacobian,
            double avg_jacobian) {
			sheet_id.push_back(id);
			sheet_energy.push_back(energy);
            sheet_on_boundary_ratio.push_back(boundary_ratio);
            sheet_on_feature_ratio.push_back(feature_ratio);
            sheet_endpoints_are_corners_or_boundary.push_back(endpoints_special);
            sheet_adjacent_feature_edges_count.push_back(adjacent_features);
            sheet_curvature_metric.push_back(curvature_metric);
            sheet_normal_variation.push_back(normal_variation);
            sheet_dihedral_angle_deviation.push_back(dihedral_deviation);
            sheet_min_scaled_jacobian.push_back(min_jacobian);
            sheet_avg_scaled_jacobian.push_back(avg_jacobian);
            sheet_num++;
		}
		
        // 重载赋值操作符
		State& operator=(const State& other) {
			if (this != &other) {
				sheet_id = other.sheet_id;
				sheet_energy = other.sheet_energy;
                sheet_on_boundary_ratio = other.sheet_on_boundary_ratio;
                sheet_on_feature_ratio = other.sheet_on_feature_ratio;
                sheet_endpoints_are_corners_or_boundary = other.sheet_endpoints_are_corners_or_boundary;
                sheet_adjacent_feature_edges_count = other.sheet_adjacent_feature_edges_count;
                sheet_curvature_metric = other.sheet_curvature_metric;
                sheet_normal_variation = other.sheet_normal_variation;
                sheet_dihedral_angle_deviation = other.sheet_dihedral_angle_deviation;
                sheet_min_scaled_jacobian = other.sheet_min_scaled_jacobian;
                sheet_avg_scaled_jacobian = other.sheet_avg_scaled_jacobian;
                sheet_num = other.sheet_num;
			}
			return *this;
		}
		
		//根据sheet_id查找sheet_energy等需要的对应位置的值
		double get_energy(int id) {
			for (size_t i = 0; i < sheet_id.size(); i++) {
				if (sheet_id[i] == id) {
					return sheet_energy[i];
				}
			}
			return 0.0;
		}

        // 根据sheet_id查找几何特征值
        double get_boundary_ratio(int id) {
            for (size_t i = 0; i < sheet_id.size(); i++) {
                if (sheet_id[i] == id) {
                    return sheet_on_boundary_ratio[i];
                }
            }
            return 0.0;
        }

        double get_feature_ratio(int id) {
            for (size_t i = 0; i < sheet_id.size(); i++) {
                if (sheet_id[i] == id) {
                    return sheet_on_feature_ratio[i];
                }
            }
            return 0.0;
        }

        bool get_endpoints_special(int id) {
            for (size_t i = 0; i < sheet_id.size(); i++) {
                if (sheet_id[i] == id) {
                    return sheet_endpoints_are_corners_or_boundary[i];
                }
            }
            return false;
        }

        int get_adjacent_features(int id) {
            for (size_t i = 0; i < sheet_id.size(); i++) {
                if (sheet_id[i] == id) {
                    return sheet_adjacent_feature_edges_count[i];
                }
            }
            return 0;
        }

        double get_curvature_metric(int id) {
            for (size_t i = 0; i < sheet_id.size(); i++) {
                if (sheet_id[i] == id) {
                    return sheet_curvature_metric[i];
                }
            }
            return 0.0;
        }

        double get_normal_variation(int id) {
            for (size_t i = 0; i < sheet_id.size(); i++) {
                if (sheet_id[i] == id) {
                    return sheet_normal_variation[i];
                }
            }
            return 0.0;
        }

        double get_dihedral_angle_deviation(int id) {
            for (size_t i = 0; i < sheet_id.size(); i++) {
                if (sheet_id[i] == id) {
                    return sheet_dihedral_angle_deviation[i];
                }
            }
            return 0.0;
        }

        double get_min_scaled_jacobian(int id) {
            for (size_t i = 0; i < sheet_id.size(); i++) {
                if (sheet_id[i] == id) {
                    return sheet_min_scaled_jacobian[i];
                }
            }
            return 1.0;
        }

        double get_avg_scaled_jacobian(int id) {
            for (size_t i = 0; i < sheet_id.size(); i++) {
                if (sheet_id[i] == id) {
                    return sheet_avg_scaled_jacobian[i];
                }
            }
            return 1.0;
        }

        // 根据sheet_id获取所有几何特征
        std::map<std::string, double> get_all_features_map(int id) {
            std::map<std::string, double> features;
            for (size_t i = 0; i < sheet_id.size(); i++) {
                if (sheet_id[i] == id) {
                    features["energy"] = sheet_energy[i];
                    features["boundary_ratio"] = sheet_on_boundary_ratio[i];
                    features["feature_ratio"] = sheet_on_feature_ratio[i];
                    features["endpoints_special"] = sheet_endpoints_are_corners_or_boundary[i] ? 1.0 : 0.0;
                    features["adjacent_features"] = static_cast<double>(sheet_adjacent_feature_edges_count[i]);
                    features["curvature_metric"] = sheet_curvature_metric[i];
                    features["normal_variation"] = sheet_normal_variation[i];
                    features["dihedral_angle_deviation"] = sheet_dihedral_angle_deviation[i];
                    features["min_scaled_jacobian"] = sheet_min_scaled_jacobian[i];
                    features["avg_scaled_jacobian"] = sheet_avg_scaled_jacobian[i];
                    return features;
                }
            }
            return features; // 返回空map，表示未找到
        }

        // 根据sheet_id获取所有几何特征作为向量
        std::vector<double> get_features_vector(int id) {
            std::vector<double> feature_vector;
            for (size_t i = 0; i < sheet_id.size(); i++) {
                if (sheet_id[i] == id) {
                    feature_vector.push_back(sheet_energy[i]);
                    feature_vector.push_back(sheet_on_boundary_ratio[i]);
                    feature_vector.push_back(sheet_on_feature_ratio[i]);
                    feature_vector.push_back(sheet_endpoints_are_corners_or_boundary[i] ? 1.0 : 0.0);
                    feature_vector.push_back(static_cast<double>(sheet_adjacent_feature_edges_count[i]));
                    feature_vector.push_back(sheet_curvature_metric[i]);
                    feature_vector.push_back(sheet_normal_variation[i]);
                    feature_vector.push_back(sheet_dihedral_angle_deviation[i]);
                    feature_vector.push_back(sheet_min_scaled_jacobian[i]);
                    feature_vector.push_back(sheet_avg_scaled_jacobian[i]);
                    return feature_vector;
                }
            }
            // 返回填充了0的向量，表示未找到
            return std::vector<double>(10, 0.0);
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
                    std::cout << "sheet id: " << sheet_id[i] << " energy: " << sheet_energy[i]
                              << " boundary_ratio: " << sheet_on_boundary_ratio[i]
                              << " feature_ratio: " << sheet_on_feature_ratio[i]
                              << " endpoints_special: " << (sheet_endpoints_are_corners_or_boundary[i] ? "true" : "false")
                              << " adjacent_features: " << sheet_adjacent_feature_edges_count[i]
                              << " curvature: " << sheet_curvature_metric[i]
                              << " normal_var: " << sheet_normal_variation[i]
                              << " dihedral_dev: " << sheet_dihedral_angle_deviation[i]
                              << " min_jacobian: " << sheet_min_scaled_jacobian[i]
                              << " avg_jacobian: " << sheet_avg_scaled_jacobian[i] << std::endl;
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
float calc_reward(int current_singularity_num, get_singularity_number<TMesh> get_singularity_num_op,
             const State& prev_state, int action_index, const State& new_state);
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
    const char* build_mode_env = std::getenv("RL_HMESH_BUILD_MODE");
    std::string build_mode = build_mode_env ? build_mode_env : "Release";
    size_t len = 0;
    _dupenv_s(&session_id_env, &len, "RL_HMESH_SESSION_ID");
    std::string session_id = session_id_env ? std::string(session_id_env) : "unknown";
    if (session_id_env) free(session_id_env);
    
    // 创建当前会话的日志目录
    std::string log_dir = "f:\\RL_HMesh\\logs\\"+ build_mode+"\\" + session_id;
    std::string log_filename = log_dir + "\\training.log";
    
    // 确保目录存在，如果不存在则创建
    std::string create_logs_dir = "mkdir \"f:\\RL_HMesh\\logs\" 2>nul";
    system(create_logs_dir.c_str());
    
    std::string create_session_dir = "mkdir \"" + log_dir + "\" 2>nul";
    system(create_session_dir.c_str());
    
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
                state.add(i, 0.0, 0.0, 0.0, false, 0, 0.0, 0.0, 0.0, 0.0, 0.0); // need_check
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
        
        // 计算几何特征
        try {
            std::cout << "Computing geometric features for sheet_id: " << sheet_id << std::endl;
            
            // 计算边界比例
            double boundary_ratio = sheet_op.get_sheet_on_boundary_ratio(sheet);
            std::cout << "  boundary_ratio: " << boundary_ratio << std::endl;
            
            // 计算特征比例
            double feature_ratio = sheet_op.get_sheet_on_feature_ratio(sheet);
            std::cout << "  feature_ratio: " << feature_ratio << std::endl;
            
            // 检查端点是否为角点或边界点
            bool endpoints_special = sheet_op.sheet_endpoints_are_corners_or_boundary(sheet);
            std::cout << "  endpoints_are_corners_or_boundary: " << (endpoints_special ? "true" : "false") << std::endl;
            
            // 计算相邻特征边的数量
            int adjacent_features = sheet_op.get_sheet_adjacent_feature_edges_count(sheet);
            std::cout << "  adjacent_feature_edges_count: " << adjacent_features << std::endl;
            
            // 计算曲率度量
            double curvature_metric = sheet_op.get_sheet_curvature_metric(sheet);
            std::cout << "  curvature_metric: " << curvature_metric << std::endl;

            // 计算法线变化
            double normal_variation = sheet_op.get_sheet_normal_variation(sheet);
            std::cout << "  normal_variation: " << normal_variation << std::endl;
            
            // 计算二面角偏差
            double dihedral_angle_deviation = sheet_op.get_sheet_dihedral_angle_deviation(sheet);
            std::cout << "  dihedral_angle_deviation: " << dihedral_angle_deviation << std::endl;
            
            // 计算相邻六面体的雅可比指标
            std::pair<double, double> jacobian_metrics = sheet_op.get_sheet_adjacent_hex_jacobian(sheet);
            std::cout << "  min_scaled_jacobian: " << jacobian_metrics.first << ", avg_scaled_jacobian: " << jacobian_metrics.second << std::endl;
            
            // 将计算出的几何特征添加到state中
            state.add(sheet_id, sheet_energy,boundary_ratio, feature_ratio, endpoints_special, adjacent_features, curvature_metric, normal_variation, dihedral_angle_deviation, jacobian_metrics.first, jacobian_metrics.second);

            
        } catch (const std::exception& e) {
            std::cerr << "Exception during feature calculation for sheet " << sheet_id << ": " << e.what() << std::endl;
            
            // 如果出现异常，仍然添加sheet但使用默认值
            state.add(sheet_id, sheet_energy, 0.0, 0.0, false, 0, 1.0, 0.0, 0.0, 1.0, 1.0);
        }
    }
    std::cout << "state size: " << state.sheet_num << std::endl;
}

// 修改play_action函数，整合日志
int play_action(int action, int done, State& state, TMesh& tmesh, sheet_operation<TMesh>& sheet_op, get_singularity_number<TMesh>& get_singularity_num_op) {
    // 不再使用固定的state size
    std::cout << "action: " << action << " state_energy: " << state.sheet_energy[action] << std::endl;
    
    /* 从环境变量获取会话ID */
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
    /* 管理logs/traininh文件 */
    
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
        std::cout << "collapse sheet: " << state.sheet_id[action]<<" done" << std::endl;

        sheet_op.get_mesh_sheet_number();
        std::cout << "get_mesh_sheet_number" << std::endl;

        get_singularity_num_op.generate_singularity_number(&tmesh);
        std::cout << "generate_singularity_number" << std::endl;

        state.clear();
        std::cout << "clear state down" << std::endl;
        calc_state(&tmesh, state, sheet_op);
        std::cout << "calc_state down" << std::endl;
        // 不再调用fill_state，让state包含所有的sheet
        if (action_log.is_open()) {
            action_log << "Done status: " << done << std::endl;
            action_log << "========== End of Action Log ==========" << std::endl << std::endl;
            action_log.close();
        }
        return done;
    }
}

float calc_reward(int current_singularity_num, get_singularity_number<TMesh> get_singularity_num_op,
             const State& prev_state, int action_index, const State& new_state) {
    // 基于奇异性数量变化的基础奖励
    float base_reward = 0.0f;
    
    std::cout << "New singularity_num: " << get_singularity_num_op.singualarity_id 
              << " current_singularity_num: " << current_singularity_num << std::endl;
              
    if (get_singularity_num_op.singualarity_id < current_singularity_num) 
        base_reward = 10.0f;
    else if (get_singularity_num_op.singualarity_id > current_singularity_num) 
        base_reward = -10.0f;
    else
        base_reward = -1.0f;
    
    // 如果没有正确的前一个状态，直接返回基础奖励
    if (action_index < 0 || action_index >= prev_state.size()) {
        return base_reward;
    }
    
    // 获取折叠的sheet的几何特征
    double boundary_ratio = prev_state.sheet_on_boundary_ratio[action_index];
    double feature_ratio = prev_state.sheet_on_feature_ratio[action_index];
    bool endpoints_special = prev_state.sheet_endpoints_are_corners_or_boundary[action_index];
    int adjacent_features = prev_state.sheet_adjacent_feature_edges_count[action_index];
    double curvature = prev_state.sheet_curvature_metric[action_index];
    double normal_variation = prev_state.sheet_normal_variation[action_index];
    double angle_deviation = prev_state.sheet_dihedral_angle_deviation[action_index];
    double min_jacobian = prev_state.sheet_min_scaled_jacobian[action_index];
    double avg_jacobian = prev_state.sheet_avg_scaled_jacobian[action_index];
    
    std::cout << "Calculating advanced reward components..." << std::endl;
    std::cout << "  boundary_ratio: " << boundary_ratio << std::endl;
    std::cout << "  feature_ratio: " << feature_ratio << std::endl;
    std::cout << "  endpoints_special: " << (endpoints_special ? "true" : "false") << std::endl;
    std::cout << "  adjacent_features: " << adjacent_features << std::endl;
    std::cout << "  curvature: " << curvature << std::endl;
    std::cout << "  normal_variation: " << normal_variation << std::endl;
    std::cout << "  angle_deviation: " << angle_deviation << std::endl;
    std::cout << "  min_jacobian: " << min_jacobian << std::endl;
    std::cout << "  avg_jacobian: " << avg_jacobian << std::endl;
    
    float geometry_reward = 0.0f;
    
    // 1. 边界和特征保护惩罚
    // 如果折叠了大量边界或特征边，给予惩罚
    if (boundary_ratio > 0.5) {
        geometry_reward -= 5.0f * boundary_ratio;
        std::cout << "  Boundary protection penalty: " << (-5.0f * boundary_ratio) << std::endl;
    }
    
    if (feature_ratio > 0.3) {
        geometry_reward -= 8.0f * feature_ratio;
        std::cout << "  Feature protection penalty: " << (-8.0f * feature_ratio) << std::endl;
    }
    
    // 2. 特殊端点保护惩罚
    if (endpoints_special) {
        geometry_reward -= 5.0f;
        std::cout << "  Special endpoints penalty: -5.0" << std::endl;
    }
    
    // 3. 曲率相关奖励
    // 鼓励折叠低曲率(接近直线)的sheet
    if (curvature > 0.8) {  // 接近直线
        geometry_reward += 3.0f * curvature;  // 增加奖励权重从2.0到3.0
        std::cout << "  High straightness bonus: " << (3.0f * curvature) << std::endl;
    }
    
    // 4. 法线变化惩罚
    // 如果sheet上法线变化大，折叠可能导致几何变形
    if (normal_variation > 20.0) {  // 降低阈值从30到20度，更严格地保护曲面特征
        geometry_reward -= 0.15f * (normal_variation - 20.0);  // 轻微增加惩罚系数
        std::cout << "  High normal variation penalty: " << (-0.15f * (normal_variation - 20.0)) << std::endl;
    }
    
    // 5. 二面角偏差惩罚
    // 如果sheet上二面角偏差大，折叠可能导致畸变
    if (angle_deviation > 15.0) {  // 降低阈值从20到15度，更严格地保护锐角特征
        geometry_reward -= 0.15f * (angle_deviation - 15.0);  // 轻微增加惩罚系数
        std::cout << "  High angle deviation penalty: " << (-0.15f * (angle_deviation - 15.0)) << std::endl;
    }
    
    // 6. Jacobian质量奖励/惩罚
    // 鼓励折叠那些周围元素已经畸变的sheet(可能改善质量)
    if (min_jacobian < 0.3) {  // 畸变元素
        geometry_reward += 4.0f * (1.0 - min_jacobian);  // 增加对于改善低质量区域的奖励
        std::cout << "  Low quality neighborhood bonus: " << (4.0f * (1.0 - min_jacobian)) << std::endl;
    } else {
        // 避免折叠高质量区域
        geometry_reward -= min_jacobian * 2.0f;
        std::cout << "  High quality neighborhood penalty: " << (-min_jacobian * 2.0f) << std::endl;
    }
    
    // 7. 网格总体质量变化奖励
    // 计算折叠前后的平均雅可比值变化
    double avg_jacobian_before = 0.0;
    double avg_jacobian_after = 0.0;
    int count_before = 0;
    int count_after = 0;
    
    for (size_t i = 0; i < prev_state.sheet_avg_scaled_jacobian.size(); i++) {
        avg_jacobian_before += prev_state.sheet_avg_scaled_jacobian[i];
        count_before++;
    }
    
    for (size_t i = 0; i < new_state.sheet_avg_scaled_jacobian.size(); i++) {
        avg_jacobian_after += new_state.sheet_avg_scaled_jacobian[i];
        count_after++;
    }
    
    if (count_before > 0 && count_after > 0) {
        avg_jacobian_before /= count_before;
        avg_jacobian_after /= count_after;
        
        float quality_change = (float)(avg_jacobian_after - avg_jacobian_before);
        geometry_reward += quality_change * 15.0f;  // 增加网格总体质量变化的重要性
        std::cout << "  Quality change reward: " << (quality_change * 15.0f) << std::endl;
    }
    
    // 组合基础奖励和几何奖励
    float total_reward = base_reward + geometry_reward;
    
    std::cout << "Base reward: " << base_reward << std::endl;
    std::cout << "Geometry reward: " << geometry_reward << std::endl;
    std::cout << "Total reward: " << total_reward << std::endl;
    
    return total_reward;
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
        // 添加 sheet_id
        row.append(state.sheet_id[i]);
        // 添加 sheet_energy
        row.append(state.sheet_energy[i]);
        // 添加所有几何特征
        row.append(state.sheet_on_boundary_ratio[i]);
        row.append(state.sheet_on_feature_ratio[i]);
        row.append(state.sheet_endpoints_are_corners_or_boundary[i] ? 1.0 : 0.0);
        row.append(static_cast<double>(state.sheet_adjacent_feature_edges_count[i]));
        row.append(state.sheet_curvature_metric[i]);
        row.append(state.sheet_normal_variation[i]);
        row.append(state.sheet_dihedral_angle_deviation[i]);
        row.append(state.sheet_min_scaled_jacobian[i]);
        row.append(state.sheet_avg_scaled_jacobian[i]);
        result.append(row);
    }
    return result;
}

#endif // STATE_FUNCTIONS_H
