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
#include "tools.h"

#include <algorithm> // 用于 std::min, std::max
#include <limits>    // 用于 std::numeric_limits
#include <cmath>     // 用于 std::abs 和其他数学函数

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

using namespace HMeshLib;

// ──────────────────── 归一化相关定义 ────────────────────

// 定义单个特征的归一化参数结构体
struct FeatureNormalizationParams {
    double min_val;
    double max_val;
};

// `sheet_energy` 的特殊标记值
const double ENERGY_CANNOT_COLLAPSE_SENTINEL = -99999; // 与 sheet_operation.h 中 predict_sheet_collapse_energy 的返回值一致
// `sheet_energy` 正值部分的缩放参数
const double ENERGY_MIN_POSITIVE_FOR_SCALING = 0.001;   // 用于缩放的最小正能量值 (避免log(0)或除以极小值)
const double ENERGY_MAX_POSITIVE_FOR_SCALING = 500.0;   // 用于缩放的最大正能量值 (作为上限，可根据数据调整，例如99百分位数)

// `sheet_energy` 归一化后的特殊目标值
const double NORMALIZED_ENERGY_COLLAPSED = 0.0;          // 已坍缩的Sheet (原始能量0.0) 映射为此值
const double NORMALIZED_ENERGY_CANNOT_COLLAPSE = -0.1;    // 不可坍缩的Sheet (原始能量-99999.0) 映射为此值
// `sheet_energy` 正值部分归一化后的目标范围
const double NORMALIZED_ENERGY_POSITIVE_TARGET_MIN = 0.0;
const double NORMALIZED_ENERGY_POSITIVE_TARGET_MAX = 1.0;

// 其他特征的归一化参数数组。顺序必须与 state_to_list 中添加特征的顺序一致。
const FeatureNormalizationParams normalization_params_others[9] = {
    // 对应特征:                                   占位符 min,占位符 max
    {0.0, 1.0},    // sheet_on_boundary_ratio
    {0.0, 1.0},    // sheet_on_feature_ratio
    {0.0, 15.0},   // sheet_adjacent_feature_edges_count
    {0.0, 1.0},    // sheet_curvature_metric
    {0.0, 180.0},  // sheet_normal_variation
    {0.0, 180.0},  // sheet_dihedral_angle_deviation
    {-1.0, 1.0},   // sheet_min_scaled_jacobian
    {-1.0, 1.0},   // sheet_avg_scaled_jacobian
    {0.0, 100.0}   // sheet_distance_to_feature
};

// State 类定义
class State {
    public:
        std::vector<int> sheet_id;
        std::vector<double> sheet_energy;
    
        // --- 几何/拓扑特征 ---
        std::vector<double> sheet_on_boundary_ratio;
        std::vector<double> sheet_on_feature_ratio;
        //std::vector<bool> sheet_endpoints_are_corners_or_boundary;
        std::vector<int> sheet_adjacent_feature_edges_count;
        std::vector<double> sheet_curvature_metric;
        std::vector<double> sheet_normal_variation;
        std::vector<double> sheet_dihedral_angle_deviation;
        std::vector<double> sheet_min_scaled_jacobian;
        std::vector<double> sheet_avg_scaled_jacobian;
        //std::vector<int> sheet_valence_along_path_count;
        // --- 新增距离属性 ---
        //std::vector<double> sheet_distance_to_boundary; // 新增
        std::vector<double> sheet_distance_to_feature;  // 新增
    
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
            //sheet_endpoints_are_corners_or_boundary.clear();
            sheet_adjacent_feature_edges_count.clear();
            sheet_curvature_metric.clear();
            sheet_normal_variation.clear();
            sheet_dihedral_angle_deviation.clear();
            sheet_min_scaled_jacobian.clear();
            sheet_avg_scaled_jacobian.clear();
            //sheet_valence_along_path_count.clear();
            //sheet_distance_to_boundary.clear(); // 新增
            sheet_distance_to_feature.clear();  // 新增
            sheet_num = 0;
        }
        // 返回状态大小
        size_t size() const {
            return sheet_id.size();
        }
        // 添加一条记录 (更新参数列表)
        void add(int id, double energy, double boundary_ratio,
                 double feature_ratio,
                // bool endpoints_special,
                 int adjacent_features,
                 double curvature_metric,
                 double normal_variation,
                 double dihedral_deviation,
                 double min_jacobian,
                 double avg_jacobian,
                 //int valence_along_path,
                 //double distance_boundary, // 新增
                 double distance_feature   // 新增
                 )
        {
            sheet_id.push_back(id);
            sheet_energy.push_back(energy);
            sheet_on_boundary_ratio.push_back(boundary_ratio);
            sheet_on_feature_ratio.push_back(feature_ratio);
            //sheet_endpoints_are_corners_or_boundary.push_back(endpoints_special);
            sheet_adjacent_feature_edges_count.push_back(adjacent_features);
            sheet_curvature_metric.push_back(curvature_metric);
            sheet_normal_variation.push_back(normal_variation);
            sheet_dihedral_angle_deviation.push_back(dihedral_deviation);
            sheet_min_scaled_jacobian.push_back(min_jacobian);
            sheet_avg_scaled_jacobian.push_back(avg_jacobian);
            //sheet_valence_along_path_count.push_back(valence_along_path);
            //sheet_distance_to_boundary.push_back(distance_boundary); // 新增
            sheet_distance_to_feature.push_back(distance_feature);   // 新增
            sheet_num++;
        }
    
        // 重载赋值操作符 (更新)
        State& operator=(const State& other) {
            if (this != &other) {
                sheet_id = other.sheet_id;
                sheet_energy = other.sheet_energy;
                sheet_on_boundary_ratio = other.sheet_on_boundary_ratio;
                sheet_on_feature_ratio = other.sheet_on_feature_ratio;
                //sheet_endpoints_are_corners_or_boundary = other.sheet_endpoints_are_corners_or_boundary;
                sheet_adjacent_feature_edges_count = other.sheet_adjacent_feature_edges_count;
                sheet_curvature_metric = other.sheet_curvature_metric;
                sheet_normal_variation = other.sheet_normal_variation;
                sheet_dihedral_angle_deviation = other.sheet_dihedral_angle_deviation;
                sheet_min_scaled_jacobian = other.sheet_min_scaled_jacobian;
                sheet_avg_scaled_jacobian = other.sheet_avg_scaled_jacobian;
                //sheet_valence_along_path_count = other.sheet_valence_along_path_count;
                //sheet_distance_to_boundary = other.sheet_distance_to_boundary; // 新增
                sheet_distance_to_feature = other.sheet_distance_to_feature;   // 新增
                sheet_num = other.sheet_num;
            }
            return *this;
        }
    
    
        // 更新 get_all_features_map (更新)
        std::map<std::string, double> get_all_features_map(int id) const {
            std::map<std::string, double> features;
            for (size_t i = 0; i < sheet_id.size(); i++) {
                if (sheet_id[i] == id) {

                        features["energy"] = sheet_energy[i];
                        features["boundary_ratio"] = sheet_on_boundary_ratio[i];
                        features["feature_ratio"] = sheet_on_feature_ratio[i];
                        //features["endpoints_special"] = sheet_endpoints_are_corners_or_boundary[i] ? 1.0 : 0.0;
                        features["adjacent_features"] = static_cast<double>(sheet_adjacent_feature_edges_count[i]);
                        features["curvature_metric"] = sheet_curvature_metric[i];
                        features["normal_variation"] = sheet_normal_variation[i];
                        features["dihedral_angle_deviation"] = sheet_dihedral_angle_deviation[i];
                        features["min_scaled_jacobian"] = sheet_min_scaled_jacobian[i];
                        features["avg_scaled_jacobian"] = sheet_avg_scaled_jacobian[i];
                        //features["valence_along_path"] = static_cast<double>(sheet_valence_along_path_count[i]);
                        //features["distance_to_boundary"] = sheet_distance_to_boundary[i]; // 新增
                        features["distance_to_feature"] = sheet_distance_to_feature[i];   // 新增
                    } 
                    return features; // Return after finding the id
                }
            }
    
       // 更新 get_features_vector (移除对应条目)
       std::vector<double> get_features_vector(int id) const {
        std::vector<double> feature_vector;
        size_t expected_remaining_features = 10; // 13 - 3 = 10

        for (size_t i = 0; i < sheet_id.size(); i++) {
            if (sheet_id[i] == id) {
                // 检查索引是否有效 (根据保留的向量)
                 if (i < sheet_energy.size() && i < sheet_on_boundary_ratio.size() &&
                     i < sheet_on_feature_ratio.size() && i < sheet_adjacent_feature_edges_count.size() &&
                     i < sheet_curvature_metric.size() && i < sheet_normal_variation.size() &&
                     i < sheet_dihedral_angle_deviation.size() && i < sheet_min_scaled_jacobian.size() &&
                     i < sheet_avg_scaled_jacobian.size() && i < sheet_distance_to_feature.size())
                 {
                    feature_vector.push_back(sheet_energy[i]);
                    feature_vector.push_back(sheet_on_boundary_ratio[i]);
                    feature_vector.push_back(sheet_on_feature_ratio[i]);
                    //feature_vector.push_back(sheet_endpoints_are_corners_or_boundary[i] ? 1.0 : 0.0); // 已移除
                    feature_vector.push_back(static_cast<double>(sheet_adjacent_feature_edges_count[i]));
                    feature_vector.push_back(sheet_curvature_metric[i]);
                    feature_vector.push_back(sheet_normal_variation[i]);
                    feature_vector.push_back(sheet_dihedral_angle_deviation[i]);
                    feature_vector.push_back(sheet_min_scaled_jacobian[i]);
                    feature_vector.push_back(sheet_avg_scaled_jacobian[i]);
                    //feature_vector.push_back(static_cast<double>(sheet_valence_along_path_count[i])); // 已移除
                    //feature_vector.push_back(sheet_distance_to_boundary[i]); // 已移除
                    feature_vector.push_back(sheet_distance_to_feature[i]);
                } else {
                    std::cerr << "Error: Index mismatch in State vectors for id " << id << " in get_features_vector (after removal)" << std::endl;
                    // 返回填充了 0 的向量
                    return std::vector<double>(expected_remaining_features, 0.0);
                }
                return feature_vector; // 找到后即可返回
            }
        }
        // 未找到则返回填充了 0 的向量
        return std::vector<double>(expected_remaining_features, 0.0);
    }
    
            // 检查state中除了当前action外是否还有其他sheet_energy > 0的情况, 有则返回true,此时为异常
        bool check_error_choice(int action) {
            for (size_t i = 0; i < sheet_id.size(); i++) {
                    // Check if the index is valid for sheet_energy
                if (i >= sheet_energy.size()) {
                        std::cerr << "Error: Index out of bounds in check_error_choice for sheet_energy" << std::endl;
                        continue; // Skip this entry
                }
                // Find the index corresponding to the action ID
                size_t action_idx = -1;
                for(size_t j=0; j < sheet_id.size(); ++j){
                    if(sheet_id[j] == action){
                        action_idx = j;
                        break;
                    }
                }
                // Check if action id was found and if the current index is not the action index
                if(action_idx != -1 && i != action_idx && sheet_energy[i] <= 0) {
                    return true;
                } else if (action_idx == -1 && sheet_energy[i] <= 0) {
                    // If action ID wasn't found (shouldn't happen if action is valid),
                    // check if any *other* energy is non-positive
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

        // 新增：计算当前状态的平均雅可比值
        double get_average_jacobian() const {
            if (sheet_avg_scaled_jacobian.empty()) {
                return 0.0; // 或者返回一个表示无效的值，例如 -1.0
            }
            double sum = 0.0;
            for (double jacobian : sheet_avg_scaled_jacobian) {
                sum += jacobian;
            }
            return sum / sheet_avg_scaled_jacobian.size();
        }

        // 更新 print (更新)
        void print() const {
            if (sheet_id.empty()) {
                std::cout << "state is empty" << std::endl;
            } else {
                std::cout << "state: (" << size() << " entries)" << std::endl;
                for (size_t i = 0; i < size(); i++) {
    
                    std::cout << "  sheet id: " << sheet_id[i]
                              << ", energy: " << sheet_energy[i]
                              << ", b_ratio: " << sheet_on_boundary_ratio[i]
                              << ", f_ratio: " << sheet_on_feature_ratio[i]
                              //<< ", end_sp: " << (sheet_endpoints_are_corners_or_boundary[i] ? "T" : "F")
                              << ", adj_f: " << sheet_adjacent_feature_edges_count[i]
                              << ", curv: " << sheet_curvature_metric[i]
                              << ", norm_var: " << sheet_normal_variation[i]
                              << ", dih_dev: " << sheet_dihedral_angle_deviation[i]
                              << ", min_jac: " << sheet_min_scaled_jacobian[i]
                              << ", avg_jac: " << sheet_avg_scaled_jacobian[i]
                              //<< ", valence: " << sheet_valence_along_path_count[i]
                              //<< ", dist_bound: " << sheet_distance_to_boundary[i] // 新增
                              << ", dist_feat: " << sheet_distance_to_feature[i]  // 新增
                              << std::endl;
                }
            }
        }
    
    }; // End State class

// ────────────────────
std::string format_matrix(const State& state);
void log(int episode, const State& state, int action, const State& next_state, float reward);
void fill_state(State& state, int origin_size);
int find_sheet_id(State state, int action);
bool check_in_state(State state, int action);
std::vector<TE*> get_sheet_byId(TMesh* mesh, int sheet_id, sheet_operation<TMesh>& sheet_op);
void calc_state(TMesh* mesh, State& state, sheet_operation<TMesh>& sheet_op);
int play_action(int action, int done, State& state, TMesh& tmesh, sheet_operation<TMesh>& sheet_op, get_singularity_number<TMesh>& get_singularity_num_op, int original_hex_count);
float calc_reward(int current_singularity_num, get_singularity_number<TMesh> get_singularity_num_op,
             const State& prev_state, int action_index, const State& new_state);
void print_state(const State& state);
void compare_diff(const State& state1, const State& state2);
py::list state_to_list(const State& state);

double calculate_average_min_jacobian(TMesh* mesh);
float evaluate_overall_mesh_quality(TMesh* mesh, get_singularity_number<TMesh>& singularity_checker);
double calculate_global_min_jacobian(TMesh* mesh);

double scale_value(double value, double source_min, double source_max, double target_min, double target_max);
double normalize_sheet_energy(double energy_value);
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

// --- Update fill_state ---
void fill_state(State& state, int origin_size) {
    const double DEFAULT_DISTANCE = 1e6; // Use a consistent default distance
   if (state.size() < static_cast<size_t>(origin_size)) {
       for (int i = 1; i <= origin_size; i++) {
           if (!state.check_sheet_id(i)) {
               // Add with default values, including the new distance and valence
               state.add(i, 0.0, 0.0, 0.0,  0, 1.0, 0.0, 0.0, 1.0, 1.0, DEFAULT_DISTANCE);
           }
       }
   } else {
      // std::cout << "fill_state: state is already full or larger." << std::endl;
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

// --- Update calc_state ---
void calc_state(TMesh* mesh, State& state, sheet_operation<TMesh>& sheet_op) {
    auto logFunc = HMeshLib::getSheetLog() ? HMeshLib::log : [](const std::string&){}; // Use sheet_operation log

    //logFunc("calc_state: Starting state calculation.");
    //logFunc("calc_state: Current sharp edge count: " + std::to_string(mesh->count_sharp_edges()));

    // 获取最大的sheet ID
    int max_sheet_id = 0;
    for (TMesh::MEIterator eite(mesh); !eite.end(); eite++) {
         TE* be = *eite;
         if(be) max_sheet_id = std::max(max_sheet_id, be->sheet());
    }
    //logFunc("calc_state: Maximum sheet ID found: " + std::to_string(max_sheet_id));

    // state.print(); // Can print initial empty state if needed

    for (int sheet_id = 1; sheet_id <= max_sheet_id; sheet_id++) {
        TE* edge_with_id = nullptr;
        for (TMesh::MEIterator eite(mesh); !eite.end(); eite++) {
             TE* be = *eite;
             if (be && be->sheet() == sheet_id) {
                 edge_with_id = be;
                 break;
             }
        }
        if (!edge_with_id){
            // logFunc("calc_state: No edge found for sheet ID " + std::to_string(sheet_id) + ". Skipping.");
             continue;
        }

        std::vector<TE*> sheet = sheet_op.get_one_sheet(edge_with_id);
        if (sheet.empty()) {
            //logFunc("calc_state: Warning - get_one_sheet returned empty for sheet ID " + std::to_string(sheet_id));
            continue;
        }

        double sheet_energy = sheet_op.predict_sheet_collapse_energy(sheet);
        // logFunc("calc_state: Calculated energy for sheet ID " + std::to_string(sheet_id) + ": " + std::to_string(sheet_energy));

        try {
            // logFunc("calc_state: Computing geometric features for sheet ID " + std::to_string(sheet_id));
            // std::cout<< "boundary count by face in state func: " << mesh->count_boundary_byF() << std::endl;
            // std::cout<< "boundary count by edge in state func: " << mesh->count_boundary_byE() << std::endl;
            // std::cout<< "sharp edge count in state func: " << mesh->count_sharp_edges() << std::endl;
            // std::cout<< "corner count in state func: " << mesh->count_corner() << std::endl;
            double boundary_ratio = sheet_op.get_sheet_on_boundary_ratio(sheet);
            double feature_ratio = sheet_op.get_sheet_on_feature_ratio(sheet);
            //bool endpoints_special = sheet_op.sheet_endpoints_are_corners_or_boundary(sheet);
            int adjacent_features = sheet_op.get_sheet_adjacent_feature_edges_count(sheet);
            double curvature_metric = sheet_op.get_sheet_curvature_metric(sheet);
            double normal_variation = sheet_op.get_sheet_normal_variation(sheet);
            double dihedral_deviation = sheet_op.get_sheet_dihedral_angle_deviation(sheet);
            std::pair<double, double> jacobian_metrics = sheet_op.get_sheet_adjacent_hex_jacobian(sheet);
            //int valence_along_path = sheet_op.get_sheet_valence_along_path(sheet);

            //double dist_bound = sheet_op.get_distance_to_boundary(sheet); 
            double dist_feat;
            if (feature_ratio >0) {
                // Sheet 本身是特征，距离为 0
                dist_feat = 0.0;
            } else {
                // Sheet 不是特征，调用修改后的距离计算函数
                // 注意：现在 get_distance_to_feature 内部不再检查 sheet 是否 sharp
                dist_feat = sheet_op.get_distance_to_feature(sheet);
            } 

            // if(feature_ratio==0 && dist_feat ==0){
            //     std::cout << "sheet_id: " << sheet_id << " feature_ratio: " << feature_ratio << " dist_feat: " << dist_feat << std::endl;
            //     // print sheet 中所有edge的is sharp情况
            //     for (TE* edge : sheet) {
            //         if (edge->sharp()) {
            //             std::cout << "Edge ID: " << edge->id() << " is sharp." << std::endl;
            //         } else {
            //             std::cout << "Edge ID: " << edge->id() << " is not sharp." << std::endl;
            //         }
            //     }
            //     std::cin.get(); // Wait for user input to pause
            // }

            // 调用更新后的 add
            state.add(sheet_id, sheet_energy, boundary_ratio, feature_ratio,
                       adjacent_features, curvature_metric,
                      normal_variation, dihedral_deviation,
                      jacobian_metrics.first, jacobian_metrics.second,
                      dist_feat); // 添加新参数

        } catch (const std::exception& e) {
            //logFunc("calc_state: Exception during feature calculation for sheet " + std::to_string(sheet_id) + ": " + std::string(e.what()));
            // 使用默认值添加 (包括新属性的默认值) - 使用大距离表示无效或未计算
             const double DEFAULT_DISTANCE = 1e6; // Or other suitable default
            state.add(sheet_id, sheet_energy, 0.0, 0.0, 0, 1.0, 0.0, 0.0, 1.0, 1.0, DEFAULT_DISTANCE);
        }
    }
    //logFunc("calc_state: State calculation finished. Final state size: " + std::to_string(state.size()));
    state.print();
}

// 修改play_action函数，整合日志
int play_action(int action, int done, State& state, TMesh& tmesh, sheet_operation<TMesh>& sheet_op, get_singularity_number<TMesh>& get_singularity_num_op,int original_hex_count) {
    // 不再使用固定的state size
    std::cout << "action: " << action << " state_energy: " << state.sheet_energy[action] << std::endl;
    
    /* 从环境变量获取会话ID */
    // char* session_id_env = nullptr;
    // size_t len = 0;
    // _dupenv_s(&session_id_env, &len, "RL_HMESH_SESSION_ID");
    // std::string session_id = session_id_env ? std::string(session_id_env) : "unknown";
    // if (session_id_env) free(session_id_env);
    
    // // 创建当前会话的日志目录
    // std::string log_dir = "f:\\RL_HMesh\\logs\\" + session_id;
    
    // // 创建特定于此操作的日志文件
    // std::string action_log_filename = log_dir + "\\action_" + 
    //                                 std::to_string(action) + "_sheet_" + 
    //                                 std::to_string(state.sheet_id[action]) + ".log";
    
    // std::ofstream action_log(action_log_filename, std::ios::app);
    
    //if (action_log.is_open()) {
    //     auto now = std::chrono::system_clock::now();
    //     std::time_t now_time = std::chrono::system_clock::to_time_t(now);
    //     std::tm now_tm;
    //     localtime_s(&now_tm, &now_time);
        
    //     char timestamp[64];
    //     std::strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", &now_tm);
        
    //     action_log << "==== Action Log: " << timestamp << " ====" << std::endl;
    //     action_log << "Session ID: " << session_id << std::endl;
    //     action_log << "Action: " << action << " Sheet ID: " << state.sheet_id[action] << " Energy: " << state.sheet_energy[action] << std::endl;
    // }
    /* 管理logs/traininh文件 */
    
    if (state.sheet_energy[action] <= 0 && state.check_error_choice(action)) {
        std::cout << "sheet_energy is 0, can't collapse" << std::endl;
        done = 0; // 0表示 因为模型选择错误导致的异常，应赋予极大惩罚
        // if (action_log.is_open()) {
        //     action_log << "Done status: " << done << std::endl;
        //     action_log << "========== End of Action Log ==========" << std::endl << std::endl;
        //     action_log.close();
        // }
        return done;
    } 
    else {

        int before_action_mesh_size = tmesh.hs.size(); // 记录折叠前的网格大小

        sheet_op.collapse_one_sheet2(get_sheet_byId(&tmesh, state.sheet_id[action], sheet_op));
        std::cout << "collapse sheet: " << state.sheet_id[action]<<" done" << std::endl;

        sheet_op.get_mesh_sheet_number();
        std::cout << "get_mesh_sheet_number" << std::endl;

        get_singularity_num_op.generate_singularity_number(&tmesh);
        std::cout << "generate_singularity_number" << std::endl;

        tmesh.compute_features();

        int current_mesh_size = tmesh.hs.size(); // 记录折叠后的网格大小

        if(current_mesh_size < 0.6*original_hex_count) {   //0.x 是一个阈值，可以根据需要进行调整
            std::cout << "mesh size is less than 0.5*original" << std::endl;
            done = 1; // 0表示 因为模型选择错误导致的异常，应赋予极大惩罚
            return done;
        }

        done=2;

        state.clear();
        std::cout << "clear state down" << std::endl;
        calc_state(&tmesh, state, sheet_op);
        std::cout << "calc_state down" << std::endl;
        // 不再调用fill_state，让state包含所有的sheet
        // if (action_log.is_open()) {
        //     action_log << "Done status: " << done << std::endl;
        //     action_log << "========== End of Action Log ==========" << std::endl << std::endl;
        //     action_log.close();
        // }
        return done;
    }
}

// --- Update calc_reward if needed ---
float calc_reward(int current_singularity_num, get_singularity_number<TMesh> get_singularity_num_op,
    const State& prev_state, int action_index, const State& new_state) {

    // auto logFunc = HMeshLib::getSheetLog() ? HMeshLib::log : [](const std::string&){};

    // logFunc("calc_reward: Calculating reward...");
    // logFunc("calc_reward: Prev singularity count = " + std::to_string(current_singularity_num) +
    // ", New singularity count = " + std::to_string(get_singularity_num_op.singualarity_id));

    // 基于奇异性数量变化的基础奖励
    float base_reward = 0.0f;
    if (get_singularity_num_op.singualarity_id < current_singularity_num)
    base_reward = 10.0f; // 正奖励: 减少奇异线
    else if (get_singularity_num_op.singualarity_id > current_singularity_num)
    base_reward = -10.0f; // 负奖励: 增加奇异线 (不期望发生)
    else
    base_reward = -1.0f; // 小惩罚: 奇异线数量不变 (鼓励进步)
    //logFunc("calc_reward: base_reward = " + std::to_string(base_reward));


    float geometry_reward = 0.0f;

    // 检查 action_index 对于 prev_state 是否有效
    // 需要检查所有将要访问的 prev_state 向量
    if (action_index < 0 || action_index >= prev_state.size() ||
    action_index >= prev_state.sheet_on_boundary_ratio.size() ||
    action_index >= prev_state.sheet_on_feature_ratio.size() ||
    //action_index >= prev_state.sheet_endpoints_are_corners_or_boundary.size() ||
    action_index >= prev_state.sheet_adjacent_feature_edges_count.size() ||
    action_index >= prev_state.sheet_curvature_metric.size() ||
    action_index >= prev_state.sheet_normal_variation.size() ||
    action_index >= prev_state.sheet_dihedral_angle_deviation.size() ||
    action_index >= prev_state.sheet_min_scaled_jacobian.size() ||
    action_index >= prev_state.sheet_avg_scaled_jacobian.size() ||
    //action_index >= prev_state.sheet_valence_along_path_count.size() || // 检查新向量
    //action_index >= prev_state.sheet_distance_to_boundary.size() ||
    action_index >= prev_state.sheet_distance_to_feature.size())
    {
    //logFunc("calc_reward: Warning - Invalid action_index (" + std::to_string(action_index) + ") or inconsistent prev_state size (" + std::to_string(prev_state.size()) + "). Using only base reward.");
    return base_reward; // 如果索引无效，只返回基础奖励
    }


    // 获取被折叠的 sheet 的几何和拓扑特征 (来自 prev_state)
    double boundary_ratio = prev_state.sheet_on_boundary_ratio[action_index];
    double feature_ratio = prev_state.sheet_on_feature_ratio[action_index];
    //bool endpoints_special = prev_state.sheet_endpoints_are_corners_or_boundary[action_index];
    int adjacent_features = prev_state.sheet_adjacent_feature_edges_count[action_index]; // 可选
    double curvature = prev_state.sheet_curvature_metric[action_index];
    double normal_variation = prev_state.sheet_normal_variation[action_index];
    double angle_deviation = prev_state.sheet_dihedral_angle_deviation[action_index];
    double min_jacobian = prev_state.sheet_min_scaled_jacobian[action_index];
    //double dist_boundary = prev_state.sheet_distance_to_boundary[action_index];
    double dist_feature = prev_state.sheet_distance_to_feature[action_index];
    //int valence = prev_state.sheet_valence_along_path_count[action_index]; // 获取 valence 值


    // --- 计算 geometry_reward ---

    // 1. 边界/特征保护 (惩罚项)
    if (boundary_ratio > 0.5) geometry_reward -= 5.0f * boundary_ratio;
    if (feature_ratio > 0.3) geometry_reward -= 8.0f * feature_ratio;
    //if (endpoints_special) geometry_reward -= 5.0f;

    // 2. 形状保持 (奖惩项)
    if (curvature > 0.8) geometry_reward += 3.0f * curvature;
    if (normal_variation > 20.0) geometry_reward -= 0.15f * (normal_variation - 20.0);
    if (angle_deviation > 15.0) geometry_reward -= 0.15f * (angle_deviation - 15.0);

    // 3. 质量考虑 (奖惩项)
    if (min_jacobian < 0.3) geometry_reward += 4.0f * (1.0 - min_jacobian);
    else geometry_reward -= min_jacobian * 2.0f;

    // 4. 距离惩罚
    //const double boundary_distance_threshold = 0.1;
    const double feature_distance_threshold = 0.05;
    const float distance_penalty_factor = 3.0f;
    if (boundary_ratio < 1e-6) { // 仅对内部 sheet 应用距离惩罚
    // if (dist_boundary < boundary_distance_threshold) {
    //     geometry_reward -= distance_penalty_factor * (1.0 - dist_boundary / boundary_distance_threshold);
    // }
    if (dist_feature < feature_distance_threshold) {
        geometry_reward -= distance_penalty_factor * (1.0 - dist_feature / feature_distance_threshold);
    }
    }

    // --- 5. 添加 Valence 惩罚逻辑 ---
    // 惩罚折叠作为主要连接点的 sheet (高 valence)
    // 阈值和惩罚因子可能需要调整
    // const int VALENCE_THRESHOLD = 2;
    // const float VALENCE_PENALTY_FACTOR = 1.0f; // 例如，每个超过阈值的连接点扣 1 分
    // if (valence > VALENCE_THRESHOLD) {
    // geometry_reward -= static_cast<float>(valence - VALENCE_THRESHOLD) * VALENCE_PENALTY_FACTOR; // 可以惩罚超出阈值的部分
    // logFunc("  calc_reward: Penalty added for high valence (" + std::to_string(valence) + ", Threshold=" + std::to_string(VALENCE_THRESHOLD) + "): " + std::to_string(-static_cast<float>(valence - VALENCE_THRESHOLD) * VALENCE_PENALTY_FACTOR));
    // }
    // --- 结束 Valence 惩罚 ---

    // 6. 整体网格质量提升奖励
    double avg_jacobian_before = prev_state.get_average_jacobian();
    double avg_jacobian_after = new_state.get_average_jacobian();
    // 检查平均雅可比值是否有效（例如，是否已计算或非零）
    if (avg_jacobian_before > -1e9 && avg_jacobian_after > -1e9 && prev_state.size() > 0 && new_state.size() > 0) {
    float quality_change = (float)(avg_jacobian_after - avg_jacobian_before);
    geometry_reward += quality_change * 15.0f; // 调整权重因子以平衡影响
    //logFunc("  calc_reward: Quality change reward: " + std::to_string(quality_change * 15.0f) + " (Before: " + std::to_string(avg_jacobian_before) + ", After: " + std::to_string(avg_jacobian_after) + ")");
    } else {
    //logFunc("  calc_reward: Skipping quality change reward due to invalid average Jacobian values or empty states.");
    }

    // 网格大小减小奖励


    // --- 总奖励 ---
    float total_reward = base_reward + geometry_reward;
    //logFunc("calc_reward: Final Reward Calculation -> Base=" + std::to_string(base_reward) + ", Geo=" + std::to_string(geometry_reward) + ", Total=" + std::to_string(total_reward));

    return total_reward;
}

void print_state(const State& state) {
    state.print(); // Use the updated print method in the State class
}

// --- 更新 state_to_list 函数以包含归一化 ---
inline py::list state_to_list(const State& state) {
    py::list result;

    for (size_t i = 0; i < state.size(); i++) {
        py::list row;
        row.append(state.sheet_id[i]); // ID 通常不参与归一化

        // 特征 0: sheet_energy (特殊处理)
        row.append(normalize_sheet_energy(state.sheet_energy[i]));

        // 其他特征 (索引从0到8，对应 normalization_params_others 数组)
        row.append(scale_value(state.sheet_on_boundary_ratio[i], normalization_params_others[0].min_val, normalization_params_others[0].max_val, 0.0, 1.0));
        row.append(scale_value(state.sheet_on_feature_ratio[i], normalization_params_others[1].min_val, normalization_params_others[1].max_val, 0.0, 1.0));
        row.append(scale_value(static_cast<double>(state.sheet_adjacent_feature_edges_count[i]), normalization_params_others[2].min_val, normalization_params_others[2].max_val, 0.0, 1.0));
        row.append(scale_value(state.sheet_curvature_metric[i], normalization_params_others[3].min_val, normalization_params_others[3].max_val, 0.0, 1.0));
        row.append(scale_value(state.sheet_normal_variation[i], normalization_params_others[4].min_val, normalization_params_others[4].max_val, 0.0, 1.0));
        row.append(scale_value(state.sheet_dihedral_angle_deviation[i], normalization_params_others[5].min_val, normalization_params_others[5].max_val, 0.0, 1.0));
        row.append(scale_value(state.sheet_min_scaled_jacobian[i], normalization_params_others[6].min_val, normalization_params_others[6].max_val, 0.0, 1.0));
        row.append(scale_value(state.sheet_avg_scaled_jacobian[i], normalization_params_others[7].min_val, normalization_params_others[7].max_val, 0.0, 1.0));
        row.append(scale_value(state.sheet_distance_to_feature[i], normalization_params_others[8].min_val, normalization_params_others[8].max_val, 0.0, 1.0));

        result.append(row);
    }
    return result;
}

double calculate_average_min_jacobian(TMesh* mesh) {
    if (!mesh || mesh->hs.empty()) {
        return 0.0; // 或者一个表示无效的值，例如负数
    }
    CMeshQuality<TMesh> quality_evaluator(mesh);
    double total_min_jacobian = 0.0;
    int valid_hex_count = 0;

    for (TMesh::MHIterator hi(mesh); !hi.end(); hi++) {
        TMesh::H* h = *hi;
        if (!h) continue;
        std::vector<double> jacobians = quality_evaluator.get_JacobianMatricesDet(h); // 获取9个雅可比值
        if (jacobians.empty()) continue;

        double min_hex_jacobian = jacobians[0]; // 假设第一个是有效的
        for (size_t i = 1; i < jacobians.size(); ++i) {
            if (jacobians[i] < min_hex_jacobian) {
                min_hex_jacobian = jacobians[i];
            }
        }
        // 排除非常小的或者负的雅可比值（这些可能表示退化或翻转的单元）
        // 一个更鲁棒的方法是只考虑正的雅可比值
        if (min_hex_jacobian > 1e-6) { // 阈值可以调整
             total_min_jacobian += min_hex_jacobian;
             valid_hex_count++;
        }
    }

    if (valid_hex_count == 0) {
        return 0.0; // 没有有效的六面体来评估质量
    }
    return total_min_jacobian / valid_hex_count;
}

// (可选) 计算网格的最小雅可比行列式 (全局最小)
double calculate_global_min_jacobian(TMesh* mesh) {
    if (!mesh || mesh->hs.empty()) {
        return 0.0;
    }
    CMeshQuality<TMesh> quality_evaluator(mesh);
    double global_min_jacobian = std::numeric_limits<double>::max();
    bool found_valid = false;

    for (TMesh::MHIterator hi(mesh); !hi.end(); hi++) {
        TMesh::H* h = *hi;
        if (!h) continue;
        std::vector<double> jacobians = quality_evaluator.get_JacobianMatricesDet(h);
        if (jacobians.empty()) continue;

        double min_hex_jacobian = jacobians[0];
        for (size_t i = 1; i < jacobians.size(); ++i) {
            if (jacobians[i] < min_hex_jacobian) {
                min_hex_jacobian = jacobians[i];
            }
        }
        if (min_hex_jacobian > 1e-6) { // 只考虑有效的正雅可比
            if (min_hex_jacobian < global_min_jacobian) {
                global_min_jacobian = min_hex_jacobian;
            }
            found_valid = true;
        }
    }
    return found_valid ? global_min_jacobian : 0.0;
}


// 综合网格质量评估函数
// 返回一个综合得分，越高越好
float evaluate_overall_mesh_quality(TMesh* mesh, get_singularity_number<TMesh>& singularity_checker) {
    if (!mesh) return -std::numeric_limits<float>::infinity();

    // 1. 奇异线数量 (越少越好)
    singularity_checker.generate_singularity_number(mesh); // 确保获取最新的奇异线数量
    float singularity_score = -static_cast<float>(singularity_checker.singualarity_id) * 10.0f; // 每条奇异线扣10分

    // 2. 最小雅可比行列式 (越大越好，阈值0.1)
    double min_jacobian = calculate_global_min_jacobian(mesh);
    float jacobian_score = 0.0f;
    if (min_jacobian < 0.1) {
        jacobian_score = static_cast<float>(min_jacobian - 0.1) * 50.0f; // 低于0.1则大力惩罚
    } else {
        jacobian_score = static_cast<float>(min_jacobian) * 5.0f; // 高于0.1则给予正向奖励
    }
    
    // 3. （可选）平均雅可比行列式 (越大越好)
    // double avg_jacobian = calculate_average_min_jacobian(mesh);
    // float avg_jacobian_score = static_cast<float>(avg_jacobian) * 2.0f;

    // 综合得分，可以调整权重
    float overall_quality = singularity_score + jacobian_score; // + avg_jacobian_score;
    // std::cout << "  Quality Eval: S_Score=" << singularity_score << ", J_Score=" << jacobian_score << ", Overall=" << overall_quality << std::endl;
    return overall_quality;
}



// 通用值缩放辅助函数: 将 value 从 [source_min, source_max] 线性映射到 [target_min, target_max]
inline double scale_value(double value, double source_min, double source_max, double target_min, double target_max) {
    if (std::abs(source_max - source_min) < 1e-9) { // 处理 source_min 和 source_max 几乎相等的情况
        return target_min; // 或者 (target_min + target_max) / 2.0
    }
    double normalized_01 = (value - source_min) / (source_max - source_min);
    double scaled = target_min + normalized_01 * (target_max - target_min);
    return std::max(target_min, std::min(target_max, scaled)); // 裁剪到目标范围
}

// `sheet_energy` 的专属归一化函数
inline double normalize_sheet_energy(double energy_value) {
    if (std::abs(energy_value - 0.0) < 1e-9) { // 已坍缩
        return NORMALIZED_ENERGY_COLLAPSED;
    } else if (energy_value < (ENERGY_CANNOT_COLLAPSE_SENTINEL + 1.0) && energy_value > (ENERGY_CANNOT_COLLAPSE_SENTINEL -1.0) ) { // 不可坍缩, 比较浮点数时增加容差
        return NORMALIZED_ENERGY_CANNOT_COLLAPSE;
    } else if (energy_value > 0) { // 正能量，可坍缩
        double clipped_energy = std::max(ENERGY_MIN_POSITIVE_FOR_SCALING, std::min(energy_value, ENERGY_MAX_POSITIVE_FOR_SCALING));
        return scale_value(clipped_energy,
                           ENERGY_MIN_POSITIVE_FOR_SCALING, ENERGY_MAX_POSITIVE_FOR_SCALING,
                           NORMALIZED_ENERGY_POSITIVE_TARGET_MIN, NORMALIZED_ENERGY_POSITIVE_TARGET_MAX);
    } else {
        // 其他意外的负值 (非标记值)
        // std::cerr << "Warning: Unexpected negative energy value encountered during normalization: " << energy_value << std::endl;
        return NORMALIZED_ENERGY_CANNOT_COLLAPSE; // 暂时也映射为不可坍缩，或定义新的错误映射值
    }
}

#endif // STATE_FUNCTIONS_H
