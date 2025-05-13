### 训练配置文件 1

这个配置下训练的模型较为简单，能够实现同样终止结果下（all energy<0),但是mesh压缩比更好

{
  "training_loop_params": {
    "episodes_per_file": 50,
    "max_steps_per_episode": 10,
    "reward_threshold_to_stop": 100.0,
    "no_improvement_steps_threshold": 5,
    "improvement_tolerance": 1.0
  },
  "agent_params": {
    "eta": 0.5,
    "gamma": 0.99,
    "use_geometric_features": true,
    "replay_memory_capacity": 100000,
    "batch_size": 32,
    "target_update_freq": 10,
    "error_filter": {
      "initial_threshold": 0.05,
      "max_threshold": 0.95,
      "increase_rate_per_episode": 0.018
    },
    "epsilon_greedy": {
      "min_epsilon": 0.05,
      "decay_rate": 0.945,
      "start_epsilon": 0.9
    }
  },
  "optimizer_params": {
    "learning_rate": 0.001
  },
  "scheduler_params": {
    "step_size": 3,
    "gamma": 0.5
  },
  "network_params": {
    "hexmesh_qnet": {
      "node_feat_dim": 10,
      "hidden_dim": 64,
      "num_gnn_layers": 3
    },
    "ulinear_qnet": {
      "input_dim": 2,
      "output_dim": 1,
      "hidden_dim1": 10,
      "hidden_dim2": 10,
      "hidden_dim3": 10
    }
  },
  "normalization_params_comment": "这些参数目前在 state_functions.h 中硬编码。如果在此处更改，请同时更新头文件。",
  "normalization_params": {
    "energy_cannot_collapse_sentinel": -99999.0,
    "energy_min_positive_for_scaling": 0.001,
    "energy_max_positive_for_scaling": 500.0,
    "normalized_energy_collapsed": 0.0,
    "normalized_energy_cannot_collapse": -0.1,
    "normalized_energy_positive_target_min": 0.0,
    "normalized_energy_positive_target_max": 1.0,
    "other_features": [
        {"name": "sheet_on_boundary_ratio", "min_val": 0.0, "max_val": 1.0},
        {"name": "sheet_on_feature_ratio", "min_val": 0.0, "max_val": 1.0},
        {"name": "sheet_adjacent_feature_edges_count", "min_val": 0.0, "max_val": 15.0},
        {"name": "sheet_curvature_metric", "min_val": 0.0, "max_val": 1.0},
        {"name": "sheet_normal_variation", "min_val": 0.0, "max_val": 180.0},
        {"name": "sheet_dihedral_angle_deviation", "min_val": 0.0, "max_val": 180.0},
        {"name": "sheet_min_scaled_jacobian", "min_val": -1.0, "max_val": 1.0},
        {"name": "sheet_avg_scaled_jacobian", "min_val": -1.0, "max_val": 1.0},
        {"name": "sheet_distance_to_feature", "min_val": 0.0, "max_val": 100.0}
    ]
  },
  "reward_params_comment": "奖励系数在 state_functions.h 的 calc_reward 中硬编码。如果在此处更改，请同时更新头文件。",
  "reward_params": {
      "base_reward_decrease_singularity": 10.0,
      "base_reward_increase_singularity": -10.0,
      "base_reward_no_change_singularity": -1.0,
      "penalty_boundary_ratio_factor": 5.0,
      "penalty_feature_ratio_factor": 8.0,
      "reward_curvature_factor": 3.0,
      "penalty_normal_variation_factor": 0.15,
      "penalty_angle_deviation_factor": 0.15,
      "reward_low_min_jacobian_factor": 4.0,
      "penalty_high_min_jacobian_factor": 2.0,
      "penalty_distance_factor": 3.0,
      "reward_quality_change_factor": 15.0
  }
}
