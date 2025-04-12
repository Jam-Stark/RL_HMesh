# Intelligent Method for Hexahedral Mesh Topology Optimization Based on Sheet Operation with Reinforcement Learning

This project explores the integration of **sheet operations** and **deep reinforcement learning (RL)** for optimizing the topology of hexahedral meshes. The goal is to represent mesh optimization as a reinforcement learning problem by using the geometric and topological properties of the mesh as the environment's state, and applying RL to perform mesh modifications like insertions and collapses. This README provides an overview of the project, its methodology, key code components, and the workflow.

## Project Overview

The project aims to optimize the topology of hexahedral meshes using reinforcement learning methods, focusing on layer insertion and collapse operations. The primary objectives are:

* **Mesh Optimization** : Improve the quality of hexahedral mesh topology by reducing singularities and optimizing angle distributions.
* **Geometric Feature Preservation** : Ensure that mesh geometry, such as edge length consistency and smooth angles, is preserved during optimization.
* **Model Training Efficiency** : Enhance the training efficiency of the RL model while ensuring the model generalizes well to different mesh structures.
* **Reward Function Design** : Design a reward function based on mesh singularity distribution and angle features, ensuring the model optimizes according to real-world engineering needs.

## Goals

* **Enhance Mesh Quality** : Use RL for efficient optimization of the hexahedral mesh topology, significantly reducing singularities while ensuring mesh quality metrics (e.g., angle distribution, element shape) meet or exceed current optimization methods.
* **Preserve Geometric Features** : Ensure that mesh geometry (shape, boundary curvature) remains stable during the optimization process.
* **Improve Training and Generalization** : Develop an RL framework that can efficiently train and generalize across various mesh scenarios, adapting to complex topologies and geometries.

## Technical Metrics

* **Singularity Optimization** : Significantly reduce the number of singularities in the mesh.
* **Angle Distribution** : Increase the proportion of elements with angles close to 90°.
* **Geometric Feature Stability** : Control geometric errors (e.g., edge length errors, surface deviation) within a specified range (e.g., below 1%).
* **Training Efficiency** : Achieve faster convergence and reduced training time compared to traditional methods.
* **Model Generalization** : Ensure stable performance in various mesh scenarios with complex topologies.

## Components

### Python and C++ Integration with pybind11

* **Python Modules** : The neural network portion is implemented in Python, and pybind11 is used to interface between Python and C++.
* **C++ Mesh Operations** : Mesh topology operations (e.g., collapse, insert) are performed in C++.

### Key Code Directories and Files

1. **C++ Core Library for Hexahedral Mesh** :

* Location: `TopoHexMeshLib_sheet_singularity/HMeshLib/core` and `TopoHexMeshLib_sheet_singularity/HMeshLib/algorithm`.
* Contains core mesh elements, geometric calculations, and algorithms for mesh manipulation.

1. **Python Modules** :

* Location: `python_modules/`.
* The Python modules define the reinforcement learning agent (`agent.py`), reward function, and Q-network for decision-making.

1. **Mesh Topology Operations** :

* Insert and collapse operations are implemented in C++:
  * **Collapse Operation** : Implemented in `TopoHexMeshLib_sheet_singularity/src/sheet_operation.h` by a team member.
  * **Insert Operation** : A partial and unverified implementation is present in `TopoHexMeshLib_sheet_singularity/src/insert_mesh.h`.

1. **Main Program (C++)** :

* Location: `TopoHexMeshLib_sheet_singularity/src/main.cpp`.
* The main program handles the training loop, including interaction with Python via pybind11 and logging mechanisms.

1. **State and Reward Functions** :

* State calculation: `state_functions.h` defines the `State` class and methods for state calculation and action execution.
* Reward calculation: Based on singularity reduction, the reward function incentivizes the agent for reducing singularities in the mesh.

## Workflow and Execution

### Main Training Loop

1. **Initialization** :

* Logs are set up to track the training process.
* A Python interpreter is initialized, and the RL agent (`Agent`) is loaded.
* Mesh data is loaded, and mesh operations (insert/collapse) are prepared.

1. **Training Loop** :

* The outer loop runs for 50 episodes (configurable).
* Each episode processes a mesh file and optimizes it using RL.
* The inner loop processes steps within each episode:
  1. **State Calculation** : The current mesh state is calculated.
  2. **Action Selection** : The RL agent chooses an action based on the current state.
  3. **Action Execution** : The corresponding mesh operation (insert/collapse) is performed in C++.
  4. **Reward Calculation** : A reward is calculated based on mesh improvements (e.g., singularity reduction).
  5. **Experience Storage** : The state-action-reward transition is stored for learning.
  6. **Model Training** : The agent’s model is trained using experience replay.

1. **Model Saving and Logging** :

* After training, the RL model is saved.
* Logs provide insights into the training process, actions taken, and rewards received.

### Python Agent

* **Agent Class** : Manages RL operations, including choosing actions (`choose_action`), storing experiences (`remember`), and updating the Q-network (`replay`).
* **Q-Network** : Uses a fully connected neural network (`U_linear_QNet`) to approximate Q-values for each action.

### C++ Mesh Operations

* **State Calculation** : The state is represented as a list of sheets, where each sheet has an ID and energy.
* **Action Execution** : Actions modify the mesh by collapsing or inserting sheets, updating the mesh topology.
* **Energy Calculation** : Each sheet's collapse energy is calculated, and the reward is based on changes in mesh singularity.

### Key Files

1. **`main.cpp`** : The C++ main entry, responsible for initializing and running the training loop.
2. **`agent.py`** : The Python agent, defining the RL model and learning logic.
3. **`state_functions.h`** : Defines the `State` class and functions for interacting with mesh data.
4. **`sheet_operation.h`** : Contains the logic for performing mesh collapse operations.
5. **`insert_mesh.h`** : Contains the partial and unverified insert operation logic.

### Logs and Debugging

* Logs are stored in `logs/` with timestamps to track each session.
* Key logs include state transitions, action choices, rewards, and mesh changes.

## Future Improvements

1. **Fix Insert Operation** : The insert operation algorithm needs further validation and refinement.
2. **Enhance Reward Function** : Fine-tuning the reward function to better align with mesh quality metrics.
3. **Optimize Model Architecture** : Investigate the possibility of using more advanced RL algorithms for improved optimization.

## Conclusion

This project presents an innovative approach to mesh optimization using deep reinforcement learning. By combining mesh topology manipulation with RL, we aim to develop a robust solution for mesh optimization that balances efficiency and quality. The integration of Python and C++ via pybind11 allows for seamless interaction between the two languages, enabling powerful mesh manipulations while leveraging the decision-making capabilities of RL.
