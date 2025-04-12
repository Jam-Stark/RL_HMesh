/*#include "sheet_operation.h"
#include "topoMesh.h"
#include <unordered_set>
#include <iostream>

namespace HMeshLib {

template<typename M>
void sheet_operation<M>::insert_sheet(E* target_edge) {
    // Step 1: 获取目标边相关的面
    std::vector<F*> sheet_faces = get_one_sheet(target_edge);
    if (sheet_faces.empty()) {
        std::cerr << "Error: No faces found for the target sheet." << std::endl;
        return;
    }

    // Step 2: 标记边界和特征线
    std::unordered_set<E*> boundary_edges;
    for (auto& face : sheet_faces) {
        for (auto& edge : face->edges()) {
            if (edge->is_boundary()) {
                boundary_edges.insert(edge);
            }
        }
    }

    // Step 3: 创建新的顶点
    std::vector<V*> new_vertices;
    for (auto& face : sheet_faces) {
        for (auto& vertex : face->vertices()) {
            V* new_vertex = construct_new_vertex();
            new_vertex->set_position(vertex->get_position());  // 复制位置
            new_vertices.push_back(new_vertex);
        }
    }

    // Step 4: 创建新的面
    std::vector<F*> new_faces;
    for (auto& face : sheet_faces) {
        F* new_face = constrcut_new_face(new_vertices);  // 传入新的顶点
        new_faces.push_back(new_face);
    }

    // Step 5: 创建新的六面体
    std::vector<H*> new_hexes;
    for (auto& face : sheet_faces) {
        H* new_hex = new H();
        new_hex->set_faces(new_faces);  // 将新面连接到六面体
        new_hexes.push_back(new_hex);
    }

    // Step 6: 更新邻接关系
    for (auto& new_face : new_faces) {
        for (auto& edge : new_face->edges()) {
            edge->add_face(new_face);  // 将新面添加到边的邻接关系中
        }
    }

    // Step 7: 重新计算法向量
    mesh->recompute_normals();

    // Step 8: 检查并标记奇异点
    mesh->mark_singularity();
    for (auto& new_hex : new_hexes) {
        new_hex->check_singularity();  // 检查并标记新六面体中的奇异点
    }

    // Step 9: 更新网格的边界条件
    mesh->update_boundary_conditions();

    std::cout << "Sheet insertion completed successfully." << std::endl;
}

template<typename M>
void sheet_operation<M>::update_boundary_conditions() {
    // 更新所有边的边界条件
    for (auto& edge : es) {
        if (edge->is_boundary()) {
            edge->update_boundary();  // 更新边的边界条件
        }
    }

    // 更新所有顶点的边界条件
    for (auto& vertex : vs) {
        if (vertex->is_boundary()) {
            vertex->update_boundary();  // 更新顶点的边界条件
        }
    }

    // 更新所有面的边界条件
    for (auto& face : fs) {
        if (face->is_boundary()) {
            face->update_boundary();  // 更新面的边界条件
        }
    }

    // 输出边界更新的状态
    std::cout << "Boundary conditions updated." << std::endl;
}

template<typename M>
void sheet_operation<M>::recompute_normals() {
    // 遍历所有面，重新计算法向量
    for (auto& face : fs) {
        face->compute_normal();  // 使用面的顶点计算法向量
    }

    // 输出法向量计算状态
    std::cout << "Normals recomputed for all faces." << std::endl;
}

template<typename M>
void sheet_operation<M>::mark_singularity() {
    // 遍历网格中的所有顶点，检查是否为奇异性
    for (auto& vertex : vs) {
        if (vertex->degree() > 4) {  // 假设度数大于4的节点为奇异节点
            vertex->mark_singular();  // 标记为奇异节点
        }
    }

    // 遍历网格中的所有边，检查是否为奇异性
    for (auto& edge : es) {
        if (edge->degree() > 2) {  // 假设度数大于2的边为奇异边
            edge->mark_singular();  // 标记为奇异边
        }
    }

    // 遍历网格中的所有面，检查是否为奇异性
    for (auto& face : fs) {
        if (face->degree() > 4) {  // 假设度数大于4的面为奇异面
            face->mark_singular();  // 标记为奇异面
        }
    }

    std::cout << "Singularities marked for vertices, edges, and faces." << std::endl;
}

}  // namespace HMeshLib

*/