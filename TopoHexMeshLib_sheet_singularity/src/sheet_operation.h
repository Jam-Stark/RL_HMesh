#pragma once
#include "topoMesh.h"
#include "TopoHexMeshLib_sheet_singularity\HMeshLib\core\Topology\TopoM.h"
#include <fstream>
#include <cstdlib>
#include <chrono>
#include <ctime>
#include <sstream>
#include <iomanip>

namespace HMeshLib
{
    // 新增日志捕捉代码，自动写入 F:/RL_HMesh/logs/<build_mode>/<session_id>/sheet_operation.log
    inline std::ofstream& getSheetLog() {
        static std::ofstream ofs; // 将 ofstream 声明为 static
        if (!ofs.is_open()) {
            const char* session_env = std::getenv("RL_HMESH_SESSION_ID");
            const char* build_mode_env = std::getenv("RL_HMESH_BUILD_MODE");
            
            std::string session = session_env ? session_env : "default";
            std::string build_mode = build_mode_env ? build_mode_env : "Release";
            
            std::string logPath = "F:/RL_HMesh/logs/" + build_mode + "/" + session + "/sheet_operation.log";
            ofs.open(logPath, std::ios::app);
        }
        return ofs;
    }

    inline void log(const std::string &msg) {
        std::ofstream &ofs = getSheetLog();
        auto now = std::chrono::system_clock::now();
        std::time_t now_time = std::chrono::system_clock::to_time_t(now);
        std::tm* tm_ptr = std::localtime(&now_time);
        char buffer[10];
        std::strftime(buffer, sizeof(buffer), "[%H:%M:%S]", tm_ptr);
        ofs << buffer << " " << msg << std::endl;
    }
#define MPI 3.141592653589793238

	class MeshFeatures
	{
	};

	template<typename M>
	class sheet_operation
	{
		typedef typename M::V V;
		typedef typename M::E E;
		typedef typename M::F F;
		typedef typename M::H H;
	public:
		sheet_operation(M* input_mesh);
		sheet_operation();
		~sheet_operation();

		V* construct_new_vertex();
		E* construct_new_edge(std::vector<V*> vs);
		F* constrcut_new_face(std::vector<V*> vs);
		void vertex_position_policy(V* v1, V* v2, V* target);
		void edge_sharp_policy(E* e1, E* e2, E* target);
		void edge_total_angle_policy_collapse(E* e1, E* e2, E* target);
		void mark_elements_attribute_collapse(std::vector<F*> fs);

        // void reset_sheet_attributes() {
        //     for (M::MEIterator eite(mesh); !eite.end(); eite++) {
        //         E* be = *eite;
        //         be->sheet() = 0;
        //     }
        // }
		
		std::vector<E*> get_one_sheet(E* e);
		int get_mesh_sheet_number();
		void collapse_one_sheet2(std::vector<E*> sheet);

		double vector_angle(CPoint a, CPoint b);
		void edge_angle(E* e);
		void edge_ideal_degree(E* e);
		void compute_edge_energy();
		std::vector<std::vector<E*>> get_sheet_parallel_edges(std::vector<E*>sheet);
		double predict_sheet_collapse_energy(std::vector<E*> sheet);
		
		/* get new info from each sheet*/
		double get_sheet_on_boundary_ratio(std::vector<E*> sheet);
		double get_sheet_on_feature_ratio(std::vector<E*> sheet);
		
		bool sheet_endpoints_are_corners_or_boundary(std::vector<E*> sheet);
		int get_sheet_adjacent_feature_edges_count(std::vector<E*> sheet);

		//通过计算sheet两端点之间的欧氏距离与sheet的总路径长度的比值来度量曲率
		double get_sheet_curvature_metric(std::vector<E*> sheet);

		// 计算sheet中的边相邻面之间的法线角度变化的平均值
		double get_sheet_normal_variation(std::vector<E*> sheet);

		// 计算sheet中边的实际二面角与理想二面角（基于边的属性）之间的平均偏差
		double get_sheet_dihedral_angle_deviation(std::vector<E*> sheet);

		// 添加计算sheet邻近六面体缩放雅可比指标的函数声明
		std::pair<double, double> get_sheet_adjacent_hex_jacobian(std::vector<E*> sheet);

		// 计算mesh的各项属性(点 边 面数量 & 几何 & 拓扑各项信息),方便对比算法处理后的mesh和原mesh的对比,检验算法的效果
		MeshFeatures get_mesh_features();

		int get_sheet_valence_along_path(const std::vector<E*>& sheet);

        double get_distance_to_boundary(const std::vector<E*>& sheet);
		double get_distance_to_feature(const std::vector<E*>& sheet);
		
	private:
		M* mesh;
		std::string filename;
	};

	template<typename M>
	sheet_operation<M>::sheet_operation(M* input_mesh)
	{
		mesh = input_mesh;

	}

	template<typename M>
	sheet_operation<M>::sheet_operation()
	{

	}

	template<typename M>
	sheet_operation<M>::~sheet_operation()
	{

	}

	/*collaspe*/
	template<typename M>
	std::vector<typename sheet_operation<M>::E*> sheet_operation<M>::get_one_sheet(E* e)
	{
		std::vector<E*> sheet;
		std::queue<E*> eQueue;
		std::unordered_set<int> visited_edge_ids; // 使用局部的 set 跟踪访问过的边 ID

		if (!e) return sheet; // 处理空指针输入

		eQueue.push(e);
		visited_edge_ids.insert(e->id());

		while (!eQueue.empty())
		{
			E* be = eQueue.front();
			sheet.push_back(be);
			eQueue.pop();

			if (!be) continue; // 安全检查

			for (int ehIndex = 0; ehIndex < be->neighbor_hs.size(); ehIndex++)
			{
				H* h = mesh->idHexs(be->neighbor_hs[ehIndex]);
				if (!h) continue; // 安全检查

				std::vector<E*> parallel_es = mesh->e_parallel_e_in_hex(h, be);
				for (int peIndex = 0; peIndex < parallel_es.size(); peIndex++)
				{
					E* pe = parallel_es[peIndex];
					if (!pe) continue; // 安全检查

					// 检查是否已在当前搜索中访问过
					if (visited_edge_ids.find(pe->id()) == visited_edge_ids.end())
					{
						eQueue.push(pe);
						visited_edge_ids.insert(pe->id());
					}
				}
			}
		}
		// 不需要重置 mark()，因为它没有被修改
		return sheet;
	}

	template<typename M>
	int sheet_operation<M>::get_mesh_sheet_number()
	{
		int sheet_id = 0;
		for (M::MEIterator eite(mesh); !eite.end(); eite++)
		{
			E* be = *eite;
			if (be->sheet()!=0) continue;

			std::vector<E*>sheet = get_one_sheet(be);
			//mark sheet
			sheet_id++;
			for (int sheetIndex = 0; sheetIndex < sheet.size(); sheetIndex++)
			{
				E* sheete = sheet[sheetIndex];
				sheete->sheet() = sheet_id;
			}
		}
		return sheet_id;
	}

	template<typename M>
	typename sheet_operation<M>::V* sheet_operation<M>::construct_new_vertex()
	{
		V* newv = new V();
		mesh->maxVid()++;
		newv->id() = mesh->maxVid();
		mesh->vs.push_back(newv);
		mesh->m_map_vertices.insert(std::pair<int, V*>(newv->id(), newv));
		return newv;
	}

	template<typename M>
	typename sheet_operation<M>::E* sheet_operation<M>::construct_new_edge(std::vector<V*> vs)
	{
		E* newe = new E();
		mesh->maxEid()++;
		newe->id() = mesh->maxEid();
		mesh->es.push_back(newe);
		mesh->m_map_edges.insert(std::pair<int, E*>(newe->id(), newe));
		for (int i = 0; i < 2; i++)
		{
			newe->vs.push_back(vs[i]->id());
			vs[i]->neighbor_es.push_back(newe->id());
		}
		return newe;
	}

	template<typename M>
	typename sheet_operation<M>::F* sheet_operation<M>::constrcut_new_face(std::vector<V*> vs)
	{
		F* newf = new F();
		mesh->maxFid()++;
		newf->id() = mesh->maxFid();
		mesh->fs.push_back(newf);
		mesh->m_map_faces.insert(std::pair<int, F*>(newf->id(), newf));

		for (int i = 0; i < 4; i++)
		{
			newf->vs.push_back(vs[i]->id());
			vs[i]->push_back_neighbor_f(newf->id());

			E* fe = mesh->VerticesEdge(vs[i], vs[(i + 1) % 4]);
			assert(fe != NULL);
			newf->es.push_back(fe->id());
			fe->neighbor_fs.push_back(newf->id());
		}
		return newf;
	}

	template<typename M>
	void sheet_operation<M>::vertex_position_policy(V* v1, V* v2, V* target)
	{
		if (v1->feature_vertex() < 0 && v2->feature_vertex() >= 0)
		{
			target->position() = v1->position();
			target->feature_vertex() = v1->feature_vertex();
			target->corner() = v1->corner();
		}
		else if (v1->feature_vertex() >= 0 && v2->feature_vertex() < 0)
		{
			target->position() = v2->position();
			target->feature_vertex() = v2->feature_vertex();
			target->corner() = v2->corner();
		}
		else if (v1->feature_vertex() && !v2->feature_vertex())
		{
			target->position() = v1->position();
			target->feature_vertex() = v1->feature_vertex();
		}
		else if (v2->feature_vertex() && !v1->feature_vertex())
		{
			target->position() = v2->position();
			target->feature_vertex() = v2->feature_vertex();
		}
		else if (v1->boundary() && !v2->boundary())
		{
			target->position() = v1->position();
		}
		else if (!v1->boundary() && v2->boundary())
		{
			target->position() = v2->position();
		}
		else
		{
			target->position() = (v1->position() + v2->position()) * 0.5;
			target->feature_vertex() = v2->feature_vertex();
		}
	}

	template<typename M>
	void sheet_operation<M>::edge_sharp_policy(E* e1, E* e2, E* target)
	{
		if (e1->sharp() && !e2->sharp())
		{
			target->sharp() = e1->sharp();
		}
		else if (!e1->sharp() && e2->sharp())
		{
			target->sharp() = e2->sharp();
		}

	}

	template<typename M>
	void sheet_operation<M>::edge_total_angle_policy_collapse(E* e1, E* e2, E* target)
	{
		double predict_angle = 0;
		if (e1->boundary() && e2->boundary())
		{
			if (e1->sharp() && !e2->sharp())
			{
				predict_angle = e1->total_angle();
			}
			else if (e2->sharp() && !e1->sharp())
			{
				predict_angle = e2->total_angle();
			}
			else
			{
				predict_angle = (e1->total_angle() + e2->total_angle()) * 0.5;
			}
		}
		else if (e1->boundary() && !e2->boundary())
		{
			predict_angle = e1->total_angle();
		}
		else if (!e1->boundary() && e2->boundary())
		{
			predict_angle = e2->total_angle();
		}
		else
		{
			predict_angle = (e1->total_angle() + e2->total_angle()) * 0.5;
		}
		target->total_angle() = predict_angle;
	}

	template<typename M>
	void sheet_operation<M>::mark_elements_attribute_collapse(std::vector<F*> fs)
	{
		try {
			log("Starting to mark element attributes, processing faces count: " + std::to_string(fs.size()));
	
			// 检查输入参数
			if (fs.empty()) {
				// log("Warning: Input face list is empty");
				return;
			}
	
			// 第一次遍历：标记边界并修正法线
			for (int fIndex = 0; fIndex < fs.size(); fIndex++)
			{
				try {
					F* f = fs[fIndex];
	
					// 检查面是否有效
					if (f == nullptr) {
						// log("Error: Detected null face pointer, index: " + std::to_string(fIndex));
						continue;
					}
	
					// log("Processing face ID: " + std::to_string(f->id()));
	
					// if (f->id() == 465) {
					// 	log("Special debug for face 465:");
					// 	log("  neighbor_hs count: " + std::to_string(f->neighbor_hs.size()));
					// 	log("  vertices count: " + std::to_string(f->vs.size()));
					// 	log("  edges count: " + std::to_string(f->es.size()));
					// }
	
					// 检查邻接六面体列表是否有效
					if (f->neighbor_hs.empty()) {
						// log("Warning: Face " + std::to_string(f->id()) + " has no adjacent hexahedra");
					}
	
					//mark boundary
					if (f->neighbor_hs.size() == 1)
					{
						f->boundary() = true;
						// log("Face " + std::to_string(f->id()) + " marked as boundary face");
					}
					else
					{
						f->boundary() = false;
					}
	
					// 修正六面体面法线
					for (int fhIndex = 0; fhIndex < f->neighbor_hs.size(); fhIndex++)
					{
						try {
							// 检查六面体ID是否有效
							if (f->neighbor_hs[fhIndex] < 0 || f->neighbor_hs[fhIndex] >= mesh->maxHid()) {
								// log("Error: Face " + std::to_string(f->id()) + " has invalid adjacent hexahedron ID: " + std::to_string(f->neighbor_hs[fhIndex]) + " with max Hid: "+ std::to_string(mesh->maxHid()));
								continue;
							}
	
							// log("About to get hexahedron with ID: " + std::to_string(f->neighbor_hs[fhIndex]));
							H* h = mesh->idHexs(f->neighbor_hs[fhIndex]);
							// log("Successfully got hexahedron with ID: " + std::to_string(h->id()));
	
							// 检查六面体是否存在
							if (h == nullptr) {
								// log("Error: Unable to get adjacent hexahedron ID: " + std::to_string(f->neighbor_hs[fhIndex]) + " for face " + std::to_string(f->id()));
								continue;
							}
	
							try {
								mesh->revise_hex_face_normal(h);
							}
							catch (const std::exception& e) {
								// log("Error while correcting hexahedron face normal, hex ID: " + std::to_string(h->id()) + ", exception: " + std::string(e.what()));
							}
							catch (...) {
								// log("Unknown error while correcting hexahedron face normal, hex ID: " + std::to_string(h->id()));
							}
						}
						catch (const std::exception& e) {
							// log("Exception in hexahedron processing: " + std::string(e.what()) + ", face ID: " + std::to_string(f->id()) + ", hex index: " + std::to_string(fhIndex));
						}
						catch (...) {
							// log("Unknown exception in hexahedron processing, face ID: " + std::to_string(f->id()) + ", hex index: " + std::to_string(fhIndex));
						}
					}
				}
				catch (const std::exception& e) {
					// log("Exception in first traversal at face index " + std::to_string(fIndex) + ": " + std::string(e.what()));
				}
				catch (...) {
					// log("Unknown exception in first traversal at face index " + std::to_string(fIndex));
				}
			}
	
			// log("First traversal completed, now marking vertices and edges");
	
			// 第二次遍历：标记顶点和边
			for (int fIndex = 0; fIndex < fs.size(); fIndex++)
			{
				try {
					F* f = fs[fIndex];
	
					// 检查面是否有效
					if (f == nullptr) {
						// log("Skipping null face, index: " + std::to_string(fIndex));
						continue;
					}
	
					// log("Processing face for vertices/edges ID: " + std::to_string(f->id()));
	
					// 标记面顶点
					try {
						// log("Marking " + std::to_string(f->vs.size()) + " vertices for face " + std::to_string(f->id()));
						for (int fvIndex = 0; fvIndex < f->vs.size(); fvIndex++)
						{
							try {
								// 检查顶点ID是否有效
								if (f->vs[fvIndex] < 0 || f->vs[fvIndex] >= mesh->maxVid()) {
									// log("Error: Face " + std::to_string(f->id()) + " has invalid vertex ID: " + std::to_string(f->vs[fvIndex]));
									continue;
								}
	
								V* fv = mesh->idVertices(f->vs[fvIndex]);
	
								// 检查顶点是否存在
								if (fv == nullptr) {
									// log("Error: Unable to get vertex ID: " + std::to_string(f->vs[fvIndex]) + " for face " + std::to_string(f->id()));
									continue;
								}
	
								if (f->boundary())
								{
									fv->boundary() = true;
									// log("Vertex " + std::to_string(fv->id()) + " marked as boundary vertex");
								}
								else
								{
									// 检查顶点的所有相邻面是否有边界面
									bool is_boundary = false;
									for (int vfIndex = 0; vfIndex < fv->neighbor_fs.size(); vfIndex++)
									{
										try {
											// 检查面ID是否有效
											if (fv->neighbor_fs[vfIndex] < 0 || fv->neighbor_fs[vfIndex] >= mesh->maxFid()) {
												// log("Warning: Vertex " + std::to_string(fv->id()) + " has invalid adjacent face ID: " + std::to_string(fv->neighbor_fs[vfIndex]));
												continue;
											}
	
											F* fvf = mesh->idFaces(fv->neighbor_fs[vfIndex]);
	
											// 检查面是否存在
											if (fvf == nullptr) {
												// log("Warning: Unable to get adjacent face ID: " + std::to_string(fv->neighbor_fs[vfIndex]) + " for vertex " + std::to_string(fv->id()));
												continue;
											}
	
											if (fvf->boundary())
											{
												is_boundary = true;
												// log("Vertex " + std::to_string(fv->id()) + " marked as boundary vertex (via adjacent face)");
												break;
											}
										}
										catch (const std::exception& e) {
											// log("Exception in vertex-face processing: " + std::string(e.what()) + ", vertex ID: " + std::to_string(fv->id()) + ", face index: " + std::to_string(vfIndex));
										}
										catch (...) {
											// log("Unknown exception in vertex-face processing, vertex ID: " + std::to_string(fv->id()) + ", face index: " + std::to_string(vfIndex));
										}
									}
									fv->boundary() = is_boundary;
								}
							}
							catch (const std::exception& e) {
								// log("Exception processing vertex at index " + std::to_string(fvIndex) + " for face " + std::to_string(f->id()) + ": " + std::string(e.what()));
							}
							catch (...) {
								// log("Unknown exception processing vertex at index " + std::to_string(fvIndex) + " for face " + std::to_string(f->id()));
							}
						}
					}
					catch (const std::exception& e) {
						// log("Exception in vertices processing for face " + std::to_string(f->id()) + ": " + std::string(e.what()));
					}
					catch (...) {
						// log("Unknown exception in vertices processing for face " + std::to_string(f->id()));
					}
	
					// 标记面边
					try {
						// log("Marking " + std::to_string(f->es.size()) + " edges for face " + std::to_string(f->id()));
						for (int feIndex = 0; feIndex < f->es.size(); feIndex++)
						{
							try {
								// 检查边ID是否有效
								if (f->es[feIndex] < 0 || f->es[feIndex] >= mesh->maxEid()) {
									// log("Error: Face " + std::to_string(f->id()) + " has invalid edge ID: " + std::to_string(f->es[feIndex]));
									continue;
								}
	
								E* fe = mesh->idEdges(f->es[feIndex]);
	
								// 检查边是否存在
								if (fe == nullptr) {
									// log("Error: Unable to get edge ID: " + std::to_string(f->es[feIndex]) + " for face " + std::to_string(f->id()));
									continue;
								}
	
								// 边界标记
								if (f->boundary())
								{
									fe->boundary() = true;
									// log("Edge " + std::to_string(fe->id()) + " marked as boundary edge");
								}
								else
								{
									// 检查是否有任何相邻面是边界
									bool is_boundary = false;
									for (int efIndex = 0; efIndex < fe->neighbor_fs.size(); efIndex++)
									{
										try {
											// 检查面ID是否有效
											if (fe->neighbor_fs[efIndex] < 0 || fe->neighbor_fs[efIndex] >= mesh->maxFid()) {
												// log("Warning: Edge " + std::to_string(fe->id()) + " has invalid adjacent face ID: " + std::to_string(fe->neighbor_fs[efIndex]));
												continue;
											}
	
											F* fef = mesh->idFaces(fe->neighbor_fs[efIndex]);
	
											// 检查面是否存在
											if (fef == nullptr) {
												// log("Warning: Unable to get adjacent face ID: " + std::to_string(fe->neighbor_fs[efIndex]) + " for edge " + std::to_string(fe->id()));
												continue;
											}
	
											if (fef->boundary())
											{
												is_boundary = true;
												// log("Edge " + std::to_string(fe->id()) + " marked as boundary edge (via adjacent face)");
												break;
											}
										}
										catch (const std::exception& e) {
											// log("Exception in edge-face processing: " + std::string(e.what()) + ", edge ID: " + std::to_string(fe->id()) + ", face index: " + std::to_string(efIndex));
										}
										catch (...) {
											// log("Unknown exception in edge-face processing, edge ID: " + std::to_string(fe->id()) + ", face index: " + std::to_string(efIndex));
										}
									}
									fe->boundary() = is_boundary;
								}
	
								try {
									// 奇异性标记
									if (fe->boundary())
									{
										if (fe->neighbor_hs.size() != 2)
										{
											fe->singularity() = true;
											// log("Boundary edge " + std::to_string(fe->id()) + " marked as singular edge, adjacent hexahedra count: " + std::to_string(fe->neighbor_hs.size()));
										}
										else
										{
											fe->singularity() = false;
										}
									}
									else
									{
										if (fe->neighbor_hs.size() != 4)
										{
											fe->singularity() = true;
											// log("Non-boundary edge " + std::to_string(fe->id()) + " marked as singular edge, adjacent hexahedra count: " + std::to_string(fe->neighbor_hs.size()));
										}
										else
										{
											fe->singularity() = false;
										}
									}
								}
								catch (const std::exception& e) {
									// log("Exception in edge singularity marking: " + std::string(e.what()) + ", edge ID: " + std::to_string(fe->id()));
								}
								catch (...) {
									// log("Unknown exception in edge singularity marking, edge ID: " + std::to_string(fe->id()));
								}
	
								try {
									// 简化能量计算
									if (fe->boundary())
									{
										if (fe->sharp())
										{
											// 检查总角度是否合理
											if (fe->total_angle() < 0 || fe->total_angle() > 360) {
												// log("Warning: Edge " + std::to_string(fe->id()) + " has abnormal total angle value: " + std::to_string(fe->total_angle()));
												// 使用默认值避免计算错误
												fe->ideal_degree() = 2;
											}
											else {
												int ideal_degree = round(fe->total_angle() / 90.0);
												ideal_degree = ideal_degree == 0 ? 1 : ideal_degree;
												fe->ideal_degree() = ideal_degree;
												// log("Sharp boundary edge " + std::to_string(fe->id()) + " ideal degree set to: " + std::to_string(ideal_degree) + " (total angle: " + std::to_string(fe->total_angle()) + ")");
											}
										}
										else
										{
											fe->ideal_degree() = 2;
											// log("Regular boundary edge " + std::to_string(fe->id()) + " ideal degree set to: 2");
										}
									}
									else
									{
										fe->ideal_degree() = 4;
										// log("Non-boundary edge " + std::to_string(fe->id()) + " ideal degree set to: 4");
									}
	
									// 计算简化能量
									fe->sim_energy() = (int)fe->ideal_degree() - (int)fe->neighbor_hs.size();
									// log("Edge " + std::to_string(fe->id()) + " simplification energy: " + std::to_string(fe->sim_energy()) +
									// 	" (ideal degree: " + std::to_string(fe->ideal_degree()) + ", actual adjacent hexahedra: " + std::to_string(fe->neighbor_hs.size()) + ")");
								}
								catch (const std::exception& e) {
									// log("Exception in edge energy calculation: " + std::string(e.what()) + ", edge ID: " + std::to_string(fe->id()));
								}
								catch (...) {
									// log("Unknown exception in edge energy calculation, edge ID: " + std::to_string(fe->id()));
								}
							}
							catch (const std::exception& e) {
								// log("Exception processing edge at index " + std::to_string(feIndex) + " for face " + std::to_string(f->id()) + ": " + std::string(e.what()));
							}
							catch (...) {
								// log("Unknown exception processing edge at index " + std::to_string(feIndex) + " for face " + std::to_string(f->id()));
							}
						}
					}
					catch (const std::exception& e) {
						// log("Exception in edges processing for face " + std::to_string(f->id()) + ": " + std::string(e.what()));
					}
					catch (...) {
						// log("Unknown exception in edges processing for face " + std::to_string(f->id()));
					}
				}
				catch (const std::exception& e) {
					// log("Exception in second traversal at face index " + std::to_string(fIndex) + ": " + std::string(e.what()));
				}
				catch (...) {
					// log("Unknown exception in second traversal at face index " + std::to_string(fIndex));
				}
			}
	
			log("Element attribute marking completed");
		}
		catch (const std::exception& e) {
			log("Major exception in mark_elements_attribute_collapse: " + std::string(e.what()));
		}
		catch (...) {
			log("Unknown major exception in mark_elements_attribute_collapse");
		}
	}

	template<typename M>
	void sheet_operation<M>::collapse_one_sheet2(std::vector<E*> sheet)
	{
		log("Starting to process sheet, size: " + std::to_string(sheet.size()));
		/*get all the hexs and the include elements*/
		std::set<H*> delete_hs;//get all hex in sheet
		std::set<F*> delete_fs;//delete faces
		std::set<E*> delete_es;//delete edges
		std::set<V*> delete_vs;//delete vs
	
		// 先收集需要删除的六面体及其关联的元素
		for (int eIndex = 0; eIndex < sheet.size(); eIndex++)
		{
			E* e = sheet[eIndex];
			if (!e) continue; // Sanity check
			e->sheetId() = 1;
			// log("Processing edge id: " + std::to_string(e->id())); // Already Commented out
			//find all hex
			for (int ehIndex = 0; ehIndex < e->neighbor_hs.size(); ehIndex++)
			{
				H* h = mesh->idHexs(e->neighbor_hs[ehIndex]);
				if (!h || h->is_delete()) continue; // Check for null and already marked
	
				h->is_delete() = true;
				h->sheetId() = -1; // Mark hex related to the sheet collapse
				delete_hs.insert(h);
	
				/*mark the delete face and the elements of it*/
				for (int hvIndex = 0; hvIndex < h->vs.size(); hvIndex++)
				{
					V* hv = mesh->idVertices(h->vs[hvIndex]);
					if (hv && !hv->is_delete()) { // Check for null and already marked
						hv->is_delete() = true;
						delete_vs.insert(hv);
					}
				}
				for (int heIndex = 0; heIndex < h->es.size(); heIndex++)
				{
					E* he = mesh->idEdges(h->es[heIndex]);
					if (he && !he->is_delete()) { // Check for null and already marked
						he->is_delete() = true;
						delete_es.insert(he);
					}
				}
				for (int hfIndex = 0; hfIndex < h->fs.size(); hfIndex++)
				{
					F* hf = mesh->idFaces(h->fs[hfIndex]);
					if (hf && !hf->is_delete()) { // Check for null and already marked
						hf->is_delete() = true;
						delete_fs.insert(hf);
					}
				}
			}
		}
		log("Collected elements to delete - Hexes: " + std::to_string(delete_hs.size()) +
			", Faces: " + std::to_string(delete_fs.size()) +
			", Edges: " + std::to_string(delete_es.size()) +
			", Vertices: " + std::to_string(delete_vs.size()));
	
		/*exit if there is self-intersection*/
		for (std::set<H*>::iterator hite = delete_hs.begin(); hite != delete_hs.end(); hite++)
		{
			H* h = *hite;
			if (!h) continue;
			int count = 0;
			for (int heIndex = 0; heIndex < h->es.size(); heIndex++)
			{
				E* e = mesh->idEdges(h->es[heIndex]);
				if (e && e->sheetId()) // Check for null
				{
					count++;
				}
			}
			if (count > 4)
			{
				log("Self-intersection detected at hex id: " + std::to_string(h->id()) + ". Aborting collapse.");
				std::cout << "The current sheet has self-intersection" << std::endl;
	
				// 清除之前设置的is_delete标记，因为我们不会继续折叠
				for (H* h_del : delete_hs) { if (h_del) h_del->is_delete() = false; }
				for (F* f_del : delete_fs) { if (f_del) f_del->is_delete() = false; }
				for (E* e_del : delete_es) { if (e_del) e_del->is_delete() = false; }
				for (V* v_del : delete_vs) { if (v_del) v_del->is_delete() = false; }
				// Also reset sheetId for the original sheet edges
				for (E* e_sheet : sheet) { if (e_sheet) e_sheet->sheetId() = 0; }
				// Reset sheetId for hexes that were marked
				for (H* h_del : delete_hs) { if (h_del) h_del->sheetId() = 0; }
	
				return;
			}
		}
	
		/*new elements*/
		std::vector<V*> newvs;
		std::vector<E*> newes;
		std::vector<F*> newfs;
		/*add new vertices*/
		log("Starting new vertex creation and topology update");
		for (int eIndex = 0; eIndex < sheet.size(); eIndex++)
		{
			E* e = sheet[eIndex];
			if (!e) continue;
			V* newv = construct_new_vertex();
			// log("Created new vertex id: " + std::to_string(newv->id())); // Already Commented out
			newvs.push_back(newv);
	
			V* ev1 = mesh->idVertices(e->vs[0]);
			V* ev2 = mesh->idVertices(e->vs[1]);
			if (!ev1 || !ev2) {
				log("Error: Could not find vertices for edge " + std::to_string(e->id()) + ". Skipping topology update for this edge.");
				continue; // Skip if original vertices are missing
			}
			vertex_position_policy(ev1, ev2, newv);
	
			std::vector<V*> evs = { ev1,ev2 };
			e->setNewv(newv); // Associate new vertex with the original sheet edge
			//change relations
			for (int evIndex = 0; evIndex < 2; evIndex++)
			{
				V* oldv = evs[evIndex];
				if (!oldv) continue;
	
				//change edge relations - update edges connected to oldv (if not deleted) to point to newv
				std::vector<int> old_neighbor_es = oldv->neighbor_es; // Copy because original might be modified indirectly
				for (int eid : old_neighbor_es)
				{
					E* ve = mesh->idEdges(eid);
					// Only update edges that are NOT marked for deletion
					if (!ve || ve->is_delete()) continue;
					int vevIndex = ve->vertexIndex(oldv->id());
					if (vevIndex != -1) {
						ve->vs[vevIndex] = newv->id();
						newv->neighbor_es.push_back(ve->id());
						// Remove oldv's reference to this edge? No, oldv is being deleted.
					}
				}
				//change face relations - update faces connected to oldv (if not deleted) to point to newv
				std::vector<int> old_neighbor_fs = oldv->neighbor_fs; // Copy
				for (int fid : old_neighbor_fs)
				{
					F* vf = mesh->idFaces(fid);
					// Only update faces that are NOT marked for deletion
					if (!vf || vf->is_delete()) continue;
					int vfvIndex = vf->vertexIndex(oldv->id());
					if (vfvIndex != -1)
					{
						vf->vs[vfvIndex] = newv->id();
						newv->neighbor_fs.push_back(vf->id());
					}
				}
				//change hex relations - update hexes connected to oldv (if not deleted) to point to newv
				std::vector<int> old_neighbor_hs = oldv->neighbor_hs; // Copy
				for (int hid : old_neighbor_hs)
				{
					H* vh = mesh->idHexs(hid);
					// Only update hexes that are NOT marked for deletion
					if (!vh || vh->is_delete()) continue;
					int vhvIndex = vh->vertexIndex(oldv->id());
					if (vhvIndex != -1)
					{
						vh->vs[vhvIndex] = newv->id();
						newv->neighbor_hs.push_back(vh->id());
					}
				}
			}
		}
		log("Vertex creation and topology update completed. New vertices: " + std::to_string(newvs.size()));
	
		/*create new edges and faces*/
		log("Starting new edge/face creation");
		for (std::set<H*>::iterator hite = delete_hs.begin(); hite != delete_hs.end(); hite++)
		{
			H* h = *hite;
			if (!h) continue;
			// log("Processing hex (new edge/face creation) id: " + std::to_string(h->id())); // Already Commented out
			//find the faces without mark edges (sheet edges)
			std::vector<F*> fs; // Pair of parallel faces in the hex not touching the sheet edges
			F* found_hf = nullptr;
			for (int hfIndex = 0; hfIndex < h->fs.size(); hfIndex++)
			{
				F* hf = mesh->idFaces(h->fs[hfIndex]);
				if (!hf) continue; // Should not happen if collected properly
				bool without_mark_edge = true;
				for (int feIndex = 0; feIndex < hf->es.size(); feIndex++)
				{
					E* fe = mesh->idEdges(hf->es[feIndex]);
					if (fe && fe->sheetId()) // Check if edge belongs to the input sheet
					{
						without_mark_edge = false;
						break;
					}
				}
				if (without_mark_edge)
				{
					// Check if this face or its parallel face has already been processed via sheetId()
					if (hf->sheetId() == 0) { // Use sheetId() on faces to mark processed pairs
						F* parallel_f = mesh->f_parallel_f_in_hex(h, hf);
						if (parallel_f && parallel_f->sheetId() == 0) {
							hf->sheetId() = 1; // Mark this pair as processed
							parallel_f->sheetId() = 1;
							fs = { hf, parallel_f };
							found_hf = hf; // Keep track of which face we used as base
							break; // Found a valid pair
						}
					}
				}
			}
	
			// 如果没有找到符合条件的面，继续下一个六面体
			if (fs.empty() || !found_hf) {
				// log("  No suitable faces found for hex " + std::to_string(h->id())); // Commented out
				continue;
			}
	
			/*create the new edgs*/
			F* f = fs[0];
			std::vector<V*> FNewVs;
			for (int fvIndex = 0; fvIndex < f->vs.size(); fvIndex++)
			{
				V* fv = mesh->idVertices(f->vs[fvIndex]);
				V* nfv = mesh->idVertices(f->vs[(fvIndex + 1) % 4]);

				E* fe = mesh->VerticesEdge(fv, nfv);
				F* nf = mesh->flip_f(h, f, fe);
				E* nfe1 = mesh->flip_e(nf, fe, fv);
				V* newv1 = nfe1->newv();
				FNewVs.push_back(newv1);
				E* nfe2 = mesh->flip_e(nf, fe, nfv);
				V* newv2 = nfe2->newv();

				std::vector<E*> newe_parallel_e;
				newe_parallel_e.push_back(fe);
				E* fe_parallel_e = mesh->e_parallel_e_in_face(nf, fe);
				newe_parallel_e.push_back(fe_parallel_e);
				E* newe = mesh->VerticesEdge(newv1, newv2);
				if (newe == NULL)//construct new edge
				{
					std::vector<V*> ENewVs = { newv1, newv2 };
					newe = construct_new_edge(ENewVs);
					newes.push_back(newe);
					edge_sharp_policy(newe_parallel_e[0], newe_parallel_e[1], newe);
					edge_total_angle_policy_collapse(newe_parallel_e[0], newe_parallel_e[1], newe);
				}
				else
				{
					continue;
				}
				//change relation
				std::vector<E*> oldEs;
				oldEs.push_back(fe);
				nfv = mesh->flip_v(nfe2, nfv);
				nfe2 = mesh->flip_e(nf, nfe2, nfv);
				oldEs.push_back(nfe2);

				for (int oldeIndex = 0; oldeIndex < oldEs.size(); oldeIndex++)
				{
					E* oldE = oldEs[oldeIndex];
					//change face relation
					for (int oefIndex = 0; oefIndex < oldE->neighbor_fs.size(); oefIndex++)
					{
						F* oef = mesh->idFaces(oldE->neighbor_fs[oefIndex]);
						if (oef->is_delete()) continue;
						int efIndex = oef->edgeIndex(oldE->id());
						if (efIndex == -1)
						{
							continue;
						}
						oef->es[efIndex] = newe->id();
						newe->neighbor_fs.push_back(oef->id());
					}
					//change hex relation
					for (int oehIndex = 0; oehIndex < oldE->neighbor_hs.size(); oehIndex++)
					{
						H* oeh = mesh->idHexs(oldE->neighbor_hs[oehIndex]);
						if (oeh->is_delete()) continue;
						int ehIndex = oeh->edgeIndex(oldE->id());
						if (ehIndex == -1)
						{
							continue;
						}
						oeh->es[ehIndex] = newe->id();
						newe->neighbor_hs.push_back(oeh->id());
					}
				}
			}

			/*create new faces*/
			F* newf = constrcut_new_face(FNewVs);
			newfs.push_back(newf);
			for (int fIndex = 0; fIndex < fs.size(); fIndex++)
			{
				F* oldf = fs[fIndex];
				for (int fhIndex = 0; fhIndex < oldf->neighbor_hs.size(); fhIndex++)
				{
					H* fh = mesh->idHexs(oldf->neighbor_hs[fhIndex]);
					if (fh->is_delete()) continue;
					int fhfIndex = fh->faceIndex(oldf->id());
					fh->fs[fhfIndex] = newf->id();
					newf->neighbor_hs.push_back(fh->id());
				}
			}
		}
	
	
		// 在删除元素之前先完成所有需要访问映射表的操作 (e.g., attribute marking if uncommented)
		// log("Starting marking element attributes"); // Example if needed
		// try {
		// 	// log("New faces count: " + std::to_string(newfs.size()));
		// 	// mark_elements_attribute_collapse(newfs); // Call if needed
		// 	// log("Element attribute marking completed");
		// } catch (...) {
		//      log("Exception during attribute marking.");
		// }
	
	
		// --- Optimized Deletion ---
		log("Starting deletion operations");
	
		// 1. Erase from maps using the collected pointers/IDs
		int deleted_hexes_map = 0;
		for (H* h : delete_hs) {
			if (h && mesh->m_map_hexs.erase(h->id())) {
				deleted_hexes_map++;
			}
		}
		log("Removed " + std::to_string(deleted_hexes_map) + " hex entries from map");
	
		int deleted_faces_map = 0;
		for (F* f : delete_fs) {
			if (f && mesh->m_map_faces.erase(f->id())) {
				deleted_faces_map++;
			}
		}
		 log("Removed " + std::to_string(deleted_faces_map) + " face entries from map");
	
		int deleted_edges_map = 0;
		for (E* e : delete_es) {
			if (e && mesh->m_map_edges.erase(e->id())) {
				deleted_edges_map++;
			}
		}
		log("Removed " + std::to_string(deleted_edges_map) + " edge entries from map");
	
		int deleted_vertices_map = 0;
		for (V* v : delete_vs) {
			if (v && mesh->m_map_vertices.erase(v->id())) {
				deleted_vertices_map++;
			}
		}
		log("Removed " + std::to_string(deleted_vertices_map) + " vertex entries from map");
	
	
		// 2. Remove from lists using the pointer sets and delete objects
		int deleted_hexes_list = 0;
		mesh->hs.remove_if([&delete_hs, &deleted_hexes_list](H* h) {
			if (h == nullptr) return true; // Remove null pointers if any
			if (delete_hs.count(h)) {
				delete h; // Delete the object
				deleted_hexes_list++;
				return true; // Remove pointer from list
			}
			return false;
		});
		log("Deleted " + std::to_string(deleted_hexes_list) + " hexahedra objects and list entries");
	
		int deleted_faces_list = 0;
		mesh->fs.remove_if([&delete_fs, &deleted_faces_list](F* f) {
			if (f == nullptr) return true;
			if (delete_fs.count(f)) {
				delete f;
				deleted_faces_list++;
				return true;
			}
			return false;
		});
		log("Deleted " + std::to_string(deleted_faces_list) + " face objects and list entries");
	
		int deleted_edges_list = 0;
		mesh->es.remove_if([&delete_es, &deleted_edges_list](E* e) {
			if (e == nullptr) return true;
			if (delete_es.count(e)) {
				delete e;
				deleted_edges_list++;
				return true;
			}
			return false;
		});
		 log("Deleted " + std::to_string(deleted_edges_list) + " edge objects and list entries");
	
		int deleted_vertices_list = 0;
		mesh->vs.remove_if([&delete_vs, &deleted_vertices_list](V* v) {
			if (v == nullptr) return true;
			if (delete_vs.count(v)) {
				delete v;
				deleted_vertices_list++;
				return true;
			}
			return false;
		});
		log("Deleted " + std::to_string(deleted_vertices_list) + " vertex objects and list entries");
	
		log("Deletion operations completed");
		// 检查mesh 的 vertexes, edges, faces, hexs 的数量是否正确
		log("Mesh vertex count: " + std::to_string(mesh->vs.size()));
		log("Mesh edge count: " + std::to_string(mesh->es.size()));
		log("Mesh face count: " + std::to_string(mesh->fs.size()));
		log("Mesh hexahedron count: " + std::to_string(mesh->hs.size()));
		log("Sheet collapse execution completed");
	}

	template<typename M>
	double sheet_operation<M>::vector_angle(CPoint a, CPoint b)
	{
		a = a / a.norm();
		b = b / b.norm();
		double temp = a * b;
		temp = temp > 1 ? 1 : temp;
		temp = temp < -1 ? -1 : temp;
		double angle = acos(temp);
		angle = angle / MPI * 180;
		return angle;
	}

	template<typename M>
	void sheet_operation<M>::edge_angle(E* e)
	{
		// 检查输入的边指针是否有效
		if (!e) {
			// log("edge_angle: Error - Input edge pointer is null."); // 注释掉
			return;
		}
		// log("edge_angle: Calculating for Edge ID " + std::to_string(e->id()) + ", Boundary: " + (e->boundary() ? "Yes" : "No")); // 注释掉
	
		// 重置角度计算
		e->total_angle() = 0.0;
		e->ave_angle() = 0.0;
		// log("  edge_angle: Angles reset for Edge " + std::to_string(e->id())); // 注释掉
	
		// 非边界边处理
		if (!e->boundary())
		{
			e->total_angle() = 360.0;
			// 根据实际邻接六面体数量计算平均角，避免除零
			size_t neighbor_count = e->neighbor_hs.size();
			e->ave_angle() = (neighbor_count > 0) ? (360.0 / neighbor_count) : 90.0; // 默认为90
			// log("  edge_angle: Internal Edge " + std::to_string(e->id()) + ". Set total_angle=360.0, ave_angle=" + std::to_string(e->ave_angle())); // 注释掉
			return;
		}
	
		// 边界边计算
		double total_angle_sum = 0.0;
		int valid_angles_count = 0;
		// log("  edge_angle: Boundary Edge " + std::to_string(e->id()) + ". Processing " + std::to_string(e->neighbor_hs.size()) + " neighboring hexes."); // 注释掉
	
		// 遍历相邻六面体
		for (M::EHIterator ehite(mesh, e); !ehite.end(); ehite++)
		{
			H* eh = *ehite;
			if (!eh) {
				// log("    edge_angle: Warning - Could not find Hex for Edge " + std::to_string(e->id())); // 保留潜在的警告日志？或者也注释掉
				continue;
			}
			// log("    edge_angle: Processing Hex ID " + std::to_string(eh->id())); // 注释掉
	
			try {
				// 获取相邻面
				std::vector<F*> ehfs = mesh->e_adj_f_in_hex(eh, e);
				if (ehfs.size() != 2) { // 内部边可能有多个，但边界边在此六面体内应该只有2个相邻面
					// log("    edge_angle: Warning - Expected 2 adjacent faces in Hex " + std::to_string(eh->id()) + " for Edge " + std::to_string(e->id()) + ", found " + std::to_string(ehfs.size()) + ". Skipping hex."); // 注释掉
					continue;
				}
	
				F* f1 = ehfs[0];
				F* f2 = ehfs[1];
				if (!f1 || !f2) {
					 // log("    edge_angle: Warning - Null face pointer found in Hex " + std::to_string(eh->id()) + " for Edge " + std::to_string(e->id()) + ". Skipping hex."); // 注释掉
					 continue;
				}
				// log("    edge_angle: Adjacent faces: F" + std::to_string(f1->id()) + " and F" + std::to_string(f2->id())); // 注释掉
	
				CPoint n1_raw = f1->normal();
				CPoint n2_raw = f2->normal();
				double n1_norm = n1_raw.norm();
				double n2_norm = n2_raw.norm();
				// log("      edge_angle: Raw normals - n1_norm=" + std::to_string(n1_norm) + ", n2_norm=" + std::to_string(n2_norm)); // 注释掉
	
				// 检查法向量有效性
				if (n1_norm < 1e-10 || n2_norm < 1e-10) {
					// log("      edge_angle: Warning - Invalid normal detected (norm < 1e-10). Skipping angle calculation for this hex."); // 注释掉
					continue;
				}
	
				CPoint n1 = n1_raw;
				CPoint n2 = n2_raw;
	
				// 确定法线方向（指向六面体外部）
				bool n1_flipped = false;
				bool n2_flipped = false;
				if (f1->neighbor_hs.empty()){ // 增加空检查
					 // log("      edge_angle: Warning - Face " + std::to_string(f1->id()) + " has no hex neighbors."); // 注释掉
					 continue;
				} else if (f1->neighbor_hs[0] != eh->id()) { // 假设第一个邻居是当前六面体
					n1 = -n1;
					n1_flipped = true;
				}
				if (f2->neighbor_hs.empty()){ // 增加空检查
					// log("      edge_angle: Warning - Face " + std::to_string(f2->id()) + " has no hex neighbors."); // 注释掉
					 continue;
				} else if (f2->neighbor_hs[0] != eh->id()) {
					n2 = -n2;
					n2_flipped = true;
				}
				// log("      edge_angle: Orientation check - n1_flipped=" + std::string(n1_flipped ? "Yes" : "No") + ", n2_flipped=" + std::string(n2_flipped ? "Yes" : "No")); // 注释掉
	
	
				double angle_between_normals = vector_angle(n1, n2); // vector_angle 内部会处理归一化
				// log("      edge_angle: Angle between oriented normals = " + std::to_string(angle_between_normals) + " degrees."); // 注释掉
	
				// 计算内部角度
				double internal_angle = 180.0 - angle_between_normals;
	
				// 角度范围约束和有效性检查
				internal_angle = std::max(0.0, std::min(180.0, internal_angle)); // 边界边的角度应在0-180
				// log("      edge_angle: Calculated internal angle = " + std::to_string(internal_angle) + " degrees."); // 注释掉
	
				if (!std::isfinite(internal_angle)) {
					// log("      edge_angle: Warning - Calculated internal angle is not finite. Skipping this hex."); // 注释掉
					continue;
				}
	
				// 存储六面体内的边角（如果需要的话）
				// if (eh->edgeIndex(e->id()) != -1) {
				//      eh->edge_angle(e->id()) = internal_angle;
				// } else {
				//      log("      edge_angle: Warning - Edge ID " + std::to_string(e->id()) + " not found in Hex ID " + std::to_string(eh->id()) + " edge list."); // 注释掉
				// }
	
				total_angle_sum += internal_angle;
				valid_angles_count++;
				// log("      edge_angle: Added angle to sum. Current sum = " + std::to_string(total_angle_sum) + ", valid hex count = " + std::to_string(valid_angles_count)); // 注释掉
			}
			catch (const std::exception& ex) {
				// log("Exception in edge angle calculation for hex " + std::to_string(eh->id()) + ": " + std::string(ex.what())); // 保留或注释掉异常日志
			}
			catch (...) {
				// log("Unknown exception in edge angle calculation for hex " + std::to_string(eh->id())); // 保留或注释掉异常日志
			}
		} // 结束遍历相邻六面体
	
		// 处理没有有效角度计算的情况
		if (valid_angles_count == 0) {
			 // log("  edge_angle: Warning - No valid hex neighbors contributed angles for Boundary Edge " + std::to_string(e->id()) + ". Setting default angles."); // 注释掉
			 e->total_angle() = 180.0; // 边界边默认180
			 e->ave_angle() = 90.0; // 默认平均90
			 return;
		}
	
		// 最终角度值验证
		if (!std::isfinite(total_angle_sum) || total_angle_sum < 0) {
			 // log("  edge_angle: Warning - Final total_angle_sum (" + std::to_string(total_angle_sum) + ") is invalid for Edge " + std::to_string(e->id()) + ". Resetting to 180.0."); // 注释掉
			 total_angle_sum = 180.0; // 边界边默认180
		}
	
		 e->total_angle() = total_angle_sum;
		 e->ave_angle() = total_angle_sum / valid_angles_count; // 平均值
	
		// log("  edge_angle: Final calculated angles for Edge " + std::to_string(e->id()) + ": total_angle=" + std::to_string(e->total_angle()) + ", ave_angle=" + std::to_string(e->ave_angle())); // 注释掉
	
	}

	template<typename M>
	void sheet_operation<M>::edge_ideal_degree(E* e)
	{
		// 检查输入的边指针是否有效
		if (!e) {
			// log("edge_ideal_degree: Error - Input edge pointer is null."); // 注释掉
			return;
		}
	
		try {
			// log("Computing ideal degree for edge ID: " + std::to_string(e->id()) + " (Unconditional calculation)"); // 注释掉
	
			if (e->boundary())
			{
				if (e->sharp() > 0) // 检查 sharp() 返回值是否大于 0
				{
					// 检查总角度是否合理
					// total_angle() 在边界锐利边上可能不是 0-360，取决于具体定义
					// 通常锐利边界边的 ideal_degree 由几何决定，比如 90度角对应 degree 1, 270度角对应 degree 3
					// 如果 total_angle 代表的是材料内部的角度（通常<180）
					if (!std::isfinite(e->total_angle()) || e->total_angle() < 0 || e->total_angle() > 360) { // 放宽检查范围，但仍需有限
						// log("Warning: Edge " + std::to_string(e->id()) + " has abnormal total angle value: " + std::to_string(e->total_angle()) + " for ideal degree calculation"); // 注释掉
						// 使用默认边界理想度
						e->ideal_degree() = 2;
					}
					else {
						// 基于角度计算理想度
						int ideal_degree = round(e->total_angle() / 90.0);
						// 确保理想度至少为 1 （对于边界边）
						ideal_degree = (ideal_degree <= 0) ? 1 : ideal_degree;
						// 边界锐利边的理想度上限通常是 3 (对应 > 270 度的内角)
						ideal_degree = (ideal_degree > 3) ? 3 : ideal_degree;
						e->ideal_degree() = ideal_degree;
						// log("Sharp boundary edge " + std::to_string(e->id()) + " ideal degree set to: " + std::to_string(ideal_degree) + " (total angle: " + std::to_string(e->total_angle()) + ")"); // 注释掉
	
						// 更新 feature_ideal_degree 映射表
						if (mesh && mesh->feature_ideal_degree.find(e->sharp()) == mesh->feature_ideal_degree.end())
						{
							// std::vector<E*> sharps; // 这两行似乎不需要
							// sharps.push_back(e);
							mesh->feature_ideal_degree.insert(std::pair<int, int>(e->sharp(), ideal_degree));
							// log("Added feature ideal degree for sharp ID " + std::to_string(e->sharp()) + ": " + std::to_string(ideal_degree)); // 注释掉
						}
					}
				}
				else // 非锐利边界边
				{
					e->ideal_degree() = 2; // 理想度为 2
					// log("Regular boundary edge " + std::to_string(e->id()) + " ideal degree set to: 2"); // 注释掉
				}
			}
			else // 内部边
			{
				e->ideal_degree() = 4; // 内部边理想度为 4
				// log("Non-boundary edge " + std::to_string(e->id()) + " ideal degree set to: 4"); // 注释掉
			}
	
			// log("Edge " + std::to_string(e->id()) + " final ideal degree: " + std::to_string(e->ideal_degree())); // 注释掉
		}
		catch (const std::exception& ex) {
			// log("Exception in edge_ideal_degree for edge " + std::to_string(e->id()) + ": " + std::string(ex.what())); // 保留或注释掉异常日志
			// 设置安全的默认值
			if(e) e->ideal_degree() = e->boundary() ? 2 : 4;
		}
		catch (...) {
			// log("Unknown exception in edge_ideal_degree for edge " + std::to_string(e->id())); // 保留或注释掉异常日志
			// 设置安全的默认值
			 if(e) e->ideal_degree() = e->boundary() ? 2 : 4;
		}
	}

	template<typename M>
	void sheet_operation<M>::compute_edge_energy()
	{
		// log("Starting edge energy computation"); // 注释掉
		// 确保 mesh 指针有效
		if(!mesh) {
			 // 如果需要，可以在此处添加错误处理
			 return;
		}
		mesh->computeNormal(); // 计算法线是后续计算的基础
	
		// int totalEdges = 0; // 移除调试用计数器
		// int boundaryEdges = 0;
		// int sharpEdges = 0;
		// int singularEdges = 0;
	
		// log("Computing normals and angles for each edge"); // 注释掉
		for (M::MEIterator eite(mesh); !eite.end(); eite++)
		{
			E* e = *eite;
			if (!e) continue; // 增加对空指针边的检查
			// totalEdges++;
	
			try {
				// // 记录进度日志 (已注释)
				// if (totalEdges % 1000 == 0) {
				//     log("Processing edge " + ...);
				// }
	
				// 计算角度 (假设 edge_angle 是可访问的成员函数或辅助函数)
				this->edge_angle(e);
	
				// 计算理想度数
				this->edge_ideal_degree(e);
	
				// 计算简化能量
				// 确保 neighbor_hs.size() 正确反映了边的度数
				e->sim_energy() = (int)e->ideal_degree() - (int)e->neighbor_hs.size();
	
				// // 统计特殊边的数量 (已注释)
				// if (e->boundary()) boundaryEdges++;
				// if (e->sharp() > 0) sharpEdges++; // 检查 sharp > 0
				// if (e->singularity()) singularEdges++;
	
				// // 记录特殊边的详细信息 (已注释)
				// if (e->boundary() || e->sharp() > 0 || e->singularity()) { ... }
	
			}
			catch (const std::exception& ex) {
				// 可以保留最小的异常日志记录
				// log("Exception processing edge ID " + std::to_string(e->id()) + ": " + std::string(ex.what())); // 注释掉
				// 考虑添加错误处理，例如设置默认能量
				if(e) e->sim_energy() = 0.0; // 示例：设置默认能量
			}
			catch (...) {
				// 可以保留最小的异常日志记录
				// log("Unknown exception processing edge ID " + std::to_string(e->id())); // 注释掉
				if(e) e->sim_energy() = 0.0; // 示例：设置默认能量
			}
		}
	
		// // 最终的统计日志 (已注释)
		// log("Edge energy computation completed");
		// log("Total edges: " + ...);
	}

	template<typename M>
	std::vector<std::vector<typename sheet_operation<M>::E*>> sheet_operation<M>::get_sheet_parallel_edges(std::vector<E*> sheet)
	{
		// log("Starting get_sheet_parallel_edges, sheet size: " + std::to_string(sheet.size())); // 注释掉
	
		std::vector<H*> hs; // 用于后续重置标记
		std::vector<std::vector<E*>> parallel_es;
	
		// 1. 重置与 sheet 相关的六面体和边的标记 (mark 用于防止重复处理)
		// log("Initializing hexahedra marks..."); // 注释掉
		std::set<int> processed_hex_ids; // 跟踪已处理的六面体ID
		for (int eIndex = 0; eIndex < sheet.size(); eIndex++) {
			E* e = sheet[eIndex];
			if (!e) continue; // 跳过空指针边
	
			for (int ehIndex = 0; ehIndex < e->neighbor_hs.size(); ehIndex++) {
				int hex_id = e->neighbor_hs[ehIndex];
				// 如果这个六面体ID已经被处理过，则跳过
				if (processed_hex_ids.count(hex_id)) {
					continue;
				}
	
				H* eh = mesh->idHexs(hex_id);
				if (eh) {
					eh->mark() = false; // 重置六面体标记
					processed_hex_ids.insert(hex_id); // 记录已处理
					hs.push_back(eh); // 添加到列表以便后续重置边的标记
	
					// 同时重置该六面体内所有边的标记
					for(int he_id : eh->es) {
						E* he = mesh->idEdges(he_id);
						if(he) he->mark() = false;
					}
				}
			}
		}
		// log("Finished initializing marks for " + std::to_string(processed_hex_ids.size()) + " hexahedra"); // 注释掉
	
		// 2. 查找平行边对
		// int edgesProcessed = 0; // 移除计数器
		// int hexsFound = 0;
		// int parallelPairsFound = 0;
	
		for (int eIndex = 0; eIndex < sheet.size(); eIndex++)
		{
			try {
				E* e = sheet[eIndex];
				// 跳过无效边
				if (!e || e->vs.size() != 2) {
					// log("  Warning: Null or invalid edge at index " + std::to_string(eIndex)); // 注释掉
					continue;
				}
				// edgesProcessed++;
	
				// 遍历邻接六面体
				for (int ehIndex = 0; ehIndex < e->neighbor_hs.size(); ehIndex++)
				{
					try {
						H* eh = mesh->idHexs(e->neighbor_hs[ehIndex]);
						// 跳过无效六面体或已处理过的六面体
						if (!eh || eh->mark()) { // 使用 mark() 检查是否已处理
							// log("    Skipping null or already processed hexahedron ID: " + std::to_string(e->neighbor_hs[ehIndex])); // 注释掉
							continue;
						}
	
						eh->mark() = true; // 标记当前六面体为已处理
						// hs.push_back(eh); // hs 列表在初始重置时已收集
						// hexsFound++;
	
						// 获取六面体内与边 e 邻接的面
						std::vector<F*> nfs = mesh->e_adj_f_in_hex(eh, e);
						if (nfs.size() < 1 || !nfs[0]) { // 需要至少一个有效面
							// log("    Warning: No valid adjacent faces found for edge " + std::to_string(e->id()) + " in hex " + std::to_string(eh->id())); // 注释掉
							continue;
						}
						F* nf = nfs[0]; // 使用第一个邻接面开始
	
						// 获取边的顶点
						V* ev1 = mesh->idVertices(e->vs[0]);
						V* ev2 = mesh->idVertices(e->vs[1]);
						if (!ev1 || !ev2) {
							// log("    Warning: Unable to get edge vertices for edge " + std::to_string(e->id())); // 注释掉
							continue;
						}
	
						// 获取起始的翻转边
						E* ne1 = mesh->flip_e(nf, e, ev1);
						E* ne2 = mesh->flip_e(nf, e, ev2);
						if (!ne1 || !ne2) {
							// log("    Warning: Failed to get initial flip edges for edge " + std::to_string(e->id())); // 注释掉
							continue;
						}
	
						// 获取起始的翻转面
						F* nf1 = mesh->flip_f(eh, nf, ne1);
						F* nf2 = mesh->flip_f(eh, nf, ne2);
						if (!nf1 || !nf2) {
							// log("    Warning: Failed to get initial flip faces for edges " + std::to_string(ne1->id()) + "/" + std::to_string(ne2->id())); // 注释掉
							continue;
						}
	
						// 循环翻转4次查找平行边对
						for (int peIndex = 0; peIndex < 4; peIndex++)
						{
							try {
								// 执行翻转操作，并检查每一步的返回指针
								ne1 = mesh->flip_e(nf1, ne1, ev1);
								if (!ne1) break;
								ev1 = mesh->flip_v(ne1, ev1);
								if (!ev1) break;
	
								ne2 = mesh->flip_e(nf2, ne2, ev2);
								if (!ne2) break;
								ev2 = mesh->flip_v(ne2, ev2);
								if (!ev2) break;
	
								// 检查找到的平行边之一 (ne1) 是否已被标记（表示这个平行关系已找到）
								if (ne1->mark()) {
									continue; // 如果标记了，说明这对平行边已处理
								}
	
								// 标记这对平行边（ne1 和 ne2）为已找到
								ne1->mark() = true;
								ne2->mark() = true;
	
								std::vector<E*> pair_es = { ne1, ne2 };
								parallel_es.push_back(pair_es);
								// parallelPairsFound++;
							}
							catch (...) { // 捕获翻转过程中的任何异常
								// log("      Exception during parallel edges search iteration " + ...); // 注释掉
								break; // 出现异常则停止当前六面体的查找
							}
						} // 结束翻转循环
					}
					catch (...) { // 捕获处理单个六面体时的异常
						// log("    Exception processing adjacent hex at index " + ...); // 注释掉
					}
				} // 结束遍历邻接六面体
			}
			catch (...) { // 捕获处理单个 sheet 边时的异常
				// log("  Exception processing edge at index " + ...); // 注释掉
			}
		} // 结束遍历 sheet 边
	
		// log("Processing complete. Found " + std::to_string(parallel_es.size()) + " parallel edge pairs."); // 注释掉
	
		// 3. 重置所有处理过的六面体及其内部边的标记
		// log("Resetting marks on " + std::to_string(hs.size()) + " hexahedra"); // 注释掉
		for (H* h : hs)
		{
			if (!h) continue;
			h->mark() = false; // 重置六面体标记
			for (int he_id : h->es) {
				E* he = mesh->idEdges(he_id);
				if (he) {
					he->mark() = false; // 重置边标记
				}
			}
		}
	
		return parallel_es;
	}

	template<typename M>
	double sheet_operation<M>::predict_sheet_collapse_energy(std::vector<E*> sheet)
	{
		// log("======= Starting sheet collapse energy prediction ======="); // 注释掉
		// log("Input sheet size: " + std::to_string(sheet.size()) + " edges"); // 注释掉
	
		// 检查 mesh 指针和空 sheet
		if (!mesh || sheet.empty()) {
			// log("predict_sheet_collapse_energy: Error - Mesh pointer is null or input sheet is empty."); // 可选保留
			return -99999.0; // 使用常量或明确的错误码
		}
	
		double original_total_energy = 0;
		double predict_total_energy = 0;
		const double cannot_collapse = -99999.0; // 定义不可折叠的返回值
	
		// 1. 检查所有 sheet 边的顶点约束条件
		for (int eIndex = 0; eIndex < sheet.size(); eIndex++)
		{
			E* e = sheet[eIndex];
			// 检查边和顶点数量
			if (!e || e->vs.size() != 2) {
				// log("Warning: Skipping invalid edge in vertex constraint check."); // 可选保留
				continue; // 跳过无效边
			}
	
			V* ev1 = mesh->idVertices(e->vs[0]);
			V* ev2 = mesh->idVertices(e->vs[1]);
	
			// 检查顶点指针
			if (!ev1 || !ev2) {
				 // log("Warning: Edge " + std::to_string(e->id()) + " has null vertex pointer(s). Cannot check constraints."); // 可选保留
				 return cannot_collapse; // 如果顶点无效，无法安全地检查约束
			}
	
			// // 检查并重置异常的 feature_vertex 值 (可选，如果数据可能不稳定则保留)
			// if (ev1->feature_vertex() < -1000000 || ev1->feature_vertex() > 1000000) {
			//      ev1->feature_vertex() = 0;
			// }
			// if (ev2->feature_vertex() < -1000000 || ev2->feature_vertex() > 1000000) {
			//      ev2->feature_vertex() = 0;
			// }
	
			// log("Checking vertex constraints for edge " + std::to_string(e->id())); // 注释掉
			// ... (其他详细顶点属性日志注释掉) ...
	
			// 执行约束检查
			if (ev1->corner() && ev2->corner()) {
				// log("Both vertices are corners..."); // 注释掉
				return cannot_collapse;
			}
			if (ev1->feature_vertex() < 0 && ev2->feature_vertex() < 0) {
				// log("Both vertices have feature values < 0..."); // 注释掉
				return cannot_collapse;
			}
			// 检查特征边和角点的连接性 (需要 mesh->feature_edge_corner)
			if ( (ev1->feature_vertex() < 0 && ev2->feature_vertex() > 0) ||
				 (ev1->feature_vertex() > 0 && ev2->feature_vertex() < 0) )
			{
				int sharp_id = (ev1->feature_vertex() > 0) ? ev1->feature_vertex() : ev2->feature_vertex();
				int corner_id = (ev1->feature_vertex() < 0) ? ev1->feature_vertex() : ev2->feature_vertex();
				// log("  Vertex feature mismatch (one < 0, one > 0)..."); // 注释掉
	
				try {
					 // 确保 map 访问安全
					if (mesh->feature_edge_corner.count(sharp_id)) {
						const std::vector<int>& connect_corners_id = mesh->feature_edge_corner.at(sharp_id);
						// log("  Connected corners count: " + ...); // 注释掉
						if (std::find(connect_corners_id.begin(), connect_corners_id.end(), corner_id) == connect_corners_id.end()) {
							// log("  Corner not connected to feature edge..."); // 注释掉
							return cannot_collapse;
						} // else { log("  Corner connected..."); } // 注释掉
					} else {
						 // log("  Sharp ID " + ... + " not found in map..."); // 可选保留警告
						 return cannot_collapse;
					}
				} catch (...) { // 捕获潜在异常
					 // log("  Exception while checking connection..."); // 注释掉
					 return cannot_collapse;
				}
			}
			// 检查两个特征点是否属于同一特征线
			if (ev1->feature_vertex() > 0 && ev2->feature_vertex() > 0)
			{
				// log("  Both vertices have feature values > 0"); // 注释掉
				if (ev1->feature_vertex() != ev2->feature_vertex()) {
					// log("  Feature IDs don't match..."); // 注释掉
					return cannot_collapse;
				} // else { log("  Feature IDs match..."); } // 注释掉
			}
		} // 结束顶点约束检查循环
	
		// log("Vertex constraint checks passed, continuing"); // 注释掉
	
		// 2. 获取平行边并评估能量
		try {
			// log("Starting to get parallel edges in sheet"); // 注释掉
			std::vector<std::vector<E*>> parallel_es = get_sheet_parallel_edges(sheet); // 假设此函数健壮
			// log("Retrieved " + std::to_string(parallel_es.size()) + " pairs of parallel edges"); // 注释掉
	
			int temp_less_num = 0;
			int temp_equal_num = 0;
			int temp_large_num = 0;
			int temp_boundary_num = 0; // 边界折叠估计参数
	
			for (int peIndex = 0; peIndex < parallel_es.size(); peIndex++)
			{
				// 检查平行边对是否有效
				if(parallel_es[peIndex].size() != 2) continue;
				E* e1 = parallel_es[peIndex][0];
				E* e2 = parallel_es[peIndex][1];
				if (!e1 || !e2) continue; // 空指针检查
	
				// log("\nProcessing parallel edge pair #" + ...); // 注释掉
				// ... (详细边属性日志注释掉) ...
	
				// 检查是否两条边都是锐利边
				if (e1->sharp() > 0 && e2->sharp() > 0) {
					// log("  Both edges are sharp, cannot collapse..."); // 注释掉
					return cannot_collapse;
				}
	
				double predict_angle = 0;
				double predict_degree = 0;
	
				// 计算预测角度
				// log("  Calculating predicted angle"); // 注释掉
				if (e1->boundary() && e2->boundary()) {
					if (e1->sharp() > 0 && !(e2->sharp() > 0)) predict_angle = e1->total_angle();
					else if (!(e1->sharp() > 0) && e2->sharp() > 0) predict_angle = e2->total_angle();
					else predict_angle = (e1->total_angle() + e2->total_angle()) * 0.5;
				} else if (e1->boundary()) {
					predict_angle = e1->total_angle();
				} else if (e2->boundary()) {
					predict_angle = e2->total_angle();
				} else {
					predict_angle = (e1->total_angle() + e2->total_angle()) * 0.5;
				}
				// log("    Predicted angle: " + ...); // 注释掉
	
	
				// 计算预测度数
				// log("  Calculating predicted degree"); // 注释掉
				if (e1->boundary() && e2->boundary()) {
					predict_degree = e1->neighbor_hs.size() + e2->neighbor_hs.size() - 2;
				} else {
					predict_degree = e1->neighbor_hs.size() + e2->neighbor_hs.size() - 4;
				}
				// log("    Predicted degree: " + ...); // 注释掉
	
				// 检查度数约束
				// log("  Checking if degree is too small"); // 注释掉
				if (e1->boundary() || e2->boundary()) {
					if (predict_degree < 1) {
						// log("    Contains boundary edge and predicted degree < 1..."); // 注释掉
						return cannot_collapse;
					}
				} else {
					if (predict_degree < 2) {
						// log("    Predicted degree < 2..."); // 注释掉
						return cannot_collapse;
					}
				}
	
				// 计算能量变化
				// 基于预测角度推算预测的理想度数
				double ideal_degree_pred_angle_f = predict_angle / 90.0;
				int ideal_degree_pred = static_cast<int>(round(ideal_degree_pred_angle_f));
				if (!e1->boundary() && !e2->boundary()) { // 内部边理想度为4
					 ideal_degree_pred = 4;
				} else { // 边界边理想度至少为1，通常不超过3
					 ideal_degree_pred = (ideal_degree_pred <= 0) ? 1 : ideal_degree_pred;
					 ideal_degree_pred = (ideal_degree_pred > 3) ? 3 : ideal_degree_pred; // 限制边界理想度上限
				}
	
				int temp_predict_energy = abs(ideal_degree_pred - static_cast<int>(predict_degree)); // 使用推算的理想度
				// log("  Calculating predicted energy: abs(" + ...); // 注释掉
				predict_total_energy += temp_predict_energy;
	
				int temp_original_energy1 = abs(e1->sim_energy());
				int temp_original_energy2 = abs(e2->sim_energy());
				original_total_energy += temp_original_energy1;
				original_total_energy += temp_original_energy2;
	
				// log("  Cumulative original energy: " + ...); // 注释掉
				// log("  Edge1 energy: " + ...); // 注释掉
	
				// 估计边界能量变化
				int boundary_energy1 = 0;
				int boundary_energy2 = 0;
				if (e1->boundary()) {
					boundary_energy1 = temp_original_energy1 - temp_predict_energy;
					// log("  Edge1 boundary energy calculation: " + ...); // 注释掉
				}
				if (e2->boundary()) {
					boundary_energy2 = temp_original_energy2 - temp_predict_energy;
					// log("  Edge2 boundary energy calculation: " + ...); // 注释掉
				}
	
				// 更新边界折叠估计参数
				if ((e1->boundary() || e2->boundary()) && !(e1->boundary() && e2->boundary())) {
					// log("  One edge is boundary while the other is not..."); // 注释掉
					if (temp_boundary_num <= 0) { // 仅当还不是正数时更新
						if (boundary_energy1 > 0 || boundary_energy2 > 0) temp_boundary_num = 1;
						else if (boundary_energy1 == 0 && boundary_energy2 == 0) temp_boundary_num = 0;
						else temp_boundary_num = -1;
						// log("    Updated temp_boundary_num: " + ...); // 注释掉
					}
				}
	
				// 估计折叠能量比较结果
				// log("  Estimating collapse energy"); // 注释掉
				if (temp_original_energy1 > temp_predict_energy || temp_original_energy2 > temp_predict_energy) temp_large_num++;
				else if (temp_original_energy1 == temp_predict_energy || temp_original_energy2 == temp_predict_energy) temp_equal_num++;
				else if (temp_original_energy1 < temp_predict_energy && temp_original_energy2 < temp_predict_energy) temp_less_num++;
				// log("    large/equal/less counts: " + ...); // 注释掉
	
			} // 结束遍历平行边对
	
			// log("\nCollapse energy calculation summary:"); // 注释掉
			// ... (总结日志注释掉) ...
	
			// 最终决策
			if (temp_boundary_num < 0) {
				// log("Boundary energy check: temp_boundary_num < 0..."); // 注释掉
				return cannot_collapse;
			} else if (temp_large_num > 0) {
				// log("Energy comparison check: temp_large_num > 0..."); // 注释掉
				return original_total_energy - predict_total_energy;
			} else if (temp_less_num == 0) { // 包括 temp_equal_num > 0 的情况
				// log("Energy comparison check: temp_less_num = 0..."); // 注释掉
				return original_total_energy - predict_total_energy;
			} else { // temp_less_num > 0
				// log("Other cases, not suitable for collapse..."); // 注释掉
				return cannot_collapse;
			}
		}
		catch (const std::exception& e) {
			// log("Exception occurred during collapse energy prediction: " + ...); // 保留或注释掉
			return cannot_collapse;
		}
		catch (...) {
			// log("Unknown exception occurred during collapse energy prediction"); // 保留或注释掉
			return cannot_collapse;
		}
	}

	template <typename M>
	inline double sheet_operation<M>::get_sheet_on_boundary_ratio(std::vector<E *> sheet)
	{
		// 检查输入的 sheet 是否为空
		if (sheet.empty()) {
			// 如果 sheet 为空，边界比例为 0
			return 0.0;
		}
	
		int is_boundary_count = 0;
		int processed_count = 0; // 用于统计实际处理的有效边数
	
		// 遍历 sheet 中的所有边
		for (int i = 0; i < sheet.size(); i++)
		{
			E *e = sheet[i];
	
			// 跳过空指针，提高健壮性
			if (!e) {
				continue;
			}
			processed_count++; // 计数有效的边
	
			// 检查边的 boundary 标志
			if (e->boundary())
			{
				is_boundary_count++;
			}
		}
	
		// 计算比例，使用实际处理的边数作为分母
		if (processed_count == 0) {
			 // 如果没有处理任何有效的边（例如 sheet 只包含空指针）
			return 0.0;
		}
	
		// 确保进行浮点数除法
		double ratio = static_cast<double>(is_boundary_count) / static_cast<double>(processed_count);
		return ratio;
	}

	template <typename M>
    inline double sheet_operation<M>::get_sheet_on_feature_ratio(std::vector<E*> sheet)
    {
        if (sheet.empty()) {
            log("get_sheet_on_feature_ratio: Warning - Input sheet is empty.");
            return 0.0;
        }

        // Assuming the first edge's sheet ID represents the whole sheet for logging
        int representative_sheet_id = sheet[0] ? sheet[0]->sheet() : -1;
        log("get_sheet_on_feature_ratio: Calculating for Sheet ID (representative) " + std::to_string(representative_sheet_id) + " with " + std::to_string(sheet.size()) + " edges.");

        int is_feature_count = 0;
        int processed_count = 0;
        for (int i = 0; i < sheet.size(); i++)
        {
            E* e = sheet[i];
            if (!e) {
                log("  get_sheet_on_feature_ratio: Skipping null edge pointer at index " + std::to_string(i));
                continue;
            }
            processed_count++;
            int sharp_value = e->sharp(); // Read the value once

            // Log every edge's sharp value for debugging
            log("  get_sheet_on_feature_ratio: Checking Edge ID " + std::to_string(e->id()) + ", sharp() = " + std::to_string(sharp_value));

            // The check should likely be for non-zero sharp value
            if (sharp_value > 0) // Check if sharp value indicates it's a feature
            {
                is_feature_count++;
                log("    -> Edge ID " + std::to_string(e->id()) + " counted as feature.");
            }
        }

        if (processed_count == 0) {
             log("get_sheet_on_feature_ratio: Warning - No valid edges processed in the sheet.");
             return 0.0;
        }

        // Ensure floating point division
        double ratio = static_cast<double>(is_feature_count) / static_cast<double>(processed_count);
        log("get_sheet_on_feature_ratio: Calculation complete for Sheet ID " + std::to_string(representative_sheet_id) + ". Feature count = " + std::to_string(is_feature_count) + ", Processed edges = " + std::to_string(processed_count) + ", Ratio = " + std::to_string(ratio));

        return ratio;
    	}

		template <typename M>
		bool sheet_operation<M>::sheet_endpoints_are_corners_or_boundary(std::vector<E*> sheet)
		{
			// 检查 mesh 指针和空 sheet
			if (!mesh || sheet.empty()) {
				// log("Warning: sheet_endpoints_are_corners_or_boundary called with empty sheet or null mesh.");
				return false;
			}
		
			std::map<int, int> vertex_counts;
			std::set<int> sheet_vertex_ids;
		
			// 1. 统计顶点出现次数
			for (E* e : sheet) {
				// 跳过无效边
				if (!e || e->vs.size() != 2) continue;
				for (int vid : e->vs) {
					// 可以在这里检查 vid 的有效性，但 set 会自动处理重复
					sheet_vertex_ids.insert(vid);
					vertex_counts[vid]++;
				}
			}
		
			// 2. 查找端点顶点 (计数为1的顶点)
			std::vector<V*> endpoint_vertices;
			for (int vid : sheet_vertex_ids) {
				// 使用 find 避免自动创建条目
				auto it = vertex_counts.find(vid);
				if (it != vertex_counts.end() && it->second == 1) {
					V* v = mesh->idVertices(vid);
					// 检查 idVertices 是否返回有效指针
					if (v) {
						endpoint_vertices.push_back(v);
					} else {
						// log("Warning: Could not find vertex with ID " + std::to_string(vid) + " while finding sheet endpoints.");
					}
				}
			}
		
			// 3. 判断端点情况
			//    通常非闭合 sheet 有 2 个端点
			//    闭合 sheet 或其他复杂情况 (如多段、断裂) 端点数不为 2
			if (endpoint_vertices.size() != 2) {
				 // 可以选择性地检查是否为闭环 (所有顶点计数都为2)
				// bool is_loop = true;
				// if (!sheet_vertex_ids.empty()) { // 只有在有顶点时才检查
				//     for(int vid : sheet_vertex_ids) {
				//         auto it = vertex_counts.find(vid);
				//         if(it == vertex_counts.end() || it->second != 2) { // 如果顶点不存在或计数不为2
				//             is_loop = false;
				//             break;
				//         }
				//     }
				// } else {
				//     is_loop = false; // 空sheet不是闭环
				// }
				// if(is_loop) {
				//      log("Sheet appears to be a loop. Returning false for endpoint check.");
				// } else {
				//     log("Warning: Found " + std::to_string(endpoint_vertices.size()) + " endpoints (expected 2). Returning false.");
				// }
				return false; // 对于非标准端点数的情况，保守返回 false
			}
		
			// 4. 检查两个端点的属性
			bool endpoint1_special = false;
			bool endpoint2_special = false;
		
			// 检查 endpoint_vertices[0] 是否有效，然后检查其属性
			V* ep0 = endpoint_vertices[0];
			if (ep0) {
				endpoint1_special = ep0->boundary() || ep0->corner() || (ep0->feature_vertex() != 0);
				// log("Endpoint 1 (ID: " + std::to_string(ep0->id()) + ") Special: " + (endpoint1_special ? "Y":"N"));
			}
		
			// 检查 endpoint_vertices[1] 是否有效，然后检查其属性
			V* ep1 = endpoint_vertices[1];
			if (ep1) {
				 endpoint2_special = ep1->boundary() || ep1->corner() || (ep1->feature_vertex() != 0);
				 // log("Endpoint 2 (ID: " + std::to_string(ep1->id()) + ") Special: " + (endpoint2_special ? "Y":"N"));
			}
		
			// 如果两个端点指针都有效，则返回它们是否至少有一个是 special
			// 如果有任何一个端点指针无效，则之前的检查已经返回false或在这个逻辑中保持false
			return endpoint1_special || endpoint2_special;
		}

	// --- IMPLEMENTATION for get_sheet_adjacent_feature_edges_count ---
	template <typename M>
	int sheet_operation<M>::get_sheet_adjacent_feature_edges_count(std::vector<E*> sheet)
        {
             if (sheet.empty()) {
                 log("Warning: get_sheet_adjacent_feature_edges_count called with empty sheet.");
                 return 0;
             }

             std::set<int> sheet_edge_ids;
             std::set<int> sheet_vertex_ids;
             std::set<int> counted_feature_edge_ids; // Track adjacent feature edges already counted
             int total_adjacent_feature_count = 0;

             // Populate sets of edge and vertex IDs belonging to the sheet
             for (E* e : sheet) {
                 if (!e || e->vs.size() != 2) continue;
                 sheet_edge_ids.insert(e->id());
                 sheet_vertex_ids.insert(e->vs[0]);
                 sheet_vertex_ids.insert(e->vs[1]);
             }

             // Iterate through vertices belonging to the sheet
             for (int vid : sheet_vertex_ids) {
                 V* v = mesh->idVertices(vid);
                 if (!v) continue; // Skip if vertex doesn't exist

                 // Iterate through edges neighboring this vertex
                 for (int neighbor_eid : v->neighbor_es) {
                     // Skip if the neighboring edge is part of the sheet itself
                     if (sheet_edge_ids.count(neighbor_eid)) {
                         continue;
                     }

                     // Skip if this adjacent feature edge has already been counted
                     if (counted_feature_edge_ids.count(neighbor_eid)) {
                         continue;
                     }

                     E* neighbor_e = mesh->idEdges(neighbor_eid);
                     if (!neighbor_e) continue; // Skip if edge doesn't exist

                     // Check if the neighboring edge is sharp/feature
                     if (neighbor_e->sharp() > 0) { // Assuming sharp > 0 means it's a feature edge
                         total_adjacent_feature_count++;
                         counted_feature_edge_ids.insert(neighbor_eid); // Mark as counted
                         log("Adjacent feature edge found: ID " + std::to_string(neighbor_e->id()) + " connected to sheet vertex ID " + std::to_string(vid));
                     }
                 }
             }
             log("Total adjacent feature edges count for sheet: " + std::to_string(total_adjacent_feature_count));
             return total_adjacent_feature_count;
        }
    template <typename M>
    inline double sheet_operation<M>::get_sheet_curvature_metric(std::vector<E *> sheet)
    {
        if (sheet.empty()) {
            return 0.0;
        }

        try {
            log("计算Sheet曲率度量(Curvature Metric)...");

            // 方法1：计算端点之间的欧氏距离与路径几何长度的比值
            // 获取sheet两端顶点
            E* first_edge = sheet.front();
            E* last_edge = sheet.back();

            // 检查边是否有效
            if (!first_edge || !last_edge || first_edge->vs.size() < 2 || last_edge->vs.size() < 2) {
                log("警告: Sheet端点边无效");
                return 0.0;
            }

            // 获取sheet两端的顶点
            V* start_v = mesh->idVertices(first_edge->vs[0]);
            V* end_v = mesh->idVertices(last_edge->vs[1]);

            // 确保能找到两端顶点
            if (!start_v || !end_v) {
                log("警告: 无法获取sheet端点顶点");
                return 0.0;
            }

            // 计算欧氏距离
            CPoint start_pos = start_v->position();
            CPoint end_pos = end_v->position();
            double euclidean_distance = (start_pos - end_pos).norm();

            // 计算sheet的总几何长度
            double total_length = 0.0;
            for (auto e : sheet) {
                if (!e || e->vs.size() < 2) continue;
                
                V* v1 = mesh->idVertices(e->vs[0]);
                V* v2 = mesh->idVertices(e->vs[1]);
                
                if (!v1 || !v2) continue;
                
                CPoint p1 = v1->position();
                CPoint p2 = v2->position();
                total_length += (p1 - p2).norm();
            }

            if (total_length < 1e-10) {
                log("警告: Sheet总长度接近零");
                return 0.0;
            }

            // 计算曲率度量 - 值接近1表示直线，值越小表示越弯曲
            double curvature_metric = euclidean_distance / total_length;
            
            log("Sheet曲率度量: " + std::to_string(curvature_metric) + 
                " (欧氏距离: " + std::to_string(euclidean_distance) + 
                ", 总路径长度: " + std::to_string(total_length) + ")");
            
            return curvature_metric; // 0~1之间的值，接近1表示sheet接近直线，值越小表示sheet越弯曲
        }
        catch (const std::exception& e) {
            log("计算sheet曲率度量时发生异常: " + std::string(e.what()));
            return 0.0;
        }
        catch (...) {
            log("计算sheet曲率度量时发生未知异常");
            return 0.0;
        }
    }
	template <typename M>
    inline double sheet_operation<M>::get_sheet_normal_variation(std::vector<E *> sheet)
    {
        if (sheet.empty()) {
            return 0.0;
        }

        try {
            log("Computing sheet normal variation...");

            // Calculate normal variation by measuring angle changes between adjacent face normals
            double total_normal_variation = 0.0;
            int valid_measurements = 0;

            // For each edge in the sheet
            for (auto e : sheet) {
                if (!e) continue;
                
                // Skip edges with insufficient adjacent faces
                if (e->neighbor_fs.size() < 2) continue;

                std::vector<CPoint> face_normals;
                
                // Collect normals for all adjacent faces
                for (int i = 0; i < e->neighbor_fs.size(); ++i) {
                    F* f = mesh->idFaces(e->neighbor_fs[i]);
                    if (!f) continue;
                    
                    face_normals.push_back(f->normal());
                }
                
                // Calculate angle variation between all pairs of faces
                for (size_t i = 0; i < face_normals.size(); ++i) {
                    for (size_t j = i + 1; j < face_normals.size(); ++j) {
                        CPoint n1 = face_normals[i];
                        CPoint n2 = face_normals[j];
                        
                        // Skip invalid normals
                        if (n1.norm() < 1e-8 || n2.norm() < 1e-8) continue;
                        
                        // Normalize vectors
                        n1 = n1 / n1.norm();
                        n2 = n2 / n2.norm();
                        
                        // Compute dot product and clamp to valid range
                        double dot_product = n1 * n2;
                        dot_product = std::max(-1.0, std::min(1.0, dot_product));
                        
                        // Calculate angle in degrees
                        double angle = acos(dot_product) * 180.0 / MPI;
                        
                        // Add to total
                        total_normal_variation += angle;
                        valid_measurements++;
                    }
                }
            }
            
            // Calculate average normal variation
            double average_variation = (valid_measurements > 0) ? 
                (total_normal_variation / valid_measurements) : 0.0;
            
            log("Sheet normal variation: " + std::to_string(average_variation) + 
                " degrees (based on " + std::to_string(valid_measurements) + " measurements)");
            
            return average_variation; //度数，值越大表示面法线变化越剧烈，通常意味着区域曲率越高
        }
        catch (const std::exception& e) {
            log("Exception occurred when computing normal variation: " + std::string(e.what()));
            return 0.0;
        }
        catch (...) {
            log("Unknown exception occurred when computing normal variation");
            return 0.0;
        }
    }
	template <typename M>
    inline double sheet_operation<M>::get_sheet_dihedral_angle_deviation(std::vector<E *> sheet)
    {
        if (sheet.empty()) {
            return 0.0;
        }

        try {
            log("Computing sheet dihedral angle deviation...");
            
            double total_deviation = 0.0;
            int valid_measurements = 0;

            // For each edge in the sheet
            for (auto e : sheet) {
                if (!e) continue;
                
                double ideal_angle = 0.0;
                double actual_angle = e->total_angle();
                
                // Skip edges with invalid angles
                if (!std::isfinite(actual_angle) || actual_angle <= 0) continue;
                
                // Determine ideal angle based on edge properties
                if (e->boundary()) {
                    ideal_angle = e->sharp() ? e->total_angle() : 180.0; // Sharp edges have specific angles, regular boundary edges should be 180°
                } else {
                    ideal_angle = 360.0; // Internal edges should have 360° total angle
                }
                
                // Calculate deviation as absolute difference between actual and ideal
                double deviation = std::abs(actual_angle - ideal_angle);
                
                // For internal edges, take the minimum deviation considering periodicity around 360°
                if (!e->boundary() && deviation > 180.0) {
                    deviation = 360.0 - deviation;
                }
                
                log("  Edge " + std::to_string(e->id()) + 
                    " - actual angle: " + std::to_string(actual_angle) + 
                    ", ideal angle: " + std::to_string(ideal_angle) + 
                    ", deviation: " + std::to_string(deviation));
                
                total_deviation += deviation;
                valid_measurements++;
            }
            
            // Calculate average deviation
            double average_deviation = (valid_measurements > 0) ? 
                (total_deviation / valid_measurements) : 0.0;
            
            log("Sheet dihedral angle deviation: " + std::to_string(average_deviation) + 
                " degrees (based on " + std::to_string(valid_measurements) + " measurements)");
            
            return average_deviation; //度数，值越大表示二面角偏离理想值越远，通常意味着网格畸变更严重
        }
        catch (const std::exception& e) {
            log("Exception occurred when computing dihedral angle deviation: " + std::string(e.what()));
            return 0.0;
        }
        catch (...) {
            log("Unknown exception occurred when computing dihedral angle deviation");
            return 0.0;
        }
    }
	template <typename M>
    inline std::pair<double, double> sheet_operation<M>::get_sheet_adjacent_hex_jacobian(std::vector<E *> sheet)
    {
        if (sheet.empty()) {
            return {0.0, 0.0};
        }

        try {
            log("Computing sheet adjacent hexahedra scaled Jacobian...");
            
            // 用于计算缩放雅可比行列式的辅助函数
            auto calculate_scaled_jacobian = [](H* hex, M* mesh) -> double {
                if (!hex || hex->vs.size() != 8) {
                    return 0.0;
                }
                
                // 获取八个顶点的坐标
                std::vector<CPoint> vertex_positions;
                for (int i = 0; i < 8; ++i) {
                    V* v = mesh->idVertices(hex->vs[i]);
                    if (!v) return 0.0;
                    vertex_positions.push_back(v->position());
                }
                
                // 计算8个角点的雅可比行列式
                double min_jacobian = 1.0;
                
                // 计算每个角点的雅可比行列式并取最小值
                for (int i = 0; i < 8; ++i) {
                    // 获取从当前角点出发的三条边的向量
                    CPoint v0, v1, v2;
                    
                    // 根据六面体角点索引确定相邻顶点
                    // 这里使用简化模型，实际应根据六面体的顶点排列确定正确的边
                    int neighbor_indices[3];
                    switch(i) {
                        case 0: neighbor_indices[0] = 1; neighbor_indices[1] = 3; neighbor_indices[2] = 4; break;
                        case 1: neighbor_indices[0] = 0; neighbor_indices[1] = 2; neighbor_indices[2] = 5; break;
                        case 2: neighbor_indices[0] = 1; neighbor_indices[1] = 3; neighbor_indices[2] = 6; break;
                        case 3: neighbor_indices[0] = 0; neighbor_indices[1] = 2; neighbor_indices[2] = 7; break;
                        case 4: neighbor_indices[0] = 0; neighbor_indices[1] = 5; neighbor_indices[2] = 7; break;
                        case 5: neighbor_indices[0] = 1; neighbor_indices[1] = 4; neighbor_indices[2] = 6; break;
                        case 6: neighbor_indices[0] = 2; neighbor_indices[1] = 5; neighbor_indices[2] = 7; break;
                        case 7: neighbor_indices[0] = 3; neighbor_indices[1] = 4; neighbor_indices[2] = 6; break;
                    }
                    
                    // 计算从当前角点到相邻顶点的向量
                    v0 = vertex_positions[neighbor_indices[0]] - vertex_positions[i];
                    v1 = vertex_positions[neighbor_indices[1]] - vertex_positions[i];
                    v2 = vertex_positions[neighbor_indices[2]] - vertex_positions[i];
                    
                    // 计算雅可比行列式 J = (v0 × v1)·v2
                    CPoint cross_product = v0 ^ v1;
                    double jacobian_det = cross_product * v2;
                    
                    // 计算缩放雅可比行列式
                    double len_v0 = v0.norm();
                    double len_v1 = v1.norm();
                    double len_v2 = v2.norm();
                    
                    if (len_v0 < 1e-10 || len_v1 < 1e-10 || len_v2 < 1e-10) {
                        return 0.0; // 存在零长度边，可能是退化的六面体
                    }
                    
                    double scaled_jacobian = jacobian_det / (len_v0 * len_v1 * len_v2);
                    
                    // 更新最小值
                    min_jacobian = std::min(min_jacobian, scaled_jacobian);
                }
                
                return min_jacobian;
            };
            
            // 收集sheet中所有边的相邻六面体
            std::set<int> adjacent_hex_ids;
            for (auto e : sheet) {
                if (!e) continue;
                
                for (int i = 0; i < e->neighbor_hs.size(); ++i) {
                    adjacent_hex_ids.insert(e->neighbor_hs[i]);
                }
            }
            
            log("Found " + std::to_string(adjacent_hex_ids.size()) + " adjacent hexahedra to evaluate");
            
            // 计算所有相邻六面体的缩放雅可比行列式
            std::vector<double> jacobian_values;
            for (int hex_id : adjacent_hex_ids) {
                H* hex = mesh->idHexs(hex_id);
                if (!hex) continue;
                
                double scaled_jacobian = calculate_scaled_jacobian(hex, mesh);
                if (scaled_jacobian != 0.0) { // 忽略无效值
                    jacobian_values.push_back(scaled_jacobian);
                }
            }
            
            // 计算最小值和平均值
            if (jacobian_values.empty()) {
                log("No valid Jacobian values calculated");
                return {0.0, 0.0};
            }
            
            double min_jacobian = *std::min_element(jacobian_values.begin(), jacobian_values.end());
            double avg_jacobian = 0.0;
            for (double val : jacobian_values) {
                avg_jacobian += val;
            }
            avg_jacobian /= jacobian_values.size();
            
            log("Sheet adjacent hexahedra - Minimum scaled Jacobian: " + std::to_string(min_jacobian) + 
                ", Average scaled Jacobian: " + std::to_string(avg_jacobian) + 
                " (based on " + std::to_string(jacobian_values.size()) + " hexahedra)");
            
            return {min_jacobian, avg_jacobian};
        }
        catch (const std::exception& e) {
            log("Exception occurred when computing scaled Jacobian: " + std::string(e.what()));
            return {0.0, 0.0};
        }
        catch (...) {
            log("Unknown exception occurred when computing scaled Jacobian");
            return {0.0, 0.0};
        }
    }

	template <typename M>
	inline int sheet_operation<M>::get_sheet_valence_along_path(const std::vector<E *> &sheet)
	{
		// 处理空 sheet 输入
		if (sheet.empty()) {
			return 0;
		}
	
		// --- 1. 识别端点顶点 ---
		std::map<int, int> vertex_counts;
		std::set<int> sheet_vertex_ids;
		std::set<int> endpoint_ids;
		std::set<int> sheet_edge_ids;
		int current_sheet_id = -1; // 用于比较其他边的 sheet ID
	
		for (E* e : sheet) {
			// 跳过无效边（空指针或顶点数不为2）
			if (!e || e->vs.size() != 2) {
				continue;
			}
			sheet_edge_ids.insert(e->id());
			// 获取 sheet ID (从第一个有效边获取)
			if(current_sheet_id == -1) current_sheet_id = e->sheet();
	
			for (int vid : e->vs) {
				sheet_vertex_ids.insert(vid);
				vertex_counts[vid]++;
			}
		}
	
		// 如果 sheet 中全是无效边，直接返回 0
		if (sheet_vertex_ids.empty()) {
			return 0;
		}
		// 如果无法确定 sheet ID (例如 sheet 中所有边 sheet() <= 0)，后面 is_other_sheet 会失效
		// 但基于 is_feature 的计数仍然有效
	
		for (int vid : sheet_vertex_ids) {
			// 使用 find 避免在 map 中自动创建不存在的 key
			auto it = vertex_counts.find(vid);
			// 只有在 map 中找到且计数为 1 的才是端点
			if (it != vertex_counts.end() && it->second == 1) {
				endpoint_ids.insert(vid);
			}
		}
	
		// --- 2. 统计连接数 ---
		int valence_count = 0;
		std::set<int> counted_connecting_edges; // 避免重复计数
	
		for (int vid : sheet_vertex_ids) {
			// 跳过端点
			if (endpoint_ids.count(vid)) {
				continue;
			}
	
			V* v = mesh->idVertices(vid);
			// 跳过无效顶点指针
			if (!v) {
				continue;
			}
	
			// 遍历顶点的邻接边
			for (int neighbor_eid : v->neighbor_es) {
				// 跳过属于当前 sheet 的边
				if (sheet_edge_ids.count(neighbor_eid)) {
					continue;
				}
				// 跳过已计数的连接边
				if (counted_connecting_edges.count(neighbor_eid)) {
					continue;
				}
	
				E* neighbor_e = mesh->idEdges(neighbor_eid);
				// 跳过无效邻接边指针
				if (!neighbor_e){
					continue;
				}
	
				// 检查邻接边是否是特征边或属于其他 sheet
				bool is_feature = neighbor_e->sharp() > 0;
				bool is_other_sheet = (current_sheet_id > 0) && // 仅当当前 sheet ID 有效时才比较
									(neighbor_e->sheet() > 0) &&
									(neighbor_e->sheet() != current_sheet_id);
	
				if (is_feature || is_other_sheet) {
					valence_count++;
					counted_connecting_edges.insert(neighbor_eid); // 标记为已计数
				}
			}
		} // 结束遍历 sheet 顶点
	
		return valence_count;
	}
	template <typename M>
	inline double sheet_operation<M>::get_distance_to_boundary(const std::vector<E *> &sheet)
	{
		// 检查 mesh 指针和空 sheet
		if (!mesh || sheet.empty()) {
			// 对于无效输入，返回一个表示“无限远”或错误的值
			return std::numeric_limits<double>::max();
		}
	
		std::set<int> sheet_vertex_ids;
	
		// 1. 收集 sheet 顶点，并检查 sheet 是否本身在边界上
		for (E* e : sheet) {
			// 跳过无效边 (空指针或顶点数不为2)
			if (!e || e->vs.size() != 2) {
				continue;
			}
			// 检查边的 boundary 标志
			if (e->boundary()) {
				return 0.0; // 如果 sheet 包含边界边，距离为 0
			}
			// 添加顶点到集合中
			sheet_vertex_ids.insert(e->vs[0]);
			sheet_vertex_ids.insert(e->vs[1]);
		}
	
		// 如果 sheet 为空或只包含无效边
		if (sheet_vertex_ids.empty()){
			 return std::numeric_limits<double>::max();
		}
	
	
		// 2. 收集网格中所有的边界顶点
		std::vector<V*> boundary_vertices;
		// 使用 map 迭代器访问顶点
		for (auto const& [id, v_ptr] : mesh->m_map_vertices) {
			 V* v = v_ptr;
			 // 跳过空指针顶点，并检查 boundary 标志
			 if (v && v->boundary()) {
				 boundary_vertices.push_back(v);
			 }
		}
	
		// 如果网格中没有边界顶点，返回“无限远”
		if (boundary_vertices.empty()) {
			return std::numeric_limits<double>::max();
		}
	
		// 3. 计算 sheet 顶点到边界顶点的最小距离
		double min_distance_sq = std::numeric_limits<double>::max();
	
		for (int sv_id : sheet_vertex_ids) {
			V* sv = mesh->idVertices(sv_id);
			// 跳过无效的 sheet 顶点指针
			if (!sv){
				 continue;
			}
			CPoint sp = sv->position();
	
			for (V* bv : boundary_vertices) {
				// 跳过无效的边界顶点指针
				if (!bv){
					continue;
				}
				CPoint bp = bv->position();
				CPoint diff = sp - bp;
				double dist_sq = diff * diff; // 使用平方距离避免开方，提高效率
				min_distance_sq = std::min(min_distance_sq, dist_sq);
			}
		}
	
		// 如果未能计算出任何有效距离（例如所有 sheet 顶点都无效）
		if (min_distance_sq == std::numeric_limits<double>::max()) {
			 return std::numeric_limits<double>::max();
		}
	
		// 返回实际的最小距离
		return std::sqrt(min_distance_sq);
	}

	template <typename M>
inline double sheet_operation<M>::get_distance_to_feature(const std::vector<E *> &sheet)
{
    // 函数开头不再包含检查 sheet 内部 sharp 边的循环

    // log("--- Entering MODIFIED get_distance_to_feature (assumes sheet has no internal features) ---");

    // 1. 基本输入检查
    if (!mesh || sheet.empty()) {
        // log("get_distance_to_feature: Error - Mesh pointer is null or input sheet is empty.");
        return std::numeric_limits<double>::max();
    }

    // 2. 收集 sheet 自身的顶点 ID
    std::set<int> sheet_vertex_ids;
    std::set<int> sheet_edge_ids; // 仍然需要收集 edge_ids 用于下一步排除
    for (E* e : sheet) {
        if (!e || e->vs.size() != 2) continue;
        sheet_edge_ids.insert(e->id());
        sheet_vertex_ids.insert(e->vs[0]);
        sheet_vertex_ids.insert(e->vs[1]);
        // !!! 移除了这里的 if (e->sharp() > 0) return 0.0; 检查 !!!
    }

    // 如果 sheet 无有效顶点
    if (sheet_vertex_ids.empty()){
         // log("get_distance_to_feature: Error - Sheet contains no valid vertices.");
         return std::numeric_limits<double>::max();
    }

    // 3. 收集网格中 *其他* sharp 边的顶点 ID
    // log("  Step: Collecting vertices from *other* sharp edges...");
    std::set<int> feature_vertex_ids;
    for(auto const& [id, e_ptr] : mesh->m_map_edges) {
        E* e_all = e_ptr;
        if(!e_all || e_all->vs.size() != 2) continue;

        if (e_all->sharp() > 0) {
            // 确保这条 sharp 边 不属于 输入的 sheet
            if (sheet_edge_ids.find(e_all->id()) == sheet_edge_ids.end()) {
                feature_vertex_ids.insert(e_all->vs[0]);
                feature_vertex_ids.insert(e_all->vs[1]);
            }
        }
    }

    // 如果网格中没有其他 feature edges
    if (feature_vertex_ids.empty()) {
        // log("  No other feature edges found in the mesh.");
        return std::numeric_limits<double>::max(); // 返回极大值表示距离无限远
    }

    // 4. 计算 sheet 顶点 到 其他特征顶点 的最小非零距离
    // log("  Step: Calculating minimum non-zero distance...");
    double min_distance_sq = std::numeric_limits<double>::max();
    const double ZERO_THRESHOLD_SQ = 1e-12; // 用于判断距离是否接近零的阈值（平方）

    for (int sv_id : sheet_vertex_ids) {
        V* sv = mesh->idVertices(sv_id);
        if (!sv) continue; // 跳过无效顶点
        CPoint sp = sv->position();

        for (int fv_id : feature_vertex_ids) {
            // 通常不需要比较 sv_id == fv_id，因为 feature_vertex_ids 是来自其他边的
            V* fv = mesh->idVertices(fv_id);
            if (!fv) continue; // 跳过无效顶点
            CPoint fp = fv->position();
            CPoint diff = sp - fp;
            double dist_sq = diff * diff;

            // 只考虑大于阈值的距离
            if (dist_sq > ZERO_THRESHOLD_SQ) {
                min_distance_sq = std::min(min_distance_sq, dist_sq);
            }
            // (可选) 如果需要记录发现零距离的情况可以在此添加日志
            // else { log("    Ignoring near-zero distance between sheet vertex " + ...); }
        }
    }

    // 处理没有找到有效非零距离的情况
    if (min_distance_sq == std::numeric_limits<double>::max()) {
         // log("  Warning: Could not calculate any valid *non-zero* distances.");
         // 根据策略决定返回值，这里返回极大值
         return std::numeric_limits<double>::max();
    }

    // 5. 返回最终结果
    double min_distance = std::sqrt(min_distance_sq);
    // log("get_distance_to_feature: Calculation complete. Min non-zero distance = " + std::to_string(min_distance));
    // log("--- Exiting MODIFIED get_distance_to_feature ---");
    return min_distance;
	}
}