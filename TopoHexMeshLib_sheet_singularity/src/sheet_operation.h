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
        static std::ofstream ofs;
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

        void reset_sheet_attributes() {
            for (M::MEIterator eite(mesh); !eite.end(); eite++) {
                E* be = *eite;
                be->sheet() = 0;
            }
        }
		std::vector<E*> get_one_sheet(E* e);
		int get_mesh_sheet_number();
		void collapse_one_sheet2(std::vector<E*> sheet);

		double vector_angle(CPoint a, CPoint b);
		void edge_angle(E* e);
		void edge_ideal_degree(E* e);
		void compute_edge_energy();
		std::vector<std::vector<E*>> get_sheet_parallel_edges(std::vector<E*>sheet);
		double predict_sheet_collapse_energy(std::vector<E*> sheet);
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
		eQueue.push(e);
		e->mark() = true;
		while (!eQueue.empty())
		{
			E* be = eQueue.front();
			sheet.push_back(be);
			eQueue.pop();
			for (int ehIndex = 0; ehIndex < be->neighbor_hs.size(); ehIndex++)
			{
				H* h = mesh->idHexs(be->neighbor_hs[ehIndex]);
				std::vector<E*> parallel_es = mesh->e_parallel_e_in_hex(h, be);
				for (int peIndex = 0; peIndex < parallel_es.size(); peIndex++)
				{
					E* pe = parallel_es[peIndex];
					if (pe->mark() == false)
					{
						eQueue.push(pe);
						pe->mark() = true;
					}
				}
			}
		}
		for (int eIndex = 0; eIndex < sheet.size(); eIndex++)
		{
			E* e = sheet[eIndex];
			e->mark() = false;
		}
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
				log("Warning: Input face list is empty");
				return;
			}
			
			// 第一次遍历：标记边界并修正法线
			for (int fIndex = 0; fIndex < fs.size(); fIndex++)
			{
				try {
					F* f = fs[fIndex];
					
					// 检查面是否有效
					if (f == nullptr) {
						log("Error: Detected null face pointer, index: " + std::to_string(fIndex));
						continue;
					}
					
					log("Processing face ID: " + std::to_string(f->id()));
					
					// if (f->id() == 465) {
					// 	log("Special debug for face 465:");
					// 	log("  neighbor_hs count: " + std::to_string(f->neighbor_hs.size()));
					// 	log("  vertices count: " + std::to_string(f->vs.size()));
					// 	log("  edges count: " + std::to_string(f->es.size()));
					// }
					
					// 检查邻接六面体列表是否有效
					if (f->neighbor_hs.empty()) {
						log("Warning: Face " + std::to_string(f->id()) + " has no adjacent hexahedra");
					}
					
					//mark boundary
					if (f->neighbor_hs.size() == 1)
					{
						f->boundary() = true;
						log("Face " + std::to_string(f->id()) + " marked as boundary face");
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
								log("Error: Face " + std::to_string(f->id()) + " has invalid adjacent hexahedron ID: " + std::to_string(f->neighbor_hs[fhIndex]));
								continue;
							}
							
							log("About to get hexahedron with ID: " + std::to_string(f->neighbor_hs[fhIndex]));
							H* h = mesh->idHexs(f->neighbor_hs[fhIndex]);
							log("Successfully got hexahedron with ID: " + std::to_string(h->id()));
							
							// 检查六面体是否存在
							if (h == nullptr) {
								log("Error: Unable to get adjacent hexahedron ID: " + std::to_string(f->neighbor_hs[fhIndex]) + " for face " + std::to_string(f->id()));
								continue;
							}
							
							try {
								mesh->revise_hex_face_normal(h);
							}
							catch (const std::exception& e) {
								log("Error while correcting hexahedron face normal, hex ID: " + std::to_string(h->id()) + ", exception: " + std::string(e.what()));
							}
							catch (...) {
								log("Unknown error while correcting hexahedron face normal, hex ID: " + std::to_string(h->id()));
							}
						}
						catch (const std::exception& e) {
							log("Exception in hexahedron processing: " + std::string(e.what()) + ", face ID: " + std::to_string(f->id()) + ", hex index: " + std::to_string(fhIndex));
						}
						catch (...) {
							log("Unknown exception in hexahedron processing, face ID: " + std::to_string(f->id()) + ", hex index: " + std::to_string(fhIndex));
						}
					}
				}
				catch (const std::exception& e) {
					log("Exception in first traversal at face index " + std::to_string(fIndex) + ": " + std::string(e.what()));
				}
				catch (...) {
					log("Unknown exception in first traversal at face index " + std::to_string(fIndex));
				}
			}

			log("First traversal completed, now marking vertices and edges");
			
			// 第二次遍历：标记顶点和边
			for (int fIndex = 0; fIndex < fs.size(); fIndex++)
			{
				try {
					F* f = fs[fIndex];
					
					// 检查面是否有效
					if (f == nullptr) {
						log("Skipping null face, index: " + std::to_string(fIndex));
						continue;
					}
					
					log("Processing face for vertices/edges ID: " + std::to_string(f->id()));
					
					// 标记面顶点
					try {
						log("Marking " + std::to_string(f->vs.size()) + " vertices for face " + std::to_string(f->id()));
						for (int fvIndex = 0; fvIndex < f->vs.size(); fvIndex++)
						{
							try {
								// 检查顶点ID是否有效
								if (f->vs[fvIndex] < 0 || f->vs[fvIndex] >= mesh->maxVid()) {
									log("Error: Face " + std::to_string(f->id()) + " has invalid vertex ID: " + std::to_string(f->vs[fvIndex]));
									continue;
								}
								
								V* fv = mesh->idVertices(f->vs[fvIndex]);
								
								// 检查顶点是否存在
								if (fv == nullptr) {
									log("Error: Unable to get vertex ID: " + std::to_string(f->vs[fvIndex]) + " for face " + std::to_string(f->id()));
									continue;
								}
								
								if (f->boundary())
								{
									fv->boundary() = true;
									log("Vertex " + std::to_string(fv->id()) + " marked as boundary vertex");
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
												log("Warning: Vertex " + std::to_string(fv->id()) + " has invalid adjacent face ID: " + std::to_string(fv->neighbor_fs[vfIndex]));
												continue;
											}
											
											F* fvf = mesh->idFaces(fv->neighbor_fs[vfIndex]);
											
											// 检查面是否存在
											if (fvf == nullptr) {
												log("Warning: Unable to get adjacent face ID: " + std::to_string(fv->neighbor_fs[vfIndex]) + " for vertex " + std::to_string(fv->id()));
												continue;
											}
											
											if (fvf->boundary())
											{
												is_boundary = true;
												log("Vertex " + std::to_string(fv->id()) + " marked as boundary vertex (via adjacent face)");
												break;
											}
										}
										catch (const std::exception& e) {
											log("Exception in vertex-face processing: " + std::string(e.what()) + ", vertex ID: " + std::to_string(fv->id()) + ", face index: " + std::to_string(vfIndex));
										}
										catch (...) {
											log("Unknown exception in vertex-face processing, vertex ID: " + std::to_string(fv->id()) + ", face index: " + std::to_string(vfIndex));
										}
									}
									fv->boundary() = is_boundary;
								}
							}
							catch (const std::exception& e) {
								log("Exception processing vertex at index " + std::to_string(fvIndex) + " for face " + std::to_string(f->id()) + ": " + std::string(e.what()));
							}
							catch (...) {
								log("Unknown exception processing vertex at index " + std::to_string(fvIndex) + " for face " + std::to_string(f->id()));
							}
						}
					}
					catch (const std::exception& e) {
						log("Exception in vertices processing for face " + std::to_string(f->id()) + ": " + std::string(e.what()));
					}
					catch (...) {
						log("Unknown exception in vertices processing for face " + std::to_string(f->id()));
					}
					
					// 标记面边
					try {
						log("Marking " + std::to_string(f->es.size()) + " edges for face " + std::to_string(f->id()));
						for (int feIndex = 0; feIndex < f->es.size(); feIndex++)
						{
							try {
								// 检查边ID是否有效
								if (f->es[feIndex] < 0 || f->es[feIndex] >= mesh->maxEid()) {
									log("Error: Face " + std::to_string(f->id()) + " has invalid edge ID: " + std::to_string(f->es[feIndex]));
									continue;
								}
								
								E* fe = mesh->idEdges(f->es[feIndex]);
								
								// 检查边是否存在
								if (fe == nullptr) {
									log("Error: Unable to get edge ID: " + std::to_string(f->es[feIndex]) + " for face " + std::to_string(f->id()));
									continue;
								}
								
								// 边界标记
								if (f->boundary())
								{
									fe->boundary() = true;
									log("Edge " + std::to_string(fe->id()) + " marked as boundary edge");
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
												log("Warning: Edge " + std::to_string(fe->id()) + " has invalid adjacent face ID: " + std::to_string(fe->neighbor_fs[efIndex]));
												continue;
											}
											
											F* fef = mesh->idFaces(fe->neighbor_fs[efIndex]);
											
											// 检查面是否存在
											if (fef == nullptr) {
												log("Warning: Unable to get adjacent face ID: " + std::to_string(fe->neighbor_fs[efIndex]) + " for edge " + std::to_string(fe->id()));
												continue;
											}
											
											if (fef->boundary())
											{
												is_boundary = true;
												log("Edge " + std::to_string(fe->id()) + " marked as boundary edge (via adjacent face)");
												break;
											}
										}
										catch (const std::exception& e) {
											log("Exception in edge-face processing: " + std::string(e.what()) + ", edge ID: " + std::to_string(fe->id()) + ", face index: " + std::to_string(efIndex));
										}
										catch (...) {
											log("Unknown exception in edge-face processing, edge ID: " + std::to_string(fe->id()) + ", face index: " + std::to_string(efIndex));
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
											log("Boundary edge " + std::to_string(fe->id()) + " marked as singular edge, adjacent hexahedra count: " + std::to_string(fe->neighbor_hs.size()));
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
											log("Non-boundary edge " + std::to_string(fe->id()) + " marked as singular edge, adjacent hexahedra count: " + std::to_string(fe->neighbor_hs.size()));
										}
										else
										{
											fe->singularity() = false;
										}
									}
								}
								catch (const std::exception& e) {
									log("Exception in edge singularity marking: " + std::string(e.what()) + ", edge ID: " + std::to_string(fe->id()));
								}
								catch (...) {
									log("Unknown exception in edge singularity marking, edge ID: " + std::to_string(fe->id()));
								}
								
								try {
									// 简化能量计算
									if (fe->boundary())
									{
										if (fe->sharp())
										{
											// 检查总角度是否合理
											if (fe->total_angle() < 0 || fe->total_angle() > 360) {
												log("Warning: Edge " + std::to_string(fe->id()) + " has abnormal total angle value: " + std::to_string(fe->total_angle()));
												// 使用默认值避免计算错误
												fe->ideal_degree() = 2;
											}
											else {
												int ideal_degree = round(fe->total_angle() / 90.0);
												ideal_degree = ideal_degree == 0 ? 1 : ideal_degree;
												fe->ideal_degree() = ideal_degree;
												log("Sharp boundary edge " + std::to_string(fe->id()) + " ideal degree set to: " + std::to_string(ideal_degree) + " (total angle: " + std::to_string(fe->total_angle()) + ")");
											}
										}
										else
										{
											fe->ideal_degree() = 2;
											log("Regular boundary edge " + std::to_string(fe->id()) + " ideal degree set to: 2");
										}
									}
									else
									{
										fe->ideal_degree() = 4;
										log("Non-boundary edge " + std::to_string(fe->id()) + " ideal degree set to: 4");
									}
									
									// 计算简化能量
									fe->sim_energy() = (int)fe->ideal_degree() - (int)fe->neighbor_hs.size();
									log("Edge " + std::to_string(fe->id()) + " simplification energy: " + std::to_string(fe->sim_energy()) + 
										" (ideal degree: " + std::to_string(fe->ideal_degree()) + ", actual adjacent hexahedra: " + std::to_string(fe->neighbor_hs.size()) + ")");
								}
								catch (const std::exception& e) {
									log("Exception in edge energy calculation: " + std::string(e.what()) + ", edge ID: " + std::to_string(fe->id()));
								}
								catch (...) {
									log("Unknown exception in edge energy calculation, edge ID: " + std::to_string(fe->id()));
								}
							}
							catch (const std::exception& e) {
								log("Exception processing edge at index " + std::to_string(feIndex) + " for face " + std::to_string(f->id()) + ": " + std::string(e.what()));
							}
							catch (...) {
								log("Unknown exception processing edge at index " + std::to_string(feIndex) + " for face " + std::to_string(f->id()));
							}
						}
					}
					catch (const std::exception& e) {
						log("Exception in edges processing for face " + std::to_string(f->id()) + ": " + std::string(e.what()));
					}
					catch (...) {
						log("Unknown exception in edges processing for face " + std::to_string(f->id()));
					}
				}
				catch (const std::exception& e) {
					log("Exception in second traversal at face index " + std::to_string(fIndex) + ": " + std::string(e.what()));
				}
				catch (...) {
					log("Unknown exception in second traversal at face index " + std::to_string(fIndex));
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
		log("Starting to process sheet");
		/*get all the hexs and the include elements*/
		std::set<H*> delete_hs;//get all hex in sheet
		std::set<F*> delete_fs;//delete faces
		std::set<E*> delete_es;//delete edges
		std::set<V*> delete_vs;//delete vs
		for (int eIndex = 0; eIndex < sheet.size(); eIndex++)
		{
			E* e = sheet[eIndex];
			e->sheetId() = 1;
			log("Processing edge id: " + std::to_string(e->id()));
			//find all hex
			for (int ehIndex = 0; ehIndex < e->neighbor_hs.size(); ehIndex++)
			{
				H* h = mesh->idHexs(e->neighbor_hs[ehIndex]);
				if (h->is_delete()) continue;

				h->is_delete() = true;
				h->sheetId() = 1;
				delete_hs.insert(h);

				/*mark the delete face and the elements of it*/
				for (int hvIndex = 0; hvIndex < h->vs.size(); hvIndex++)
				{
					V* hv = mesh->idVertices(h->vs[hvIndex]);
					hv->is_delete() = true;
					delete_vs.insert(hv);
				}
				for (int heIndex = 0; heIndex < h->es.size(); heIndex++)
				{
					E* he = mesh->idEdges(h->es[heIndex]);
					he->is_delete() = true;
					delete_es.insert(he);
				}
				for (int hfIndex = 0; hfIndex < h->fs.size(); hfIndex++)
				{
					F* hf = mesh->idFaces(h->fs[hfIndex]);
					hf->is_delete() = true;
					delete_fs.insert(hf);
				}
			}
		}
		/*exit if there is self-intersection*/
		for (std::set<H*>::iterator hite = delete_hs.begin(); hite != delete_hs.end(); hite++)
		{
			H* h = *hite;
			int count = 0;
			for (int heIndex = 0; heIndex < h->es.size(); heIndex++)
			{
				E* e = mesh->idEdges(h->es[heIndex]);
				if (e->sheetId())
				{
					count++;
				}
			}
			if (count > 4)
			{
				log("Self-intersection detected at hex id: " + std::to_string(h->id()));
				std::cout << "The current sheet has self-intersection" << std::endl;
				return;
			}
		}

		/*new elements*/
		std::vector<V*> newvs;
		std::vector<E*> newes;
		std::vector<F*> newfs;
		/*add new vertices*/
		log("Starting new vertex creation");
		for (int eIndex = 0; eIndex < sheet.size(); eIndex++)
		{
			E* e = sheet[eIndex];
			V* newv = construct_new_vertex();
			log("Created new vertex id: " + std::to_string(newv->id()));
			newvs.push_back(newv);

			V* ev1 = mesh->idVertices(e->vs[0]);
			V* ev2 = mesh->idVertices(e->vs[1]);
			vertex_position_policy(ev1, ev2, newv);

			std::vector<V*> evs = { ev1,ev2 };
			e->setNewv(newv);
			//change relations
			for (int evIndex = 0; evIndex < 2; evIndex++)
			{
				V* oldv = evs[evIndex];
				//change edge relations
				for (int veIndex = 0; veIndex < oldv->neighbor_es.size(); veIndex++)
				{
					E* ve = mesh->idEdges(oldv->neighbor_es[veIndex]);
					if (ve->is_delete()) continue;
					int vevIndex = ve->vertexIndex(oldv->id());
					ve->vs[vevIndex] = newv->id();
					newv->neighbor_es.push_back(ve->id());
				}
				//change face relations
				for (int vfIndex = 0; vfIndex < oldv->neighbor_fs.size(); vfIndex++)
				{
					F* vf = mesh->idFaces(oldv->neighbor_fs[vfIndex]);
					if (vf->is_delete()) continue;
					int vfvIndex = vf->vertexIndex(oldv->id());
					if (vfvIndex == -1)
					{
						continue;
					}
					vf->vs[vfvIndex] = newv->id();
					newv->neighbor_fs.push_back(vf->id());
				}
				//change hex relations
				for (int vhIndex = 0; vhIndex < oldv->neighbor_hs.size(); vhIndex++)
				{
					H* vh = mesh->idHexs(oldv->neighbor_hs[vhIndex]);
					if (vh->is_delete()) continue;
					int vhvIndex = vh->vertexIndex(oldv->id());
					if (vhvIndex == -1)
					{
						continue;
					}
					vh->vs[vhvIndex] = newv->id();
					newv->neighbor_hs.push_back(vh->id());
				}
			}
		}
		log("Vertex creation completed");
		/*create new edges and faces*/
		for (std::set<H*>::iterator hite = delete_hs.begin(); hite != delete_hs.end(); hite++)
		{
			H* h = *hite;
			log("Processing hex (new edge/face creation) id: " + std::to_string(h->id()));
			//find the faces without mark edges
			std::vector<F*> fs;
			for (int hfIndex = 0; hfIndex < h->fs.size(); hfIndex++)
			{
				F* hf = mesh->idFaces(h->fs[hfIndex]);
				bool without_mark_edge = true;
				for (int feIndex = 0; feIndex < hf->es.size(); feIndex++)
				{
					E* fe = mesh->idEdges(hf->es[feIndex]);
					if (fe->sheetId())
					{
						without_mark_edge = false;
						break;
					}
				}
				if (without_mark_edge)
				{
					F* parallel_f = mesh->f_parallel_f_in_hex(h, hf);
					parallel_f->sheetId() = 1;
					fs = { hf,parallel_f };
				}
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

		/*delete the elements*/
		log("Starting deletion operations");
		//delete hex
		std::list<H*>::iterator hite = mesh->hs.begin();
		while (hite != mesh->hs.end())
		{
			H* h = *hite;
			if (h == NULL)
			{
				hite = mesh->hs.erase(hite);
			}
			else if (h->is_delete())
			{
				mesh->m_map_hexs.erase(h->id());
				hite = mesh->hs.erase(hite);
			}
			else
			{
				++hite;
			}
		}

		//delete face
		std::list<F*>::iterator fite = mesh->fs.begin();
		while (fite != mesh->fs.end())
		{
			F* f = *fite;
			if (f == NULL)
			{
				fite = mesh->fs.erase(fite);
			}
			else if (f->is_delete())
			{
				mesh->m_map_faces.erase(f->id());
				fite = mesh->fs.erase(fite);
			}
			else
			{
				++fite;
			}
		}
		//delete edge
		std::list<E*>::iterator eite = mesh->es.begin();
		while (eite != mesh->es.end())
		{
			E* e = *eite;
			if (e == NULL)
			{
				eite = mesh->es.erase(eite);
			}
			else if (e->is_delete())
			{
				mesh->m_map_edges.erase(e->id());
				eite = mesh->es.erase(eite);
			}
			else
			{
				++eite;
			}
		}
		//delete vertex
		std::list<V*>::iterator vite = mesh->vs.begin();
		while (vite != mesh->vs.end())
		{
			V* v = *vite;
			if (v == NULL)
			{
				vite = mesh->vs.erase(vite);
			}
			else if (v->is_delete())
			{
				mesh->m_map_vertices.erase(v->id());
				vite = mesh->vs.erase(vite);
			}
			else
			{
				++vite;
			}
		}
		log("Deletion operations completed");

		log("Starting marking element attributes");
		try {
			log("New faces count: " + std::to_string(newfs.size()));
			mark_elements_attribute_collapse(newfs);
			log("Element attribute marking completed");
			
			log("Starting marking singularities");
			mesh->singularities = mesh->mark_singularity();
			log("Singularities marking completed: " + std::to_string(mesh->singularities.size()) + " singularities found");
		}
		catch (const std::exception& e) {
			log("Exception occurred: " + std::string(e.what()));
		}
		catch (...) {
			log("Unknown exception occurred during post-processing");
		}
		
		log("Sheet collapse execution completed");

		for (int i = 0; i < newfs.size(); i++) {
			F* f = newfs[i];
			log("Checking face " + std::to_string(f->id()));
			log("  vs size: " + std::to_string(f->vs.size()));
			log("  es size: " + std::to_string(f->es.size()));
			log("  neighbor_hs size: " + std::to_string(f->neighbor_hs.size()));
			// 检查顶点是否都有效
			for (int j = 0; j < f->vs.size(); j++) {
				if (f->vs[j] < 0 || f->vs[j] >= mesh->maxVid()) {
					log("  Warning: Invalid vertex ID: " + std::to_string(f->vs[j]));
				}
			}
			// 检查边是否都有效
			for (int j = 0; j < f->es.size(); j++) {
				if (f->es[j] < 0 || f->es[j] >= mesh->maxEid()) {
					log("  Warning: Invalid edge ID: " + std::to_string(f->es[j]));
				}
			}
		}
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
		try {
			log("Computing angle for edge ID: " + std::to_string(e->id()));
			
			// 如果已经有值，检查是否为有效值（避免使用未初始化的内存垃圾值）
			if (e->total_angle() != 0)
			{
				// 添加有效性检查，确保 total_angle 的值在合理范围内
				if (std::isfinite(e->total_angle()) && e->total_angle() > 0 && e->total_angle() <= 1080) {
					log("Edge " + std::to_string(e->id()) + " already has valid total angle: " + std::to_string(e->total_angle()));
					return;
				} else {
					log("Edge " + std::to_string(e->id()) + " has invalid total angle: " + std::to_string(e->total_angle()) + ", recalculating...");
					// 重置异常值
					e->total_angle() = 0;
				}
			}
			
			// 非边界边处理
			if (!e->boundary())
			{
				e->total_angle() = 360.0;
				e->ave_angle() = 90.0; // 一般为90度
				log("Non-boundary edge " + std::to_string(e->id()) + " setting default total angle: 360.0");
				return;
			}
			
			// 边界边需要计算
			double total_angle = 0;
			int valid_angles_count = 0;
			
			// 遍历相邻六面体
			for (M::EHIterator ehite(mesh, e); !ehite.end(); ehite++)
			{
				H* eh = *ehite;
				if (!eh) {
					log("Warning: Null hexahedron found for edge " + std::to_string(e->id()));
					continue;
				}
				
				try {
					// 获取相邻面
					std::vector<F*> ehfs = mesh->e_adj_f_in_hex(eh, e);
					if (ehfs.size() < 2) {
						log("Warning: Edge " + std::to_string(e->id()) + " has insufficient adjacent faces in hex " + std::to_string(eh->id()) + ", faces: " + std::to_string(ehfs.size()));
						continue;
					}
					
					// 获取法向量
					std::vector<CPoint> fn;
					bool valid_normals = true;
					
					for (int fIndex = 0; fIndex < ehfs.size(); fIndex++)
					{
						F* f = ehfs[fIndex];
						if (!f) {
							log("Warning: Null face found for edge " + std::to_string(e->id()) + " in hex " + std::to_string(eh->id()));
							valid_normals = false;
							break;
						}
						
						CPoint n = f->normal();
						
						// 检查法向量是否有效
						if (fabs(n.norm()) < 1e-10) {
							log("Warning: Face " + std::to_string(f->id()) + " has invalid normal (zero length)");
							valid_normals = false;
							break;
						}
						
						if (f->neighbor_hs.empty()) {
							log("Warning: Face " + std::to_string(f->id()) + " has empty neighbor_hs array");
							valid_normals = false;
							break;
						}
						
						if (f->neighbor_hs[0] == eh->id())
						{
							fn.push_back(n);
						}
						else
						{
							fn.push_back(-n);
						}
					}
					
					if (!valid_normals || fn.size() < 2) {
						log("Warning: Skipping angle calculation for edge " + std::to_string(e->id()) + " in hex " + std::to_string(eh->id()) + " due to invalid normals");
						continue;
					}
					
					// 计算角度
					double angle_between = vector_angle(fn[0], fn[1]);
					
					// 检查角度是否有效
					if (!std::isfinite(angle_between) || angle_between < 0 || angle_between > 180) {
						log("Warning: Invalid angle between normals: " + std::to_string(angle_between) + " for edge " + std::to_string(e->id()));
						continue;
					}
					
					double temp_angle = 180 - angle_between;
					
					// 范围检查
					if (temp_angle < 0) temp_angle = 0;
					if (temp_angle > 180) temp_angle = 180;
					
					eh->edge_angle(e->id()) = temp_angle;
					total_angle += temp_angle;
					valid_angles_count++;
					
					log("  Hex " + std::to_string(eh->id()) + " contributes angle: " + std::to_string(temp_angle) + " to edge " + std::to_string(e->id()));
				}
				catch (const std::exception& ex) {
					log("Exception in edge angle calculation for hex " + std::to_string(eh->id()) + ": " + std::string(ex.what()));
				}
				catch (...) {
					log("Unknown exception in edge angle calculation for hex " + std::to_string(eh->id()));
				}
			}
			
			// 检查是否有有效角度
			if (valid_angles_count == 0) {
				log("Warning: No valid angles calculated for edge " + std::to_string(e->id()) + ", setting default angle");
				e->total_angle() = e->boundary() ? 180.0 : 360.0;
				e->ave_angle() = e->boundary() ? 90.0 : 90.0;
				return;
			}
			
			// 验证最终的 total_angle
			if (!std::isfinite(total_angle) || total_angle < 0 || total_angle > 720) {
				log("Error: Invalid total angle calculated: " + std::to_string(total_angle) + " for edge " + std::to_string(e->id()) + ", resetting to default");
				total_angle = e->boundary() ? 180.0 : 360.0;
			}
			
			e->total_angle() = total_angle;
			e->ave_angle() = valid_angles_count > 0 ? (total_angle / valid_angles_count) : 90.0;
			
			log("Edge " + std::to_string(e->id()) + " final total angle: " + std::to_string(e->total_angle()) + 
				", average angle: " + std::to_string(e->ave_angle()) + 
				", based on " + std::to_string(valid_angles_count) + " valid measurements");
		}
		catch (const std::exception& ex) {
			log("Major exception in edge_angle for edge " + std::to_string(e->id()) + ": " + std::string(ex.what()));
			// 设置安全的默认值
			e->total_angle() = e->boundary() ? 180.0 : 360.0;
			e->ave_angle() = 90.0;
		}
		catch (...) {
			log("Unknown major exception in edge_angle for edge " + std::to_string(e->id()));
			// 设置安全的默认值
			e->total_angle() = e->boundary() ? 180.0 : 360.0;
			e->ave_angle() = 90.0;
		}
	}

	template<typename M>
	void sheet_operation<M>::edge_ideal_degree(E* e)
	{
		try {
			log("Computing ideal degree for edge ID: " + std::to_string(e->id()));
			
			// 如果已经有值，检查是否为有效值（避免使用未初始化的内存垃圾值）
			if (e->ideal_degree() != 0)
			{
				// 添加有效性检查，确保 ideal_degree 的值在合理范围内
				if (e->ideal_degree() >= 1 && e->ideal_degree() <= 12) {
					log("Edge " + std::to_string(e->id()) + " already has valid ideal degree: " + std::to_string(e->ideal_degree()));
					return;
				} else {
					log("Edge " + std::to_string(e->id()) + " has abnormal ideal degree: " + std::to_string(e->ideal_degree()) + ", recalculating...");
					// 重置异常值
					e->ideal_degree() = 0;
				}
			}

			if (e->boundary())
			{
				if (e->sharp())
				{
					// 检查总角度是否合理
					if (e->total_angle() < 0 || e->total_angle() > 360) {
						log("Warning: Edge " + std::to_string(e->id()) + " has abnormal total angle value: " + std::to_string(e->total_angle()) + " for ideal degree calculation");
						// 使用默认值避免计算错误
						e->ideal_degree() = 2;
					}
					else {
						int ideal_degree = round(e->total_angle() / 90.0);
						ideal_degree = ideal_degree == 0 ? 1 : ideal_degree;
						e->ideal_degree() = ideal_degree;
						log("Sharp boundary edge " + std::to_string(e->id()) + " ideal degree set to: " + std::to_string(ideal_degree) + " (total angle: " + std::to_string(e->total_angle()) + ")");
						
						if (mesh->feature_ideal_degree.find(e->sharp()) == mesh->feature_ideal_degree.end())
						{
							std::vector<E*> sharps;
							sharps.push_back(e);
							mesh->feature_ideal_degree.insert(std::pair<int, int>(e->sharp(), ideal_degree));
							log("Added feature ideal degree for sharp ID " + std::to_string(e->sharp()) + ": " + std::to_string(ideal_degree));
						}
					}
				}
				else
				{
					e->ideal_degree() = 2;
					log("Regular boundary edge " + std::to_string(e->id()) + " ideal degree set to: 2");
				}
			}
			else
			{
				e->ideal_degree() = 4;
				log("Non-boundary edge " + std::to_string(e->id()) + " ideal degree set to: 4");
			}
			
			log("Edge " + std::to_string(e->id()) + " final ideal degree: " + std::to_string(e->ideal_degree()));
		}
		catch (const std::exception& ex) {
			log("Exception in edge_ideal_degree for edge " + std::to_string(e->id()) + ": " + std::string(ex.what()));
			// 设置安全的默认值
			e->ideal_degree() = e->boundary() ? 2 : 4;
		}
		catch (...) {
			log("Unknown exception in edge_ideal_degree for edge " + std::to_string(e->id()));
			// 设置安全的默认值
			e->ideal_degree() = e->boundary() ? 2 : 4;
		}
	}

	template<typename M>
	void sheet_operation<M>::compute_edge_energy()
	{
		log("Starting edge energy computation");
		mesh->computeNormal();
		
		int totalEdges = 0;
		int boundaryEdges = 0;
		int sharpEdges = 0;
		int singularEdges = 0;
		
		log("Computing normals and angles for each edge");
		for (M::MEIterator eite(mesh); !eite.end(); eite++)
		{
			E* e = *eite;
			totalEdges++;
			
			try {
				// 记录当前处理的边ID
				if (totalEdges % 1000 == 0) {
					log("Processing edge " + std::to_string(totalEdges) + ", current edge ID: " + std::to_string(e->id()));
				}
				
				// 计算角度
				edge_angle(e);
				
				// 计算理想度数
				edge_ideal_degree(e);
				
				// 计算简化能量
				e->sim_energy() = (int)e->ideal_degree() - (int)e->neighbor_hs.size();
				
				// 统计特殊边的数量
				if (e->boundary()) boundaryEdges++;
				if (e->sharp()) sharpEdges++;
				if (e->singularity()) singularEdges++;
				
				// 记录特殊边的详细信息
				if (e->boundary() || e->sharp() || e->singularity()) {
					std::ostringstream details;
					details << "Edge " << e->id() 
						<< " (boundary: " << (e->boundary() ? "yes" : "no") 
						<< ", sharp: " << (e->sharp() ? std::to_string(e->sharp()) : "no")
						<< ", singularity: " << (e->singularity() ? "yes" : "no") << ")";
					details << " - total angle: " << e->total_angle()
						<< ", ideal degree: " << e->ideal_degree()
						<< ", actual degree: " << e->neighbor_hs.size()
						<< ", energy: " << e->sim_energy();
					log(details.str());
				}
			}
			catch (const std::exception& ex) {
				log("Exception processing edge ID " + std::to_string(e->id()) + ": " + std::string(ex.what()));
			}
			catch (...) {
				log("Unknown exception processing edge ID " + std::to_string(e->id()));
			}
		}
		
		log("Edge energy computation completed");
		log("Total edges: " + std::to_string(totalEdges) + 
			", boundary edges: " + std::to_string(boundaryEdges) + 
			", sharp edges: " + std::to_string(sharpEdges) + 
			", singular edges: " + std::to_string(singularEdges));
	}

	template<typename M>
	std::vector<std::vector<typename sheet_operation<M>::E*>> sheet_operation<M>::get_sheet_parallel_edges(std::vector<E*>sheet)
	{
		log("Starting get_sheet_parallel_edges, sheet size: " + std::to_string(sheet.size()));
		std::vector<H*> hs;
		std::vector<std::vector<E*>> parallel_es;

		// 首先重置sheet中所有边相邻的六面体的标记
		log("Initializing hexahedra marks...");
		std::set<int> processed_hex_ids; // 用于避免重复处理同一个六面体
		for (int eIndex = 0; eIndex < sheet.size(); eIndex++) {
			E* e = sheet[eIndex];
			if (!e) continue;
			
			for (int ehIndex = 0; ehIndex < e->neighbor_hs.size(); ehIndex++) {
				int hex_id = e->neighbor_hs[ehIndex];
				if (processed_hex_ids.find(hex_id) != processed_hex_ids.end()) {
					continue; // 已经处理过这个六面体
				}
				
				H* eh = mesh->idHexs(hex_id);
				if (eh) {
					eh->mark() = false; // 重置标记
					processed_hex_ids.insert(hex_id);
				}
			}
		}
		log("Finished initializing marks for " + std::to_string(processed_hex_ids.size()) + " hexahedra");
		
		// 记录初始处理信息
		int edgesProcessed = 0;
		int hexsFound = 0;
		int parallelPairsFound = 0;
		
		for (int eIndex = 0; eIndex < sheet.size(); eIndex++)
		{
			try {
				E* e = sheet[eIndex];
				if (!e) {
					log("  Warning: Null edge at index " + std::to_string(eIndex));
					continue;
				}
				
				log("  Processing edge ID: " + std::to_string(e->id()) + 
					", Adjacent hexahedra: " + std::to_string(e->neighbor_hs.size()));
				edgesProcessed++;
				
				for (int ehIndex = 0; ehIndex < e->neighbor_hs.size(); ehIndex++)
				{
					try {
						//检查六面体ID是否有效
						if (e->neighbor_hs[ehIndex] < 0 || e->neighbor_hs[ehIndex] >= mesh->maxHid()) {
							log("    Warning: Edge " + std::to_string(e->id()) + " has invalid adjacent hexahedron ID: " + std::to_string(e->neighbor_hs[ehIndex]));
							continue;
						}
						
						//获取邻接六面体
						H* eh = mesh->idHexs(e->neighbor_hs[ehIndex]);
						if (!eh) {
							log("    Warning: Unable to get adjacent hexahedron ID: " + std::to_string(e->neighbor_hs[ehIndex]));
							continue;
						}
						
						log("    Processing adjacent hexahedron ID: " + std::to_string(eh->id()));
						
						if (eh->mark()!=0) {
							log("    Hexahedron already processed, skipping");
							continue;
						}
						
						eh->mark() = true;
						hs.push_back(eh);
						hexsFound++;

						//获取平行边
						std::vector<F*> nfs;
						try {
							nfs = mesh->e_adj_f_in_hex(eh, e);
							log("    Found " + std::to_string(nfs.size()) + " adjacent faces in hexahedron");
							
							if (nfs.size() < 1) {
								log("    Warning: No adjacent faces found for edge in hexahedron");
								continue;
							}
						}
						catch(const std::exception& ex) {
							log("    Exception getting adjacent faces: " + std::string(ex.what()));
							continue;
						}
						
						F* nf = nfs[0];
						if (!nf) {
							log("    Warning: Null adjacent face");
							continue;
						}
						
						log("    Using adjacent face ID: " + std::to_string(nf->id()));

						//获取边的顶点
						if (e->vs.size() < 2) {
							log("    Warning: Edge " + std::to_string(e->id()) + " has insufficient vertices: " + std::to_string(e->vs.size()));
							continue;
						}
						
						// 检查顶点ID是否有效
						if (e->vs[0] < 0 || e->vs[0] >= mesh->maxVid() || e->vs[1] < 0 || e->vs[1] >= mesh->maxVid()) {
							log("    Warning: Edge has invalid vertex IDs: " + std::to_string(e->vs[0]) + ", " + std::to_string(e->vs[1]));
							continue;
						}
						
						V* ev1 = mesh->idVertices(e->vs[0]);
						V* ev2 = mesh->idVertices(e->vs[1]);
						
						if (!ev1 || !ev2) {
							log("    Warning: Unable to get edge vertices");
							continue;
						}
						
						log("    Edge vertices: " + std::to_string(ev1->id()) + " and " + std::to_string(ev2->id()));

						//获取第一个翻转边
						E* ne1 = nullptr;
						E* ne2 = nullptr;
						try {
							ne1 = mesh->flip_e(nf, e, ev1);
							if (!ne1) {
								log("    Warning: Failed to get first flip edge for vertex " + std::to_string(ev1->id()));
								continue;
							}
							
							ne2 = mesh->flip_e(nf, e, ev2);
							if (!ne2) {
								log("    Warning: Failed to get first flip edge for vertex " + std::to_string(ev2->id()));
								continue;
							}
							
							log("    First flip edges: " + std::to_string(ne1->id()) + " and " + std::to_string(ne2->id()));
						}
						catch(const std::exception& ex) {
							log("    Exception during edge flipping: " + std::string(ex.what()));
							continue;
						}
						
						//获取翻转面
						F* nf1 = nullptr;
						F* nf2 = nullptr;
						try {
							nf1 = mesh->flip_f(eh, nf, ne1);
							if (!nf1) {
								log("    Warning: Failed to get first flip face for edge " + std::to_string(ne1->id()));
								continue;
							}
							
							nf2 = mesh->flip_f(eh, nf, ne2);
							if (!nf2) {
								log("    Warning: Failed to get first flip face for edge " + std::to_string(ne2->id()));
								continue;
							}
							
							log("    First flip faces: " + std::to_string(nf1->id()) + " and " + std::to_string(nf2->id()));
						}
						catch(const std::exception& ex) {
							log("    Exception during face flipping: " + std::string(ex.what()));
							continue;
						}

						//找平行边对
						log("    Starting to find parallel edge pairs by flipping 4 times...");
						for (int peIndex = 0; peIndex < 4; peIndex++)
						{
							try {
								log("      Flip iteration " + std::to_string(peIndex+1));
								
								ne1 = mesh->flip_e(nf1, ne1, ev1);
								if (!ne1) {
									log("      Warning: Failed to flip edge for vertex " + std::to_string(ev1->id()) + " at iteration " + std::to_string(peIndex+1));
									break;
								}
								
								ev1 = mesh->flip_v(ne1, ev1);
								if (!ev1) {
									log("      Warning: Failed to flip vertex " + std::to_string(ev1->id()) + " at iteration " + std::to_string(peIndex+1));
									break;
								}
								
								ne2 = mesh->flip_e(nf2, ne2, ev2);
								if (!ne2) {
									log("      Warning: Failed to flip edge for vertex " + std::to_string(ev2->id()) + " at iteration " + std::to_string(peIndex+1));
									break;
								}
								
								ev2 = mesh->flip_v(ne2, ev2);
								if (!ev2) {
									log("      Warning: Failed to flip vertex " + std::to_string(ev2->id()) + " at iteration " + std::to_string(peIndex+1));
									break;
								}
								
								log("      After flipping: ne1=" + std::to_string(ne1->id()) + ", ne2=" + std::to_string(ne2->id()) +
									", ev1=" + std::to_string(ev1->id()) + ", ev2=" + std::to_string(ev2->id()));

								if (ne1->mark()) {
									log("      Edge " + std::to_string(ne1->id()) + " already processed, skipping");
									continue;
								}

								ne1->mark() = true;
								ne2->mark() = true;
								std::vector<E*> pair_es = { ne1, ne2 };
								parallel_es.push_back(pair_es);
								parallelPairsFound++;
								
								log("      Found parallel edge pair: " + std::to_string(ne1->id()) + " and " + std::to_string(ne2->id()));
							}
							catch(const std::exception& ex) {
								log("      Exception during parallel edges search iteration " + std::to_string(peIndex+1) + ": " + std::string(ex.what()));
								break;
							}
						}
					}
					catch(const std::exception& ex) {
						log("    Exception processing adjacent hex at index " + std::to_string(ehIndex) + ": " + std::string(ex.what()));
					}
					catch(...) {
						log("    Unknown exception processing adjacent hex at index " + std::to_string(ehIndex));
					}
				}
			}
			catch(const std::exception& ex) {
				log("  Exception processing edge at index " + std::to_string(eIndex) + ": " + std::string(ex.what()));
			}
			catch(...) {
				log("  Unknown exception processing edge at index " + std::to_string(eIndex));
			}
		}

		log("Processing complete. Edges processed: " + std::to_string(edgesProcessed) + 
			", hexahedra found: " + std::to_string(hexsFound) + 
			", parallel pairs found: " + std::to_string(parallelPairsFound));
		
		//重置标记
		log("Resetting marks on " + std::to_string(hs.size()) + " hexahedra");
		for (int hIndex = 0; hIndex < hs.size(); hIndex++)
		{
			try {
				H* h = hs[hIndex];
				if (!h) {
					log("  Warning: Null hexahedron at index " + std::to_string(hIndex));
					continue;
				}
				
				h->mark() = false;
				
				log("  Resetting mark on hexahedron " + std::to_string(h->id()) + " with " + std::to_string(h->es.size()) + " edges");
				for (int heIndex = 0; heIndex < h->es.size(); heIndex++)
				{
					try {
						if (h->es[heIndex] < 0 || h->es[heIndex] >= mesh->maxEid()) {
							log("    Warning: Hexahedron has invalid edge ID: " + std::to_string(h->es[heIndex]));
							continue;
						}
						
						E* he = mesh->idEdges(h->es[heIndex]);
						if (!he) {
							log("    Warning: Unable to get edge ID: " + std::to_string(h->es[heIndex]));
							continue;
						}
						
						he->mark() = false;
					}
					catch(const std::exception& ex) {
						log("    Exception resetting edge mark at index " + std::to_string(heIndex) + ": " + std::string(ex.what()));
					}
				}
			}
			catch(const std::exception& ex) {
				log("  Exception resetting hexahedron at index " + std::to_string(hIndex) + ": " + std::string(ex.what()));
			}
		}
		
		return parallel_es;
	}

	template<typename M>
	double sheet_operation<M>::predict_sheet_collapse_energy(std::vector<E*> sheet)
	{
		log("======= Starting sheet collapse energy prediction =======");
		log("Input sheet size: " + std::to_string(sheet.size()) + " edges");
		
		/*energy = ideal-real，energy越小越好*/
		double original_total_energy = 0;
		double predict_total_energy = 0;
		double cannot_collapse = -99999;//99999表示该sheet不能被折叠

		//不能折叠两个角点
		for (int eIndex = 0; eIndex < sheet.size(); eIndex++)
		{
			E* e = sheet[eIndex];
			V* ev1 = mesh->idVertices(e->vs[0]);
			V* ev2 = mesh->idVertices(e->vs[1]);
			
			// 检查feature_vertex值是否异常（内存垃圾值），如果是则重置为安全值0
			if (ev1->feature_vertex() < -1000000 || ev1->feature_vertex() > 1000000) {
				log("Warning: Vertex " + std::to_string(ev1->id()) + " has abnormal feature_vertex value: " + 
					std::to_string(ev1->feature_vertex()) + ", resetting to 0");
				ev1->feature_vertex() = 0;
			}
			
			if (ev2->feature_vertex() < -1000000 || ev2->feature_vertex() > 1000000) {
				log("Warning: Vertex " + std::to_string(ev2->id()) + " has abnormal feature_vertex value: " + 
					std::to_string(ev2->feature_vertex()) + ", resetting to 0");
				ev2->feature_vertex() = 0;
			}
			
			log("Checking vertex constraints for edge " + std::to_string(e->id()));
			log("  Vertex1 ID: " + std::to_string(ev1->id()) + 
				", Feature vertex: " + std::to_string(ev1->feature_vertex()) + 
				", Is corner: " + (ev1->corner() ? "yes" : "no"));
			log("  Vertex2 ID: " + std::to_string(ev2->id()) + 
				", Feature vertex: " + std::to_string(ev2->feature_vertex()) + 
				", Is corner: " + (ev2->corner() ? "yes" : "no"));

			if (ev1->corner() && ev2->corner())
			{
				log("Both vertices are corners, cannot collapse, returning: " + std::to_string(cannot_collapse));
				return cannot_collapse;
			}
			if (ev1->feature_vertex() < 0 && ev2->feature_vertex() < 0)
			{
				log("Both vertices have feature values < 0, cannot collapse, returning: " + std::to_string(cannot_collapse));
				return cannot_collapse;
			}
			if (ev1->feature_vertex() < 0 && ev2->feature_vertex() > 0)
			{
				int sharp_id = ev2->feature_vertex();
				int corner_id = ev1->feature_vertex();
				log("  Vertex1 feature < 0, Vertex2 feature > 0");
				log("  Checking connection between sharp_id " + std::to_string(sharp_id) + " and corner_id " + std::to_string(corner_id));
				
				try {
					std::vector<int>connect_corners_id = mesh->feature_edge_corner[sharp_id];
					log("  Connected corners count: " + std::to_string(connect_corners_id.size()));
					
					std::vector<int>::iterator ite = std::find(connect_corners_id.begin(), connect_corners_id.end(), corner_id);
					if (ite == connect_corners_id.end())
					{
						log("  Corner not connected to feature edge, cannot collapse, returning: " + std::to_string(cannot_collapse));
						return cannot_collapse;
					}
					else {
						log("  Corner connected to feature edge, can proceed with checks");
					}
				} catch (const std::exception& e) {
					log("  Exception while checking connection: " + std::string(e.what()));
					log("  Exception handling: assuming cannot collapse, returning: " + std::to_string(cannot_collapse));
					return cannot_collapse;
				}
			}
			if (ev1->feature_vertex() > 0 && ev2->feature_vertex() < 0)
			{
				int sharp_id = ev1->feature_vertex();
				int corner_id = ev2->feature_vertex();
				log("  Vertex1 feature > 0, Vertex2 feature < 0");
				log("  Checking connection between sharp_id " + std::to_string(sharp_id) + " and corner_id " + std::to_string(corner_id));
				
				try {
					std::vector<int>connect_corners_id = mesh->feature_edge_corner[sharp_id];
					log("  Connected corners count: " + std::to_string(connect_corners_id.size()));
					
					std::vector<int>::iterator ite = std::find(connect_corners_id.begin(), connect_corners_id.end(), corner_id);
					if (ite == connect_corners_id.end())
					{
						log("  Corner not connected to feature edge, cannot collapse, returning: " + std::to_string(cannot_collapse));
						return cannot_collapse;
					}
					else {
						log("  Corner connected to feature edge, can proceed with checks");
					}
				} catch (const std::exception& e) {
					log("  Exception while checking connection: " + std::string(e.what()));
					log("  Exception handling: assuming cannot collapse, returning: " + std::to_string(cannot_collapse));
					return cannot_collapse;
				}
			}

			if (ev1->feature_vertex() > 0 && ev2->feature_vertex() > 0)
			{
				log("  Both vertices have feature values > 0");
				log("  Checking if feature IDs match: " + std::to_string(ev1->feature_vertex()) + " vs " + std::to_string(ev2->feature_vertex()));
				
				if (ev1->feature_vertex() != ev2->feature_vertex())
				{
					log("  Feature IDs don't match, cannot collapse, returning: " + std::to_string(cannot_collapse));
					return cannot_collapse;
				}
				else {
					log("  Feature IDs match, can proceed with checks");
				}
			}
		}
		
		log("Vertex constraint checks passed, continuing");

		//获取sheet中的所有边
		try {
			log("Starting to get parallel edges in sheet");
			std::vector<std::vector<E*>>parallel_es = get_sheet_parallel_edges(sheet);
			log("Retrieved " + std::to_string(parallel_es.size()) + " pairs of parallel edges");

			//评估是否可以折叠
			int temp_less_num = 0;
			int temp_equal_num = 0;
			int temp_large_num = 0;
			int temp_boundary_num = 0;//估计边界折叠参数

			for (int peIndex = 0; peIndex < parallel_es.size(); peIndex++)
			{
				std::vector<E*> one_pes = parallel_es[peIndex];
				E* e1 = one_pes[0];
				E* e2 = one_pes[1];
				log("\nProcessing parallel edge pair #" + std::to_string(peIndex+1) + ":");
				log("  Edge1 ID: " + std::to_string(e1->id()) + 
					", Boundary: " + (e1->boundary() ? "yes" : "no") + 
					", Sharp: " + (e1->sharp() ? std::to_string(e1->sharp()) : "no") + 
					", Angle: " + std::to_string(e1->total_angle()) + 
					", Adjacent hexahedra: " + std::to_string(e1->neighbor_hs.size()) +
					", Energy: " + std::to_string(e1->sim_energy()));
				log("  Edge2 ID: " + std::to_string(e2->id()) + 
					", Boundary: " + (e2->boundary() ? "yes" : "no") + 
					", Sharp: " + (e2->sharp() ? std::to_string(e2->sharp()) : "no") + 
					", Angle: " + std::to_string(e2->total_angle()) + 
					", Adjacent hexahedra: " + std::to_string(e2->neighbor_hs.size()) +
					", Energy: " + std::to_string(e2->sim_energy()));
				
				if (e1->sharp() && e2->sharp())
				{
					log("  Both edges are sharp, cannot collapse, returning: " + std::to_string(cannot_collapse));
					return cannot_collapse;
				}
				
				double predict_angle = 0;
				double predict_degree = 0;

				//计算预测角度
				log("  Calculating predicted angle");
				if (e1->boundary() && e2->boundary())
				{
					if (e1->sharp() && !e2->sharp())
					{
						predict_angle = e1->total_angle();
						log("    Both edges are boundary, Edge1 is sharp, using Edge1 angle: " + std::to_string(predict_angle));
					}
					else if (e2->sharp() && !e1->sharp())
					{
						predict_angle = e2->total_angle();
						log("    Both edges are boundary, Edge2 is sharp, using Edge2 angle: " + std::to_string(predict_angle));
					}
					else
					{
						predict_angle = (e1->total_angle() + e2->total_angle()) * 0.5;
						log("    Both edges are boundary, using average angle: " + std::to_string(predict_angle));
					}
				}
				else if (e1->boundary() && !e2->boundary())
				{
					predict_angle = e1->total_angle();
					log("    Edge1 is boundary, using Edge1 angle: " + std::to_string(predict_angle));
				}
				else if (!e1->boundary() && e2->boundary())
				{
					predict_angle = e2->total_angle();
					log("    Edge2 is boundary, using Edge2 angle: " + std::to_string(predict_angle));
				}
				else
				{
					predict_angle = (e1->total_angle() + e2->total_angle()) * 0.5;
					log("    Neither edge is boundary, using average angle: " + std::to_string(predict_angle));
				}

				//计算预测度数
				log("  Calculating predicted degree");
				if (e1->boundary() && e2->boundary())
				{
					predict_degree = e1->neighbor_hs.size() + e2->neighbor_hs.size() - 2;
					log("    Both edges are boundary, degree calculation: " + std::to_string(e1->neighbor_hs.size()) + " + " + std::to_string(e2->neighbor_hs.size()) + " - 2 = " + std::to_string(predict_degree));
				}
				else
				{
					predict_degree = e1->neighbor_hs.size() + e2->neighbor_hs.size() - 4;
					log("    Degree calculation: " + std::to_string(e1->neighbor_hs.size()) + " + " + std::to_string(e2->neighbor_hs.size()) + " - 4 = " + std::to_string(predict_degree));
				}
				
				//度数过小的边禁止折叠
				log("  Checking if degree is too small");
				if (e1->boundary() || e2->boundary())
				{
					if (predict_degree < 1)
					{
						log("    Contains boundary edge and predicted degree < 1, cannot collapse, returning: " + std::to_string(cannot_collapse));
						return cannot_collapse;
					}
				}
				else
				{
					if (predict_degree < 2)
					{
						log("    Predicted degree < 2, cannot collapse, returning: " + std::to_string(cannot_collapse));
						return cannot_collapse;
					}
				}

				//计算能量
				int temp_predict_energy = abs(round(predict_angle / 90.0) - predict_degree);
				log("  Calculating predicted energy: abs(round(" + std::to_string(predict_angle) + " / 90.0) - " + std::to_string(predict_degree) + ") = " + std::to_string(temp_predict_energy));
				predict_total_energy += temp_predict_energy;

				int temp_original_energy1 = abs(e1->sim_energy());
				int temp_original_energy2 = abs(e2->sim_energy());
				original_total_energy += temp_original_energy1;
				original_total_energy += temp_original_energy2;
				
				log("  Cumulative original energy: " + std::to_string(original_total_energy) + ", cumulative predicted energy: " + std::to_string(predict_total_energy));
				log("  Edge1 energy: " + std::to_string(temp_original_energy1) + ", Edge2 energy: " + std::to_string(temp_original_energy2));

				//边界能量		
				int boundary_energy1 = 0;
				int boundary_energy2 = 0;
				if (e1->boundary())
				{
					boundary_energy1 = temp_original_energy1 - temp_predict_energy;
					log("  Edge1 boundary energy calculation: " + std::to_string(temp_original_energy1) + " - " + std::to_string(temp_predict_energy) + " = " + std::to_string(boundary_energy1));
				}
				if (e2->boundary())
				{
					boundary_energy2 = temp_original_energy2 - temp_predict_energy;
					log("  Edge2 boundary energy calculation: " + std::to_string(temp_original_energy2) + " - " + std::to_string(temp_predict_energy) + " = " + std::to_string(boundary_energy2));
				}
				
				//估计边界折叠能量
				if ((e1->boundary() || e2->boundary()) && !(e1->boundary() && e2->boundary()))
				{
					log("  One edge is boundary while the other is not, estimating boundary collapse energy");
					if (temp_boundary_num <= 0)
					{
						if (boundary_energy1 > 0 || boundary_energy2 > 0)
						{
							temp_boundary_num = 1;
							log("    Boundary energy > 0, setting temp_boundary_num = 1");
						}
						else if (boundary_energy1 == 0 && boundary_energy2 == 0)
						{
							log("    Boundary energy = 0, temp_boundary_num remains: " + std::to_string(temp_boundary_num));
						}
						else
						{
							temp_boundary_num = -1;
							log("    Boundary energy < 0, setting temp_boundary_num = -1");
						}
					}
				}

				//估计折叠能量
				log("  Estimating collapse energy");
				if (temp_original_energy1 > temp_predict_energy || temp_original_energy2 > temp_predict_energy)
				{
					temp_large_num++;
					log("    Original energy > predicted energy, temp_large_num increased to: " + std::to_string(temp_large_num));
				}
				else if (temp_original_energy1 == temp_predict_energy || temp_original_energy2 == temp_predict_energy)
				{
					temp_equal_num++;
					log("    Original energy = predicted energy, temp_equal_num increased to: " + std::to_string(temp_equal_num));
				}
				else if (temp_original_energy1 < temp_predict_energy && temp_original_energy2 < temp_predict_energy)
				{
					temp_less_num++;
					log("    Original energy < predicted energy, temp_less_num increased to: " + std::to_string(temp_less_num));
				}
			}
			
			log("\nCollapse energy calculation summary:");
			log("  temp_large_num: " + std::to_string(temp_large_num) + 
				", temp_equal_num: " + std::to_string(temp_equal_num) + 
				", temp_less_num: " + std::to_string(temp_less_num) + 
				", temp_boundary_num: " + std::to_string(temp_boundary_num));
			log("  Total original energy: " + std::to_string(original_total_energy) + 
				", total predicted energy: " + std::to_string(predict_total_energy) + 
				", energy difference: " + std::to_string(original_total_energy - predict_total_energy));
			
			if (temp_boundary_num < 0)
			{
				log("Boundary energy check: temp_boundary_num < 0, not suitable for collapse, returning: -999");
				return -999;
			}
			else if (temp_large_num > 0)
			{
				log("Energy comparison check: temp_large_num > 0, suitable for collapse, returning final energy difference: " + std::to_string(original_total_energy - predict_total_energy));
				return original_total_energy - predict_total_energy;
			}
			else if (temp_less_num == 0)
			{
				log("Energy comparison check: temp_less_num = 0, suitable for collapse, returning final energy difference: " + std::to_string(original_total_energy - predict_total_energy));
				return original_total_energy - predict_total_energy;
			}
			else
			{
				log("Other cases, not suitable for collapse, returning: -999");
				return -999;
			}
		}
		catch (const std::exception& e) {
			log("Exception occurred during collapse energy prediction: " + std::string(e.what()));
			return cannot_collapse;
		}
		catch (...) {
			log("Unknown exception occurred during collapse energy prediction");
			return cannot_collapse;
		}
	}



}