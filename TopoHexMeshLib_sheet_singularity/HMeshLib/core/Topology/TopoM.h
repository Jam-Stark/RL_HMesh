#ifndef TOPOLOGY_MESH_H
#define TOPOLOGY_MESH_H
#include<vector>
#include<map>
#include<unordered_map>
#include<list>
#include "../Parser/StrUtil.h"

#include <cmath> // for acos, abs, round
#include <limits> // for DBL_MIN if needed later, though not directly for angle

#ifndef M_PI // Define PI if not available
#define M_PI 3.14159265358979323846
#endif

namespace HMeshLib
{
	/*class topoV;
	class topoE;
	class topoF;
	class topoH;*/

	template<typename V,typename E, typename F,typename H>
	class topoM
	{
	public:
		topoM() 
		{
			m_maxVertexId = 0;
			m_maxEdgeId = 0;
			m_maxFaceId = 0;
			m_maxHexId = 0;
		};
		~topoM() 
		{

		};

		/*value*/
		std::list<V*> vs;
		std::list<E*> es;
		std::list<F*> fs;
		std::list<H*> hs;
		
		std::unordered_map<int, V*> m_map_vertices;
		std::unordered_map<int, E*> m_map_edges;
		std::unordered_map<int, F*> m_map_faces;
		std::unordered_map<int, H*> m_map_hexs;

		/*singularities*/
		std::vector<E*> singularities;

		/*funcation*/
		V* idVertices(int vid) { 
			//return m_map_vertices[vid];
			std::unordered_map<int, V*>::iterator iter = m_map_vertices.find(vid);
			if (iter!= m_map_vertices.end())
			{
				return iter->second;
			}
			else 
				return NULL;
		};
		E* idEdges(int eid) { 
			//return m_map_edges[eid];
			std::unordered_map<int, E*>::iterator iter = m_map_edges.find(eid);
			if (iter != m_map_edges.end())
			{
				return iter->second;
			}
			else
				return NULL;
		};
		F* idFaces(int fid) { 
			//return m_map_faces[fid];
			std::unordered_map<int, F*>::iterator iter = m_map_faces.find(fid);
			if (iter != m_map_faces.end())
			{
				return iter->second;
			}
			else
				return NULL;
		};
		H* idHexs(int hid) { 
			//return m_map_hexs[hid];
			std::unordered_map<int, H*>::iterator iter = m_map_hexs.find(hid);
			if (iter != m_map_hexs.end())
			{
				return iter->second;
			}
			else
				return NULL; 
		};
		E* VerticesEdge(V* v1, V* v2);
		F* VerticesFace(std::vector<int> vid);

		/*normalize*/
		void normalize();
		void construct_hex(H* h, int hid, int* v);
		void construct_hex(H* h, int hid, std::vector<int> v);
		void load_Qhex(const char* input);
		void faceNormal();
		void vertexNormal();
		void computeNormal();
		void write_Qhex(const char* output);
		

		/*flip funcation*/
		V* flip_v(E* edge, V* v);
		E* flip_e(F* face, E* e, V* v);
		F* flip_f(H* hex, F* face, E* edge);
		H* flip_H(H* hex, F* face);
		/*Element neighbor relationship*/
		V* v_adj_v_in_face(F* face, E* e);
		E* e_parallel_e_in_face(F* face, E* e);
		std::vector<E*> e_parallel_e_in_hex(H* hex, E* e);
		std::vector<F*> e_adj_f_in_hex(H* hex, E* e);
		std::vector<E*> v_adj_e_in_hex(H* hex, V* v);
		std::vector<F*> v_adj_f_in_hex(H* hex, V* v);
		F* f_parallel_f_in_hex(H* hex ,F* face);
		std::vector<E*> f_adj_e_in_hex(H* hex, F* face);

		/*mark attributes*/
		std::vector<E*> mark_singularity();
		void push_singularity(E* e);
		void delete_singularity(E* e);
		void mark_singularity_node();
		void mark_boundary();
		int& maxVid() { return m_maxVertexId; };
		int& maxEid() { return m_maxEdgeId; };
		int& maxFid() { return m_maxFaceId; };
		int& maxHid() { return m_maxHexId; };
		
		/*delete element*/
		void delete_vertex(V* v);
		void delete_edge(E* e);
		void delete_face(F* f);
		void delete_hex(H* h);

		/*change a element in Hex*/
		bool change_Hex_V(H* h, V* source, V* target);
		bool change_Hex_E(H* h, E* source, E* target);
		bool change_Hex_F(H* h, F* source, F* target);

		/*change a element in face*/
		bool change_Face_V(F* f, V* source, V* target);
		bool change_Face_E(F* f, E* source, E* target);
		
		/*change a element in edge*/
		bool change_Edge_V(E* e, V* source, V* target);

		/*revise normal*/
		void revise_hex_face_normal(H* h);
		void revise_all_face_normal();
		/*output mesh*/
		void output_Quad_mesh(const char* fileName);

		/*clear*/
		void _clear();

		// 打印map信息的辅助函数
		void print_vertices_map();
		void print_edges_map();
		void print_faces_map();
		void print_hexs_map();

		// 打印特定ID的顶点、边、面、六面体信息
		void print_vertex(int id);
		void print_edge(int id);
		void print_face(int id);
		void print_hex(int id);

		double vector_angle(CPoint a, CPoint b);
		void edge_angle(E* e);
		void compute_sharp_edges(double threshold_degrees);
		void compute_features_and_corners(double sharp_threshold_degrees);

		void compute_features()
		{
					//mark boundary
		mark_boundary();
		singularities = mark_singularity();
		/*compute the face normal*/
		computeNormal();

		for (std::list<E*>::iterator eite = this->es.begin(); eite != this->es.end(); eite++) {

			if (*eite)
            	this->edge_angle(*eite);
            
		}
		compute_sharp_edges(30.0);

		compute_features_and_corners(30.0);
		}

		int count_sharp_edges()
		{
			int sharp_count = 0;
			// 遍历网格中的所有边
			for (typename std::list<E*>::iterator eite = this->es.begin(); eite != this->es.end(); eite++)
			{
				E* e = *eite;
				// 检查边指针是否有效，以及 sharp() 属性是否大于 0 (假设 > 0 表示 sharp)
				if (e && e->sharp() > 0)
				{
					sharp_count++;
				}
			}
			return sharp_count;
		}
		int count_corner()
		{
			int corner_count = 0;
			for(typename std::list<V*>::iterator eite = this->vs.begin(); eite != this->vs.end(); ++eite) {
				V* v = *eite;

				if(v&& v->corner()==true)
					corner_count++;
			}
			return corner_count;
		}
		int count_boundary_byF()
		{
			int boundary_count = 0;
			for(typename std::list<V*>::iterator eite = this->vs.begin(); eite != this->vs.end(); ++eite) {
				V* v = *eite;

				if(v&& v->boundary()==true)
					boundary_count++;
			}
			return boundary_count;
		}
		int count_boundary_byE()
		{
			int boundary_count = 0;
			for(typename std::list<E*>::iterator eite = this->es.begin(); eite != this->es.end(); ++eite) {
				E* e = *eite;

				if(e&& e->boundary()==true)
					boundary_count++;
			}
			return boundary_count;
		}

		topoM& operator=(const topoM& other) {
			if (this == &other) {
				return *this; // 处理自赋值
			}
	
			// 1. 清理当前对象的状态 (释放旧内存)
			this->_clear(); // 你已有的清理函数
	
			// 2. 复制简单成员变量 (同拷贝构造函数)
			m_nVertices = other.m_nVertices;
			m_nHexs = other.m_nHexs;
			m_nSharpEdges = other.m_nSharpEdges;
			m_nCorners = other.m_nCorners;
			m_nFaces = other.m_nFaces;
			m_maxVertexId = other.m_maxVertexId;
			m_maxEdgeId = other.m_maxEdgeId;
			m_maxFaceId = other.m_maxFaceId;
			m_maxHexId = other.m_maxHexId;
	
			// 3. 创建拓扑元素的深拷贝并建立映射 (同拷贝构造函数)
			std::unordered_map<int, V*> old_v_to_new_v;
			std::unordered_map<int, E*> old_e_to_new_e;
			std::unordered_map<int, F*> old_f_to_new_f;
			std::unordered_map<int, H*> old_h_to_new_h;
	
			// 3a. 拷贝顶点
			for (const auto& pair : other.m_map_vertices) {
				const V* old_v = pair.second;
				if (!old_v) continue;
				V* new_v = new V(*old_v);
				new_v->neighbor_es.clear(); 
				new_v->neighbor_fs.clear();
				new_v->neighbor_hs.clear();
				this->vs.push_back(new_v);
				this->m_map_vertices[new_v->id()] = new_v;
				old_v_to_new_v[old_v->id()] = new_v;
			}
			// 3b. 拷贝边
			for (const auto& pair : other.m_map_edges) {
				const E* old_e = pair.second;
				if (!old_e) continue;
				E* new_e = new E(*old_e);
				new_e->neighbor_fs.clear();
				new_e->neighbor_hs.clear();
				this->es.push_back(new_e);
				this->m_map_edges[new_e->id()] = new_e;
				old_e_to_new_e[old_e->id()] = new_e;
			}
			// 3c. 拷贝面
			for (const auto& pair : other.m_map_faces) {
				const F* old_f = pair.second;
				if (!old_f) continue;
				F* new_f = new F(*old_f);
				new_f->neighbor_hs.clear();
				this->fs.push_back(new_f);
				this->m_map_faces[new_f->id()] = new_f;
				old_f_to_new_f[old_f->id()] = new_f;
			}
			// 3d. 拷贝六面体
			for (const auto& pair : other.m_map_hexs) {
				const H* old_h = pair.second;
				if (!old_h) continue;
				H* new_h = new H(*old_h);
				this->hs.push_back(new_h);
				this->m_map_hexs[new_h->id()] = new_h;
				old_h_to_new_h[old_h->id()] = new_h;
			}
	
			// 4. 重建邻接关系 (同拷贝构造函数中的步骤 3)
			// 4a. 更新顶点的邻接列表
			for (V* new_v : this->vs) {
				const V* old_v = nullptr;
				auto it_old_v = other.m_map_vertices.find(new_v->id());
				if (it_old_v != other.m_map_vertices.end()) {
					old_v = it_old_v->second;
					 if (old_v) {
						for (int old_eid : old_v->neighbor_es) {
							if (old_e_to_new_e.count(old_eid)) {
							   new_v->neighbor_es.push_back(old_eid);
							}
						}
						for (int old_fid : old_v->neighbor_fs) {
							if (old_f_to_new_f.count(old_fid)) {
							   new_v->neighbor_fs.push_back(old_fid);
							}
						}
						for (int old_hid : old_v->neighbor_hs) {
							if (old_h_to_new_h.count(old_hid)) {
							   new_v->neighbor_hs.push_back(old_hid);
							}
						}
					}
				}
			}
			// 4b. 更新边的邻接列表
			 for (E* new_e : this->es) {
				const E* old_e = nullptr;
				auto it_old_e = other.m_map_edges.find(new_e->id());
				if (it_old_e != other.m_map_edges.end()) {
					old_e = it_old_e->second;
					 if (old_e) {
						for (int old_fid : old_e->neighbor_fs) {
							 if (old_f_to_new_f.count(old_fid)) {
								new_e->neighbor_fs.push_back(old_fid);
							 }
						}
						for (int old_hid : old_e->neighbor_hs) {
							if (old_h_to_new_h.count(old_hid)) {
								new_e->neighbor_hs.push_back(old_hid);
							}
						}
					}
				}
			}
			// 4c. 更新面的邻接列表
			for (F* new_f : this->fs) {
				const F* old_f = nullptr;
				auto it_old_f = other.m_map_faces.find(new_f->id());
				if (it_old_f != other.m_map_faces.end()) {
					old_f = it_old_f->second;
					if (old_f) {
						for (int old_hid : old_f->neighbor_hs) {
							if (old_h_to_new_h.count(old_hid)) {
								new_f->neighbor_hs.push_back(old_hid);
							}
						}
					}
				}
			}
	
			// 5. 复制奇异线列表 (同拷贝构造函数)
			singularities.clear();
			for (E* old_sing_e : other.singularities) {
				 if (old_sing_e && old_e_to_new_e.count(old_sing_e->id())) {
					this->singularities.push_back(old_e_to_new_e[old_sing_e->id()]);
				}
			}
	
			return *this;
		}

	protected:
		/*! number of vertices */
		int m_nVertices;
		/*! number of Hexs */
		int m_nHexs;
		/*! num of sharp edges*/
		int m_nSharpEdges;
		/*! num of corner*/
		int m_nCorners;
		/*! num of faces */
		int m_nFaces;

		/*! max vertex id */
		int m_maxVertexId;
		/* max Edge id*/
		int m_maxEdgeId;
		/*max Face id */
		int m_maxFaceId;
		/*max hex id*/
		int m_maxHexId;		
	};

	template<typename V, typename E, typename F, typename H>
	void topoM<V, E, F, H>::compute_features_and_corners(double sharp_threshold_degrees = 30.0)
	{
		//topoMeshLog("--- Starting Feature and Corner Computation ---"); // 使用您选择的日志函数
	
		// --- 1. Ensure necessary pre-computations are done ---
		// Mark boundary first, as it affects sharpness criteria
		// //topoMeshLog("Step 1: Marking boundaries...");
		// this->mark_boundary(); // 确保已调用
	
		// // Compute face normals, needed for dihedral angles
		// topoMeshLog("Step 2: Computing normals...");
		// this->computeNormal(); // 确保已调用
	
		// Compute angles around edges - 必须在 compute_sharp_edges 之前调用!
		// topoMeshLog("Step 3: Computing edge angles...");
		// for (auto it = this->es.begin(); it != this->es.end(); ++it) {
		// 	 E* e = *it;
		// 	 if (e) {
		// 		 // 假设 edge_angle 函数存在于 topoM 或可以被 topoM 访问
		// 		 // 您可能需要从 sheet_operation 移动或复制 edge_angle 函数
		// 		 this->edge_angle(e); // <--- 确保这个函数被调用
		// 	 }
		// }
		//  topoMeshLog("Edge angles computed.");
	
	
		// // --- 2. Compute Sharp Edges ---
		// topoMeshLog("Step 4: Computing sharp edges (threshold: " + std::to_string(sharp_threshold_degrees) + " degrees)...");
		// // 调用或实现 compute_sharp_edges 函数
		// this->compute_sharp_edges(sharp_threshold_degrees); // 确保此函数正确设置 e->sharp()
		// int sharp_count = this->count_sharp_edges(); // 假设有此函数
		// topoMeshLog("Sharp edges computed. Found " + std::to_string(sharp_count) + " sharp edges.");
	
	
		// --- 3. Compute Corner Vertices ---
		//topoMeshLog("Step 5: Computing corner vertices (threshold >= 3 sharp edges)...");
		int corner_count = 0;
		for (auto it = this->vs.begin(); it != this->vs.end(); ++it)
		{
			V* v = *it;
			if (!v) continue;
	
			int incident_sharp_edges = 0;
			for (int neighbor_eid : v->neighbor_es)
			{
				E* neighbor_e = this->idEdges(neighbor_eid);
				// 检查边是否存在及其 sharp 属性
				if (neighbor_e && neighbor_e->sharp() > 0)
				{
					incident_sharp_edges++;
				}
			}
	
			// 定义角点条件（例如，>= 3 条锐利边）
			if (incident_sharp_edges >= 3)
			{
				v->corner() = true; // 设置 corner 标志
				corner_count++;
				 // topoMeshLog("  Vertex " + std::to_string(v->id()) + " marked as corner (" + std::to_string(incident_sharp_edges) + " sharp edges)."); // Optional detailed topoMeshLog
			}
			else
			{
				v->corner() = false; // 确保非角点被明确标记为 false
			}
		}
		//topoMeshLog("Corner vertices computed. Found " + std::to_string(corner_count) + " corner vertices.");
		//topoMeshLog("--- Feature and Corner Computation Finished ---");
	}

	template<typename V, typename E, typename F, typename H>
	void topoM< V, E, F, H>::compute_sharp_edges(double threshold_degrees)
    {
        //topoMeshLog("compute_sharp_edges: Starting automatic sharp edge computation with threshold: " + std::to_string(threshold_degrees));
        int sharp_count = 0;
        int edge_processed_count = 0;

        for (std::list<E*>::iterator eite = this->es.begin(); eite != this->es.end(); eite++)
        {
            edge_processed_count++;
            E* e = *eite;
            if (!e) {
                //topoMeshLog("  compute_sharp_edges: Skipping null edge pointer at index " + std::to_string(edge_processed_count -1));
                continue;
            }
            // //topoMeshLog("  compute_sharp_edges: Processing Edge ID " + std::to_string(e->id()));

            // Ensure edge_angle has been computed (critical pre-requisite)
            // Check if total_angle seems reasonable (e.g., not default 0 if boundary)
            if (e->boundary() && e->total_angle() == 0.0) {
                 //topoMeshLog("    compute_sharp_edges: Warning - Edge " + std::to_string(e->id()) + " is boundary but total_angle is 0. Angle might not have been computed correctly. Recalculating...");
                 edge_angle(e); // Attempt to calculate angle now
                 if(e->total_angle() == 0.0) { // Check again
                     //topoMeshLog("    compute_sharp_edges: Error - Failed to calculate angle for boundary edge " + std::to_string(e->id()) + ". Cannot determine sharpness.");
                     e->sharp() = 0; // Default to not sharp if angle calculation failed
                     continue;
                 }
            }


            if (e->boundary())
            {
                 // //topoMeshLog("    compute_sharp_edges: Edge " + std::to_string(e->id()) + " is boundary. total_angle = " + std::to_string(e->total_angle()));
                double deviation = abs(e->total_angle() - 180.0);
                 // Clamp deviation for robustness, although angles should ideally be within [0, 360]
                 deviation = std::min(deviation, 180.0);
                 // //topoMeshLog("    compute_sharp_edges: Deviation from 180 = " + std::to_string(deviation));

                if (deviation > threshold_degrees)
                {
                    e->sharp() = 1;
                    sharp_count++;
                    //topoMeshLog("    compute_sharp_edges: Marked edge " + std::to_string(e->id()) + " as SHARP (deviation=" + std::to_string(deviation) + ")");
                }
                else
                {
                    e->sharp() = 0;
                    // //topoMeshLog("    compute_sharp_edges: Marked edge " + std::to_string(e->id()) + " as NOT sharp.");
                }
            }
            else // Internal Edge
            {
                // //topoMeshLog("    compute_sharp_edges: Edge " + std::to_string(e->id()) + " is internal.");
                // Current logic: Internal edges are not marked sharp based on dihedral angle alone
                 e->sharp() = 0;
                 // //topoMeshLog("    compute_sharp_edges: Marked edge " + std::to_string(e->id()) + " as NOT sharp (internal).");

                // --- Optional: Uncomment below to mark internal edges based on 360 deviation ---
                /*
                double deviation = abs(e->total_angle() - 360.0);
                if (deviation > 180.0) deviation = 360.0 - deviation; // Handle wrap-around
                deviation = std::min(deviation, 180.0);
                //topoMeshLog("      compute_sharp_edges: Internal edge deviation from 360 = " + std::to_string(deviation));
                if (deviation > threshold_degrees) {
                    e->sharp() = 1; // Or a different value for internal sharp?
                    sharp_count++;
                    //topoMeshLog("      compute_sharp_edges: Marked internal edge " + std::to_string(e->id()) + " as SHARP (deviation=" + std::to_string(deviation) + ")");
                } else {
                    e->sharp() = 0;
                }
                */
                // --- End Optional ---
            }
        } // End edge loop
        //topoMeshLog("compute_sharp_edges: Sharp edge computation finished. Processed " + std::to_string(edge_processed_count) + " edges. Marked " + std::to_string(sharp_count) + " edges as sharp.");
    }
	template<typename V, typename E, typename F, typename H>
	double topoM< V, E, F, H>::vector_angle(CPoint a, CPoint b)
    {
        // //topoMeshLog("  vector_angle: Input a(" + std::to_string(a[0]) + "," + std::to_string(a[1]) + "," + std::to_string(a[2]) + ")");
        // //topoMeshLog("  vector_angle: Input b(" + std::to_string(b[0]) + "," + std::to_string(b[1]) + "," + std::to_string(b[2]) + ")");

        double norm_a = a.norm();
        double norm_b = b.norm();
        // //topoMeshLog("  vector_angle: norm_a = " + std::to_string(norm_a) + ", norm_b = " + std::to_string(norm_b));

        if (norm_a < 1e-10 || norm_b < 1e-10) {
             //topoMeshLog("  vector_angle: Warning - Zero vector detected. Returning 0.0");
             return 0.0;
        }

        CPoint norm_vec_a = a / norm_a;
        CPoint norm_vec_b = b / norm_b;
        double temp = norm_vec_a * norm_vec_b;
        // //topoMeshLog("  vector_angle: Dot product (before clamp) = " + std::to_string(temp));

        temp = std::max(-1.0, std::min(1.0, temp));
        // //topoMeshLog("  vector_angle: Dot product (after clamp) = " + std::to_string(temp));

        double angle_rad = acos(temp);
        // //topoMeshLog("  vector_angle: Angle (radians) = " + std::to_string(angle_rad));

        if (!std::isfinite(angle_rad)) {
            //topoMeshLog("  vector_angle: Warning - Invalid angle from acos. Returning 0.0");
            return 0.0;
        }

        double angle_deg = angle_rad / M_PI * 180.0;
        // //topoMeshLog("  vector_angle: Angle (degrees) = " + std::to_string(angle_deg));
        return angle_deg;
    }

	template<typename V, typename E, typename F, typename H>
    void topoM< V, E, F, H>::edge_angle(E* e)
    {
        if (!e) {
            //topoMeshLog("edge_angle: Error - Input edge pointer is null.");
            return;
        }
        //topoMeshLog("edge_angle: Calculating for Edge ID " + std::to_string(e->id()) + ", Boundary: " + (e->boundary() ? "Yes" : "No"));

        // Reset angle calculation for this edge
        e->total_angle() = 0.0;
        e->ave_angle() = 0.0;
        // //topoMeshLog("  edge_angle: Angles reset for Edge " + std::to_string(e->id()));

        if (!e->boundary())
        {
            // Internal edge, ideally 360 degrees total angle
            e->total_angle() = 360.0;
            // Calculate average based on actual neighbors, handle division by zero
            size_t neighbor_count = e->neighbor_hs.size();
            e->ave_angle() = (neighbor_count > 0) ? (360.0 / neighbor_count) : 90.0; // Default to 90 if no hex neighbors?
            //topoMeshLog("  edge_angle: Internal Edge " + std::to_string(e->id()) + ". Set total_angle=360.0, ave_angle=" + std::to_string(e->ave_angle()));
            return;
        }

        // Boundary edge calculation
        double total_angle_sum = 0.0;
        int valid_hex_count = 0;
        //topoMeshLog("  edge_angle: Boundary Edge " + std::to_string(e->id()) + ". Processing " + std::to_string(e->neighbor_hs.size()) + " neighboring hexes.");

        for (int ehIndex = 0; ehIndex < e->neighbor_hs.size(); ehIndex++)
        {
            int hex_id = e->neighbor_hs[ehIndex];
            H* eh = idHexs(hex_id);
            if (!eh) {
                //topoMeshLog("    edge_angle: Warning - Could not find Hex ID " + std::to_string(hex_id) + " for Edge " + std::to_string(e->id()));
                continue;
            }
            //topoMeshLog("    edge_angle: Processing Hex ID " + std::to_string(eh->id()));

            // Get adjacent faces within this hex
            std::vector<F*> ehfs = e_adj_f_in_hex(eh, e);
            if (ehfs.size() != 2) {
                //topoMeshLog("    edge_angle: Warning - Expected 2 adjacent faces in Hex " + std::to_string(eh->id()) + " for Edge " + std::to_string(e->id()) + ", found " + std::to_string(ehfs.size()) + ". Skipping hex.");
                continue;
            }

            F* f1 = ehfs[0];
            F* f2 = ehfs[1];
            if (!f1 || !f2) {
                 //topoMeshLog("    edge_angle: Warning - Null face pointer found in Hex " + std::to_string(eh->id()) + " for Edge " + std::to_string(e->id()) + ". Skipping hex.");
                 continue;
            }
            //topoMeshLog("    edge_angle: Adjacent faces: F" + std::to_string(f1->id()) + " and F" + std::to_string(f2->id()));

            CPoint n1_raw = f1->normal();
            CPoint n2_raw = f2->normal();
            double n1_norm = n1_raw.norm();
            double n2_norm = n2_raw.norm();
            //topoMeshLog("      edge_angle: Raw normals - n1_norm=" + std::to_string(n1_norm) + ", n2_norm=" + std::to_string(n2_norm));

            // Ensure normals are valid before proceeding
            if (n1_norm < 1e-10 || n2_norm < 1e-10) {
                //topoMeshLog("      edge_angle: Warning - Invalid normal detected (norm < 1e-10). Skipping angle calculation for this hex.");
                continue;
            }

            CPoint n1 = n1_raw;
            CPoint n2 = n2_raw;

            // Ensure correct normal orientation relative to the hex 'eh'
            bool n1_flipped = false;
            bool n2_flipped = false;
            if (f1->neighbor_hs.empty()){
                //topoMeshLog("      edge_angle: Warning - Face " + std::to_string(f1->id()) + " has no hex neighbors.");
                 continue; // Should not happen if it's adjacent to eh, indicates data inconsistency
            } else if (f1->neighbor_hs[0] != eh->id()) {
                n1 = -n1;
                n1_flipped = true;
            }
            if (f2->neighbor_hs.empty()){
                //topoMeshLog("      edge_angle: Warning - Face " + std::to_string(f2->id()) + " has no hex neighbors.");
                 continue; // Should not happen if it's adjacent to eh
            } else if (f2->neighbor_hs[0] != eh->id()) {
                n2 = -n2;
                n2_flipped = true;
            }
            //topoMeshLog("      edge_angle: Orientation check - n1_flipped=" + std::string(n1_flipped ? "Yes" : "No") + ", n2_flipped=" + std::string(n2_flipped ? "Yes" : "No"));


            double angle_between_normals = vector_angle(n1, n2); // vector_angle handles normalization now
            //topoMeshLog("      edge_angle: Angle between oriented normals = " + std::to_string(angle_between_normals) + " degrees.");

            // The angle *inside* the material is 180 - angle_between_normals
            double internal_angle = 180.0 - angle_between_normals;

            // Clamp internal angle to a reasonable range (e.g., 0 to 360, though typically < 360 for single hex)
            internal_angle = std::max(0.0, std::min(360.0, internal_angle));
            //topoMeshLog("      edge_angle: Calculated internal angle = " + std::to_string(internal_angle) + " degrees.");

            // Check if calculated angle is finite
            if (!std::isfinite(internal_angle)) {
                //topoMeshLog("      edge_angle: Warning - Calculated internal angle is not finite. Skipping this hex.");
                continue;
            }

            // Store angle for this hex contribution in the Hex object
            // Check if edgeIndex is valid before using it
            if (eh->edgeIndex(e->id()) != -1) {
                 eh->edge_angle(e->id()) = internal_angle;
            } else {
                 //topoMeshLog("      edge_angle: Warning - Edge ID " + std::to_string(e->id()) + " not found in Hex ID " + std::to_string(eh->id()) + " edge list. Cannot store hex-specific angle.");
            }

            total_angle_sum += internal_angle;
            valid_hex_count++;
            //topoMeshLog("      edge_angle: Added angle to sum. Current sum = " + std::to_string(total_angle_sum) + ", valid hex count = " + std::to_string(valid_hex_count));
        } // end loop over neighboring hexes

        if (valid_hex_count > 0) {
            // Final validation on total angle sum
            if (!std::isfinite(total_angle_sum) || total_angle_sum < 0) {
                 //topoMeshLog("  edge_angle: Warning - Final total_angle_sum (" + std::to_string(total_angle_sum) + ") is invalid for Edge " + std::to_string(e->id()) + ". Resetting to 180.0.");
                 total_angle_sum = 180.0; // Default boundary angle
            }
             e->total_angle() = total_angle_sum;
             e->ave_angle() = total_angle_sum / valid_hex_count;
        } else {
             //topoMeshLog("  edge_angle: Warning - No valid hex neighbors contributed angles for Boundary Edge " + std::to_string(e->id()) + ". Setting default angles.");
             e->total_angle() = 180.0; // Default boundary angle
             e->ave_angle() = 90.0; // Default average
        }
        //topoMeshLog("  edge_angle: Final calculated angles for Edge " + std::to_string(e->id()) + ": total_angle=" + std::to_string(e->total_angle()) + ", ave_angle=" + std::to_string(e->ave_angle()));
    }

	template<typename V, typename E, typename F, typename H>
	E* topoM< V, E, F, H>::VerticesEdge(V* v1, V* v2)
	{
		E* resultE = NULL;
		for (int veIndex = 0; veIndex < v1->neighbor_es.size(); veIndex++)
		{
			E* e = idEdges(v1->neighbor_es[veIndex]);
			if (e->vertexIndex(v2->id())!=-1)
			{
				resultE = e;
				break;
			}
		}
		return resultE;
	}

	template<typename V, typename E, typename F, typename H>
	F* topoM< V, E, F, H>::VerticesFace(std::vector<int> vid)
	{
		F* result = NULL;
		V* v = idVertices(vid[0]);
		if (v==NULL)
		{
			return result;
		}
		for (int vfIndex = 0; vfIndex < v->neighbor_fs.size(); vfIndex++)
		{
			F* f = idFaces(v->neighbor_fs[vfIndex]);
			if (f->is_v_equal(vid))
			{
				result = f;
				break;
			}
		}
		return result;
	}

	template<typename V, typename E, typename F, typename H>
	V* topoM< V, E, F, H>::flip_v(E* edge, V* v)
	{
		int vIndex = edge->vertexIndex(v->id());
		if (vIndex ==-1)
		{
			return NULL;
		}
		if (vIndex==0)
		{
			return idVertices(edge->vs[1]);
		}
		else
		{
			return idVertices(edge->vs[0]);
		}
	}
	
	template<typename V, typename E, typename F, typename H>
	E* topoM< V, E, F, H>:: flip_e(F* face, E* e, V* v)
	{
		int eIndex = face->edgeIndex(e->id());
		if (eIndex==-1)
		{
			return NULL;
		}
		for (int i = 0; i < 2; i++)
		{
			int adjEId = face->es[face->eadje[eIndex][i]];
			
			E* adjE = idEdges(adjEId);
			if (adjE->vertexIndex(v->id())!=-1)
			{
				return adjE;
			}
		}
		return NULL;
	}

	template<typename V, typename E, typename F, typename H>
	F* topoM< V, E, F, H>::flip_f(H* hex, F* face, E* edge)
	{
		
		int fIndex = hex->faceIndex(face->id());
		if (fIndex==-1)
		{
			return NULL;
		}
		for (int i = 0; i < 4; i++)
		{
			int adjFid = hex->fs[hex->fadjf[fIndex][i]];
			F* adjF = idFaces(adjFid);
			if (adjF->edgeIndex(edge->id())!=-1)
			{		
				return adjF;
			}
		}
		return NULL;
	}	
	
	template<typename V, typename E, typename F, typename H>
	H* topoM< V, E, F, H>::flip_H(H* hex, F* face)
	{
		int hIndex = face->neighborHIndex(hex->id());
		if (hIndex ==-1 || face->neighbor_hs.size() == 1)
		{		
			return NULL;
		}
		else
		{
			H* result = hIndex == 0 ? idHexs(face->neighbor_hs[1]) : idHexs(face->neighbor_hs[0]);
			return result;
		}
	}

	template<typename V, typename E, typename F, typename H>
	E* topoM< V, E, F, H>::e_parallel_e_in_face(F* face, E* e)
	{
		int eIndex = face->edgeIndex(e->id());
		if (eIndex==-1)
		{
			return NULL;
		}
		E* resultE = idEdges(face->es[face->eparalle[eIndex][0]]);
		return resultE;
	}
	
	template<typename V, typename E, typename F, typename H>
	std::vector<E*> topoM< V, E, F, H>::v_adj_e_in_hex(H* hex, V* v)
	{
		std::vector<E*> results;
		int vIndex = hex->vertexIndex(v->id());
		if (vIndex == -1)
		{
			return results;
		}
		for (int i = 0; i < 3; i++)
		{
			E* pE = idEdges(hex->es[hex->vadje[vIndex][i]]);
			results.push_back(pE);
		}
		return results;
	}

	template<typename V, typename E, typename F, typename H>
	std::vector<F*> topoM< V, E, F, H>::v_adj_f_in_hex(H* hex, V* v)
	{
		std::vector<F*> results;
		int vIndex = hex->vertexIndex(v->id());
		if (vIndex == -1)
		{
			return results;
		}
		for (int i = 0; i < 3; i++)
		{
			F* pF = idFaces(hex->fs[hex->vadjf[vIndex][i]]);
			results.push_back(pF);
		}
		return results;
	}

	template<typename V, typename E, typename F, typename H>
	std::vector<E*> topoM< V, E, F, H>::e_parallel_e_in_hex(H* hex, E* e)
	{
		std::vector<E*> results;
		int eIndex = hex->edgeIndex(e->id());
		if (eIndex==-1)
		{
			return results;
		}
		for (int i = 0; i < 3; i++)
		{
			E* pE = idEdges(hex->es[hex->eparalle[eIndex][i]]);
			results.push_back(pE);
		}
		return results;
	}

	template<typename V, typename E, typename F, typename H>
	std::vector<F*> topoM< V, E, F, H>::e_adj_f_in_hex(H* hex, E* e)
	{
		std::vector<F*> results;
		int eIndex = hex->edgeIndex(e->id());
		if (eIndex == -1)
		{
			return results;
		}
		for (int i = 0; i < 2; i++)
		{
			F* pF = idFaces(hex->fs[hex->eadjf[eIndex][i]]);
			results.push_back(pF);
		}
		return results;
	}

	template<typename V, typename E, typename F, typename H>
	F* topoM< V, E, F, H>::f_parallel_f_in_hex(H* hex ,F* face)
	{
		int fIndex = hex->faceIndex(face->id());
		if (fIndex==-1)
		{
			return NULL;
		}
		F* resultF = idFaces(hex->fs[hex->fparallelf[fIndex][0]]);
		return resultF;
	}

	template<typename V, typename E, typename F, typename H>
	std::vector<E*> topoM< V, E, F, H>::f_adj_e_in_hex(H* hex, F* face) //edited
	{
		std::vector<E*> results;
		int fIndex = hex->faceIndex(face->id());
		if (fIndex == -1)
		{
			return NULL;
		}
		for (int i = 0; i < 4; i++)  // 使用固定值4，因为每个面有4条边
		{
			E* resultE = idEdges(hex->es[hex->fadje[fIndex][i]]);
			results.push_back(resultE);
		}	
		return results;
	}

	template<typename V, typename E, typename F, typename H>
	void topoM< V, E, F, H>::faceNormal()
	{
		for(std::list<F*>::iterator fite = fs.begin();fite!=fs.end();fite++ )
		{
			F* f = *fite;
			std::vector<CPoint> vp;
			for (int i = 0; i < 4; i++)
			{
				V* fv = idVertices(f->vs[i]);
				if (fv == NULL)
				{
					std::cout << "can't find v id: " << f->vs[i] << std::endl;
					std::cout << "can't find v of neighbor face id: " << f->id() << std::endl;
				}
				
				
				/*std::cout << "fv id: " << fv->id() << std::endl;
				std::cout << "fv position: " << fv->position() << std::endl;*/
				vp.push_back(fv->position());
			}
			CPoint n;
			for (int i = 0; i < 4; i++)
			{
				n+= (vp[(i+1)%4] - vp[i]) ^ (vp[(i+3)%4] - vp[i]);
			}			
			n = n*0.25;
			n = n / n.norm();
			f->normal() = n;	
		}
	}

	template<typename V, typename E, typename F, typename H>
	void topoM< V, E, F, H>::vertexNormal()
	{

		for (std::list<V*>::iterator vite = vs.begin();vite!=vs.end();vite++ )
		{
			V* v = *vite;
			if (!v->boundary()) continue;
			CPoint normal;
			
			for (int vfIndex = 0; vfIndex < v->neighbor_fs.size(); vfIndex++)
			{
				F* f = idFaces(v->neighbor_fs[vfIndex]);
				if (!f->boundary()) continue;
				
				CPoint fn = f->normal();
				normal = fn + normal;
				
			}
			normal /= normal.norm();
			v->normal() = normal;
		}
	}

	template<typename V, typename E, typename F, typename H>
	void topoM< V, E, F, H>::computeNormal()
	{
		faceNormal();
		vertexNormal();
	}

	template<typename V, typename E, typename F, typename H>
	void topoM< V, E, F, H>::normalize()
	{
		CPoint vmax(-1e+10, -1e+10, -1e+10);
		CPoint vmin(1e+10, 1e+10, 1e+10);

		for (std::list<V*>::iterator vite = vs.begin(); vite != vs.end(); vite++)
		{
			V* v = *vite;
			
			CPoint p = v->position();
			for (int k = 0; k < 3; k++)
			{
				vmax[k] = (vmax[k] > p[k]) ? vmax[k] : p[k];
				vmin[k] = (vmin[k] < p[k]) ? vmin[k] : p[k];
			}
		}
		
		CPoint center = (vmax + vmin) / 2.0;
		double d = 0;
		for (std::list<V*>::iterator vite = vs.begin(); vite != vs.end(); vite++)
		{
			V* pV = *vite;

			CPoint p = pV->position();
			p = p - center;
			pV->pre_position() = pV->position();
			pV->position() = p;

			for (int k = 0; k < 3; k++)
			{
				if (fabs(p[k]) > d) d = fabs(p[k]);
			}
		}
		
		for (std::list<V*>::iterator vite = vs.begin(); vite != vs.end(); vite++)
		{
			V* pV = *vite;
			CPoint p = pV->position();
			p /= d;
			pV->position() = p;
		}
		
	}

	template<typename V, typename E, typename F, typename H>
	std::vector<E*> topoM< V, E, F, H>::mark_singularity()
	{
		std::vector<E*>singularEdges;
		//mark the vertex singularity
		//clear all singualrity
		for (std::list<V*>::iterator vite = vs.begin(); vite != vs.end(); vite++)
		{
			V* v = *vite;
			v->singularity() = false;
		}
		for (std::list<E*>::iterator eite = es.begin(); eite != es.end(); eite++)
		{
			E* e = *eite;
			e->singularity() = false;
		}	

		//mark the edge singularity
		for (std::list<E*>::iterator eite = es.begin();eite!=es.end();eite++)
		{
			E* e = *eite;
			
			if (e->boundary())
			{
				if (e->neighbor_hs.size() != 2)
				{
					e->singularity() = true;
					singularEdges.push_back(e);
				}
				else
				{
					e->singularity() = false;
				}
			}
			else
			{
				if (e->neighbor_hs.size() != 4)
				{
					e->singularity() = true;
					singularEdges.push_back(e);
				}
				else
				{
					e->singularity() = false;
				}
			}
		}

		//mark the vertex singularity
		for (std::list<V*>::iterator vite = vs.begin(); vite != vs.end(); vite++)
		{
			V* v = *vite;
			for (int veIndex = 0; veIndex < v->neighbor_es.size(); veIndex++)
			{
				E* e = idEdges(v->neighbor_es[veIndex]);
				if (e->singularity())
				{
					v->singularity() = true;
					break;
				}
			}			
		}

		//mark singularity node
		mark_singularity_node();
		return singularEdges;
	}

	template<typename V, typename E, typename F, typename H>
	void topoM< V, E, F, H>::mark_singularity_node()
	{
		for (std::list<V*>::iterator vite = vs.begin(); vite != vs.end(); vite++)
		{
			V* v = *vite;
			if (!v->singularity()) continue;
			
			int count = 0;
			for (int veIndex = 0; veIndex < v->neighbor_es.size(); veIndex++)
			{
				E* e = idEdges(v->neighbor_es[veIndex]);
				if (e->singularity()) count++;				
			}
			if (count>2||count<2)
			{
				v->singularity_node() = true;
			}
		}

	}

	template<typename V, typename E, typename F, typename H>
	void topoM< V, E, F, H>::push_singularity(E* e)
	{
		std::vector<E*>::iterator ite = find(singularities.begin(), singularities.end(), e);
		if (ite==singularities.end())
		{
			singularities.push_back(e);
		}
	}

	template<typename V, typename E, typename F, typename H>
	void topoM< V, E, F, H>::delete_singularity(E* e)
	{
		std::vector<E*>::iterator ite = find(singularities.begin(), singularities.end(), e);
		if (ite != singularities.end())
		{
			singularities.erase(ite);
		}
	}

	template<typename V, typename E, typename F, typename H>
	void  topoM< V, E, F, H>::mark_boundary()
	{
		//clear all the mark
		for (std::list<V*>::iterator vite = vs.begin(); vite != vs.end(); vite++)
		{
			V* v = *vite;
			v->boundary() = false;
		}
		for (std::list<E*>::iterator eite = es.begin(); eite != es.end(); eite++)
		{
			E* e = *eite;
			e->boundary() = false;
		}
		for (std::list<F*>::iterator fite = fs.begin(); fite != fs.end(); fite++)
		{
			F* f = *fite;			
			f->boundary() = false;
			if (f->neighbor_hs.size() == 1)
			{				
				f->boundary() = true;
				for (int eIndex = 0; eIndex < 4; eIndex++)
				{
					E* e = idEdges(f->es[eIndex]);		
					e->boundary() = true;
				}
				for (int vIndex = 0; vIndex < 4; vIndex++)
				{
					V* v = idVertices(f->vs[vIndex]);
					v->boundary() = true;
				}
				H* hex = idHexs(f->neighbor_hs[0]);		
				hex->boundary() = true;
			}			
		}
	}

	template<typename V, typename E, typename F, typename H>
	void topoM<V, E, F, H>::construct_hex(H* h, int hid, int* v)
	{		
		h->id() = hid;
		
		//set the vertices of hex
		for (int i = 0; i < 8; i++)
		{
			V* hexV = this->idVertices(v[i]);
			if (hexV ==NULL)
			{
				std::cout << "can't find v:" << std::endl;
				std::cout << "V->id(): " << v[i] << std::endl;
			}
			h->vs.push_back(v[i]);
			hexV->neighbor_hs.push_back(hid);
			
		}
		
		/*set face and edge*/
		int order[6][4] = { { 3,2,1,0 }, { 0,1,5,4 }, { 4,5,6,7 }, { 7,6,2,3 },{1,2,6,5},{3,0,4,7} };
		F* hexFs[6];
		E* hexEs[12];
		for (int i = 0; i < 6; i++)
		{
			std::vector<int> fv(4);
			for (int j = 0; j < 4; j++)
			{
				fv[j] = v[order[i][j]];
			}
			
			//�?��?���?������
			V* facev = idVertices(fv[0]);
			
			bool construct_new_face = true;
			for (int adjFIndex = 0; adjFIndex < facev->neighbor_fs.size(); adjFIndex++)
			{
				F* adjF = idFaces(facev->neighbor_fs[adjFIndex]);

				if (adjF->is_v_equal(fv))
				{				
					construct_new_face = false;
					h->fs.push_back(adjF->id());
					hexFs[i] = adjF;//���?����?�face����fs��
					adjF->neighbor_hs.push_back(h->id());
					
					//�����?����?�����?����?�?���hexEs��
					if (i<4)
					{
						int edgeOrder[3][2] = { {0,1},{3,0},{1,2} };
						//���??hexEs��
						for (int feIndex = 0; feIndex < 3; feIndex++)
						{
							V* ev1 = idVertices(fv[edgeOrder[feIndex][0]]);
							V* ev2 = idVertices(fv[edgeOrder[feIndex][1]]);
							E* adjFE = VerticesEdge(ev1, ev2);
							hexEs[i * 3 + feIndex] = adjFE;

						}
						
					}
					break;
				}
			}
			
			if (construct_new_face)
			{
				F* newf = new F();
				m_maxFaceId++;
				newf->id() = m_maxFaceId;
				fs.push_back(newf);
				m_map_faces.insert(std::pair<int, F*>(newf->id(), newf));
				//��?��������?;
				for (int fvIndex = 0; fvIndex < 4; fvIndex++)
				{
					newf->vs.push_back(fv[fvIndex]);
					V* pV = this->idVertices(fv[fvIndex]);
					pV->neighbor_fs.push_back(newf->id());
				}
				//��?�����?�?
				newf->neighbor_hs.push_back(h->id());
				hexFs[i] = newf;
				h->fs.push_back(newf->id());
							
				if (i<4)//?�?��?����?��������	
				{
					//������
					//�?��?�������				
					int edgeOrder[3][2] = { {0,1},{3,0},{1,2} };
					for (int eIndex = 0; eIndex < 3; eIndex++)
					{
						//int ev[2];
						std::vector<int> ev(2);
						ev[0] = fv[edgeOrder[eIndex][0]];
						ev[1] = fv[edgeOrder[eIndex][1]];				

						V* ev1 = this->idVertices(ev[0]);
						bool construct_edge = true;
						for (int adjEIndex = 0; adjEIndex < ev1->neighbor_es.size(); adjEIndex++)
						{
							E* adjE = this->idEdges(ev1->neighbor_es[adjEIndex]);
							if (adjE->is_v_equal(ev))
							{
								construct_edge = false;
								hexEs[i * 3 + eIndex] = adjE;								
								break;
							}
						}
						
						if (construct_edge==true)
						{
							E* newE = new E();
							m_maxEdgeId++;
							newE->id() = m_maxEdgeId;
							es.push_back(newE);
							m_map_edges.insert(std::pair<int, E*>(newE->id(), newE));
							
							
							if (ev[0]<ev[1])
							{
								newE->vs.push_back(ev[0]);
								newE->vs.push_back(ev[1]);
							}
							else
							{
								newE->vs.push_back(ev[1]);
								newE->vs.push_back(ev[0]);
							}
													
							//������???��?
							for (int vIndex = 0; vIndex < 2; vIndex++)
							{
								V* Edgev = this->idVertices(ev[vIndex]);
								Edgev->neighbor_es.push_back(newE->id());
							}					
							hexEs[i * 3 + eIndex] = newE;
							
						}
					}
				}
				
			}			
		}
		

		//�������?����??�?		
		int feOrder[6][4] = { {0,2,3,1},{3, 5 ,6 ,4},{6,8,9,7},{9,11,0,10}, {2,11,8,5},{1,4,7,10} };
		for (int hexFIndex = 0; hexFIndex < 6; hexFIndex++)
		{		
			F* hexf = hexFs[hexFIndex];	
			if (hexf->es.size()==4)//?���������������
			{
				continue;
			}				
			for (int eIndex = 0; eIndex < 4; eIndex++)
			{
				E* fe = hexEs[feOrder[hexFIndex][eIndex]];
				hexf->es.push_back(fe->id());
				fe->neighbor_fs.push_back(hexf->id());
				
			}
		}
		//�����?�������?��?
		for (int hexEIndex = 0; hexEIndex < 12; hexEIndex++)
		{
			hexEs[hexEIndex]->neighbor_hs.push_back(h->id());
			h->es.push_back(hexEs[hexEIndex]->id());
		}

	}

	template<typename V, typename E, typename F, typename H>
	void topoM<V, E, F, H>::construct_hex(H* h, int hid, std::vector<int> v)
	{

		h->id() = hid;
		//set the vertices of hex
		for (int i = 0; i < 8; i++)
		{
			V* hexV = this->idVertices(v[i]);
			if (hexV == NULL)
			{
				std::cout << "can't find v:" << std::endl;
				std::cout << "V->id(): " << v[i] << std::endl;
			}
			h->vs.push_back(v[i]);
			hexV->neighbor_hs.push_back(hid);

		}
		
		/*set face and edge*/
		int order[6][4] = { { 3,2,1,0 }, { 0,1,5,4 }, { 4,5,6,7 }, { 7,6,2,3 },{1,2,6,5},{3,0,4,7} };
		F* hexFs[6];
		E* hexEs[12];
		for (int i = 0; i < 6; i++)
		{
			std::vector<int> fv(4);
			for (int j = 0; j < 4; j++)
			{
				fv[j] = v[order[i][j]];
			}

			//�?��?���?������
			V* facev = idVertices(fv[0]);

			bool construct_new_face = true;
			for (int adjFIndex = 0; adjFIndex < facev->neighbor_fs.size(); adjFIndex++)
			{
				F* adjF = idFaces(facev->neighbor_fs[adjFIndex]);

				if (adjF->is_v_equal(fv))
				{
					construct_new_face = false;
					h->fs.push_back(adjF->id());
					hexFs[i] = adjF;//���?����?�face����fs��
					adjF->neighbor_hs.push_back(h->id());

					//�����?����?�����?����?�?���hexEs��
					if (i < 4)
					{
						int edgeOrder[3][2] = { {0,1},{3,0},{1,2} };
						//���??hexEs��	
						
						for (int feIndex = 0; feIndex < 3; feIndex++)
						{
							V* ev1 = idVertices(fv[edgeOrder[feIndex][0]]);
							V* ev2 = idVertices(fv[edgeOrder[feIndex][1]]);
							E* adjFE = VerticesEdge(ev1, ev2);
							hexEs[i * 3 + feIndex] = adjFE;

						}
						//std::vector<int> ev(2);
						//ev[0] = fv[edgeOrder[2][0]];
						//ev[1] = fv[edgeOrder[2][1]];

						////���?����?����???���hexEs��

						//for (int eIndex = 0; eIndex < 4; eIndex++)
						//{
						//	/*std::cout << "fv->id(): " << facev->id() << std::endl;
						//	std::cout << "adjF->id(): " << adjF->id() << std::endl;
						//	std::cout << "adjF->es[eIndex]: " << adjF->es[eIndex] << std::endl;
						//	std::cout << "---------------------------------------"<< std::endl;*/
						//	E* adjFE = idEdges(adjF->es[eIndex]);
						//	
						//	if (adjFE->is_v_equal(ev))
						//	{
						//		hexEs[i * 3 + 2] = adjFE;
						//		adjFE = idEdges(adjF->es[(eIndex + 1) % 4]);
						//		hexEs[i * 3 + 0] = adjFE;
						//		adjFE = idEdges(adjF->es[(eIndex + 2) % 4]);
						//		hexEs[i * 3 + 1] = adjFE;
						//		break;
						//	}
						//}
					}
					break;
				}
			}

			if (construct_new_face)
			{
				F* newf = new F();
				m_maxFaceId++;
				newf->id() = m_maxFaceId;
				fs.push_back(newf);
				m_map_faces.insert(std::pair<int, F*>(newf->id(), newf));
				//��?��������?;
				for (int fvIndex = 0; fvIndex < 4; fvIndex++)
				{
					newf->vs.push_back(fv[fvIndex]);
					V* pV = this->idVertices(fv[fvIndex]);
					pV->neighbor_fs.push_back(newf->id());
				}
				//��?�����?�?
				newf->neighbor_hs.push_back(h->id());
				hexFs[i] = newf;
				h->fs.push_back(newf->id());

				if (i < 4)//?�?��?����?��������	
				{
					//������
					//�?��?�������				
					int edgeOrder[3][2] = { {0,1},{3,0},{1,2} };
					for (int eIndex = 0; eIndex < 3; eIndex++)
					{
						//int ev[2];
						std::vector<int> ev(2);
						ev[0] = fv[edgeOrder[eIndex][0]];
						ev[1] = fv[edgeOrder[eIndex][1]];

						V* ev1 = this->idVertices(ev[0]);
						bool construct_edge = true;
						for (int adjEIndex = 0; adjEIndex < ev1->neighbor_es.size(); adjEIndex++)
						{
							E* adjE = this->idEdges(ev1->neighbor_es[adjEIndex]);
							if (adjE->is_v_equal(ev))
							{
								construct_edge = false;
								hexEs[i * 3 + eIndex] = adjE;
								break;
							}
						}

						if (construct_edge == true)
						{
							E* newE = new E();
							m_maxEdgeId++;
							newE->id() = m_maxEdgeId;
							es.push_back(newE);
							m_map_edges.insert(std::pair<int, E*>(newE->id(), newE));


							if (ev[0] < ev[1])
							{
								newE->vs.push_back(ev[0]);
								newE->vs.push_back(ev[1]);
							}
							else
							{
								newE->vs.push_back(ev[1]);
								newE->vs.push_back(ev[0]);
							}

							//������???��?
							for (int vIndex = 0; vIndex < 2; vIndex++)
							{
								V* Edgev = this->idVertices(ev[vIndex]);
								Edgev->neighbor_es.push_back(newE->id());
							}
							hexEs[i * 3 + eIndex] = newE;

						}
					}
				}

			}
		}
		//�������?����??�?		
		int feOrder[6][4] = { {0,2,3,1},{3, 5 ,6 ,4},{6,8,9,7},{9,11,0,10}, {2,11,8,5},{1,4,7,10} };
		for (int hexFIndex = 0; hexFIndex < 6; hexFIndex++)
		{
			F* hexf = hexFs[hexFIndex];
			if (hexf->es.size() == 4)//?���������������
			{
				continue;
			}
			for (int eIndex = 0; eIndex < 4; eIndex++)
			{
				E* fe = hexEs[feOrder[hexFIndex][eIndex]];
				hexf->es.push_back(fe->id());
				fe->neighbor_fs.push_back(hexf->id());
				
			}
		}
		
		//�����?�������?��?
		for (int hexEIndex = 0; hexEIndex < 12; hexEIndex++)
		{
			hexEs[hexEIndex]->neighbor_hs.push_back(h->id());
			h->es.push_back(hexEs[hexEIndex]->id());
			
		}
		std::cout << "-----------------------------" << std::endl;

	}

	// Structure to hold temporary data during single-pass loading
	struct TempHexData {
		int id;
		std::vector<int> v_ids;
		std::string str_attr;
	};
	struct TempEdgeData {
		int v0_id, v1_id;
		std::string str_attr;
	};
	struct TempFaceData {
		std::vector<int> v_ids;
		std::string str_attr;
	};
	struct TempCornerData {
		int v_id;
		std::string str_attr;
	};

	template<typename V, typename E, typename F, typename H>
	void topoM<V,E,F,H>::load_Qhex(const char* input)
	{
		std::fstream is(input, std::fstream::in);

		if (is.fail())
		{
			fprintf(stderr, "Error in opening file %s\n", input);
			return;
		}

		char buffer[MAX_LINE];
		
		m_nVertices = 0;

		m_maxVertexId = 0;		
		m_maxEdgeId = 0;
		m_maxFaceId = 0;
		m_maxHexId = 0;
		m_nSharpEdges = 0;
		m_nCorners = 0;

		while (!is.eof())
		{
			is.getline(buffer, MAX_LINE);
			std::string line(buffer);
			line = strutil::trim(line);
			strutil::Tokenizer stokenizer(line, " \r\n");

			stokenizer.nextToken();
			std::string token = stokenizer.getToken();

			if (token == "Vertex") m_nVertices++;
			if (token == "Hex") m_nHexs++;
			if (token == "Face") m_nFaces++;
			if (token == "Edge") m_nSharpEdges++;
			if (token == "corner") m_nCorners++;
		}
		

		is.clear();              // forget we hit the end of file
		is.seekg(0, std::ios::beg);   // move to the start of the file
		
		for (int i = 0; i < m_nVertices && is.getline(buffer, MAX_LINE); i++)
		{			
			std::string line(buffer);
			line = strutil::trim(line);
			strutil::Tokenizer stokenizer(line, " \r\n");
			stokenizer.nextToken();
			std::string token = stokenizer.getToken();

			if (token != "Vertex")
			{
				fprintf(stderr, "File Format Error\r\n");
				return;
			}
			stokenizer.nextToken();
			token = stokenizer.getToken();
			int vid = strutil::parseString<int>(token);
			m_maxVertexId = vid > m_maxVertexId ? vid : m_maxVertexId;

			CPoint p;
			for (int k = 0; k < 3; k++)
			{
				stokenizer.nextToken();
				std::string token = stokenizer.getToken();
				p[k] = strutil::parseString<float>(token);
			}
			V* v = new V();
			v->id() = vid;
			v->position() = p;
			vs.push_back(v);
			m_map_vertices.insert(std::pair<int, V*>(vid, v));
			if (!stokenizer.nextToken("\t\r\n")) continue;
			token = stokenizer.getToken();

			int sp = (int)token.find("{");
			int ep = (int)token.find("}");

			if (sp >= 0 && ep >= 0)
			{
				v->string() = token.substr(sp + 1, ep - sp - 1);
			}
		}
		
		//read in Hex 		
		for (int id = 0; id < m_nHexs && is.getline(buffer, MAX_LINE); id++)
		{
			
			int vid[8];
			std::string line(buffer);
			line = strutil::trim(line);
			strutil::Tokenizer stokenizer(line, " \r\n");

			stokenizer.nextToken();
			std::string token = stokenizer.getToken();

			if (token != "Hex")
			{
				fprintf(stderr, "File Format Error\r\n");
				return;
			}
			//skip the first "4" in the line
			stokenizer.nextToken();
			token = stokenizer.getToken();
			int hid = strutil::parseString<int>(token);
			m_maxHexId = hid > m_maxHexId ? hid : m_maxHexId;
			for (int k = 0; k < 8; k++)
			{
				stokenizer.nextToken();
				std::string token = stokenizer.getToken();
				vid[k] = strutil::parseString<int>(token);
			}

			H* pHex = new H();
			hs.push_back(pHex);
			m_map_hexs.insert(std::pair<int, H*>(hid, pHex));
			construct_hex(pHex, hid, vid);

			// read in string
			if (!stokenizer.nextToken("\t\r\n")) continue;
			token = stokenizer.getToken();

			int sp = (int)token.find("{");
			int ep = (int)token.find("}");

			if (sp >= 0 && ep >= 0)
			{
				pHex->string() = token.substr(sp + 1, ep - sp - 1);
			}

		}
		
		//read in edge
		for (int id = 0; id < m_nSharpEdges && is.getline(buffer, MAX_LINE); id++)
		{
			std::string line(buffer);
			line = strutil::trim(line);
			strutil::Tokenizer stokenizer(line, " \r\n");

			stokenizer.nextToken();
			std::string token = stokenizer.getToken();

			if (token != "Edge")
			{
				fprintf(stderr, "File Format Error\r\n");
				return;
			}
			stokenizer.nextToken();
			token = stokenizer.getToken();
			int id0 = strutil::parseString<int>(token);

			stokenizer.nextToken();
			token = stokenizer.getToken();
			int id1 = strutil::parseString<int>(token);

			V* v0 = idVertices(id0);
			V* v1 = idVertices(id1);

			E* edge = VerticesEdge(v0, v1);
			if (edge==NULL)
			{
				std::cout << "edge is not exist: v1:" << id0 << " v2 " << id1 << std::endl;
				continue;
			}

			if (!stokenizer.nextToken("\t\r\n")) continue;
			token = stokenizer.getToken();

			int sp = (int)token.find("{");
			int ep = (int)token.find("}");

			if (sp >= 0 && ep >= 0)
			{
				edge->string() = token.substr(sp + 1, ep - sp - 1);
			}
		}

		//read in face
		for (int id = 0; id < m_nFaces && is.getline(buffer, MAX_LINE); id++)
		{
			std::string line(buffer);
			line = strutil::trim(line);
			strutil::Tokenizer stokenizer(line, " \r\n");

			stokenizer.nextToken();
			std::string token = stokenizer.getToken();

			if (token != "Face")
			{
				fprintf(stderr, "File Format Error\r\n");
				return;
			}
			std::vector<int> vid;
			for (int k = 0; k < 4; k++)
			{
				stokenizer.nextToken();
				std::string token = stokenizer.getToken();
				int v_id = strutil::parseString<int>(token);
				vid.push_back(v_id);
			}
			F* face = VerticesFace(vid);
			if (face == NULL)
			{
				std::cout << "face is not exist: ";
				for (int vIndex = 0; vIndex < 4; vIndex++)
				{
					std::cout << " " << vid[vIndex];
				}
				std::cout << std::endl;
				continue;
			}

			if (!stokenizer.nextToken("\t\r\n")) continue;
			token = stokenizer.getToken();

			int sp = (int)token.find("{");
			int ep = (int)token.find("}");

			if (sp >= 0 && ep >= 0)
			{
				face->string() = token.substr(sp + 1, ep - sp - 1);

			}

		}
		
		/*read in corner*/
		for (int id = 0; id < m_nCorners && is.getline(buffer, MAX_LINE); id++)
		{
			std::string line(buffer);
			line = strutil::trim(line);
			strutil::Tokenizer stokenizer(line, " \r\n");

			stokenizer.nextToken();
			std::string token = stokenizer.getToken();
			if (token != "Corner")
			{
				fprintf(stderr, "File Format Error\r\n");
				return;
			}
			stokenizer.nextToken();
			token = stokenizer.getToken();
			int vid = strutil::parseString<int>(token);
			V* v = idVertices(vid);

			if (!stokenizer.nextToken("\t\r\n")) continue;
			token = stokenizer.getToken();

			int sp = (int)token.find("{");
			int ep = (int)token.find("}");

			if (sp >= 0 && ep >= 0)
			{
				v->string() = token.substr(sp + 1, ep - sp - 1);
			}
		}

		is.close();

		compute_features();
		// print_hexs_map();
		// print_faces_map();
		// print_edges_map();
		// print_vertices_map();
	}

	template<typename V, typename E, typename F, typename H>
	void topoM<V, E, F, H>::write_Qhex(const char* output)
	{
		// 创建日志函数
		auto logFunc = [](const std::string& msg) {
			// 检查是否存在topoMeshLog函数
			if (HMeshLib::getTopoMeshLog()) {
				HMeshLib::topoMeshLog(msg);
			}
			// 检查是否存在sheet_operation日志函数
			else if (HMeshLib::getSheetLog()) {
				HMeshLib::log(msg);
			}
		};

		try {
			//logFunc("Starting write_Qhex to file: " + std::string(output));
			
			// 准备要写入的数据
			//logFunc("Preparing data for output - converting to string format");
			for (std::list<V*>::iterator vite = vs.begin(); vite != vs.end(); vite++)
			{
				V* v = *vite;
				v->_to_string();
			}
			//logFunc("Vertices data prepared: " + std::to_string(vs.size()) + " vertices");
			
			for (std::list<E*>::iterator eite = es.begin(); eite != es.end(); eite++)
			{
				E* e = *eite;
				e->_to_string();
			}
			//logFunc("Edges data prepared: " + std::to_string(es.size()) + " edges");
			
			for (std::list<H*>::iterator hite = hs.begin(); hite != hs.end(); hite++)
			{
				H* h = *hite;
				h->_to_string();
			}
			//logFunc("Hexes data prepared: " + std::to_string(hs.size()) + " hexes");

			// 打开输出文件
			//logFunc("Opening output file: " + std::string(output));
			std::fstream _os(output, std::fstream::out);

			if (_os.fail())
			{
				//logFunc("ERROR: Failed to open output file: " + std::string(output));
				fprintf(stderr, "Error in opening file %s\n", output);
				return;
			}
			
			// 输出顶点数据
			//logFunc("Writing vertices data to file...");
			int vertexCount = 0;
			for (std::list<V*>::iterator vite = vs.begin(); vite != vs.end();vite++)
			{
				V* v = *vite;
				CPoint p = v->position();
				_os << "Vertex " << v->id() << " " << p[0] << " " << p[1] << " " << p[2];
				if (v->string().size() > 0)
				{
					_os << " " << "{" << v->string() << "}";
				}
				_os << std::endl;
				vertexCount++;
				
				// 每写入100个顶点记录一次日志
				// if (vertexCount % 100 == 0) {
				// 	logFunc("  Wrote " + std::to_string(vertexCount) + " vertices...");
				// }
			}
			//logFunc("Completed writing " + std::to_string(vertexCount) + " vertices");
			
			// 输出六面体数据
			//logFunc("Writing hex data to file...");
			int hexCount = 0;
			for (std::list<H*>::iterator hite = hs.begin(); hite != hs.end(); hite++)
			{
				H* h = *hite;
				_os << "Hex " << h->id();
				for (int i = 0; i < 8; i++)
				{
					_os << " " << h->vs[i];
				}
				if (h->string().size() > 0)
				{
					_os << " " << "{" << h->string() << "}";
				}
				_os << std::endl;
				hexCount++;
				
				// 每写入100个六面体记录一次日志
				// if (hexCount % 100 == 0) {
				// 	logFunc("  Wrote " + std::to_string(hexCount) + " hexes...");
				// }
			}
			//logFunc("Completed writing " + std::to_string(hexCount) + " hexes");

			// 输出边数据
			//logFunc("Writing edge data to file...");
			int edgeCount = 0;
			for (std::list<E*>::iterator eite = es.begin(); eite != es.end(); eite++)
			{
				E* e = *eite;
				if (e->string().size() > 0)
				{
					_os << "Edge " << e->vs[0] << " " << e->vs[1] << " ";
					_os << "{" << e->string() << "}" << std::endl;
					edgeCount++;
				}
			}
			//logFunc("Completed writing " + std::to_string(edgeCount) + " edges with string attributes");

			// 关闭文件
			_os.close();
			//logFunc("File closed successfully: " + std::string(output));
		}
		catch (const std::exception& e) {
			logFunc("Exception occurred during write_Qhex: " + std::string(e.what()));
		}
		catch (...) {
			logFunc("Unknown exception occurred during write_Qhex");
		}
	}

	/*	
	* delete element
	*/
	template<typename V, typename E, typename F, typename H>
	void topoM<V, E, F, H>::delete_vertex(V* v)
	{
		m_map_vertices.erase(v->id());
		vs.remove(v);
		//delete v;
	}
	template<typename V, typename E, typename F, typename H>
	void topoM<V, E, F, H>::delete_edge(E* e)
	{
		for (int i = 0; i < 2; i++)
		{
			if (m_map_vertices.count(e->vs[i])>0)
			{
				V* v = idVertices(e->vs[i]);
				v->delete_neighbor_e(e->id());
			}		
		}	
		m_map_edges.erase(e->id());
		es.remove(e);
		//delete e;
	}
	template<typename V, typename E, typename F, typename H>
	void topoM<V, E, F, H>::delete_face(F* f)
	{
		for (int i = 0; i < 4; i++)
		{
			if (m_map_vertices.count(f->vs[i]) > 0)
			{
				V* v = idVertices(f->vs[i]);
				v->delete_neighbor_f(f->id());
			}
			if (m_map_edges.count(f->es[i]) > 0)
			{
				E* e = idEdges(f->es[i]);
				e->delete_neighbor_f(f->id());
			}
		}
		m_map_faces.erase(f->id());
		fs.remove(f);
		//delete f;
	}
	template<typename V, typename E, typename F, typename H>
	void topoM<V, E, F, H>::delete_hex(H* h)
	{
		for (int i = 0; i < 8; i++)
		{
			if (m_map_vertices.count(h->vs[i]) > 0)
			{
				V* v = idVertices(h->vs[i]);
				v->delete_neighbor_h(h->id());
				if (v->neighbor_hs.size() == 0)
				{
					delete_vertex(v);
				}
			}
		}
		for (int i = 0; i < 12; i++)
		{
			if (m_map_edges.count(h->es[i]) > 0)
			{
				E* e = idEdges(h->es[i]);
				e->delete_neighbor_h(h->id());
				if (e->neighbor_hs.size() == 0)
				{
					delete_edge(e);
				}
			}
		}
		for (int i = 0; i < 6; i++)
		{
			if (m_map_faces.count(h->fs[i]) > 0)
			{
				F* f = idFaces(h->fs[i]);
				f->delete_neighbor_h(h->id());
				if (f->neighbor_hs.size() == 0)
				{
					delete_face(f);
				}
			}
		}
		m_map_hexs.erase(h->id());	
		hs.remove(h);		
		delete h;
	}

	/*change element in hex*/
	template<typename V, typename E, typename F, typename H>
	bool topoM<V, E, F, H>::change_Hex_V(H* h, V* source, V* target)
	{
		//Is source in hex?
		int vIndex = h->vertexIndex(source->id());
		if (vIndex==-1)
		{
			return false;
		}
		//change the v in hex
		h->vs[vIndex] = target->id();
		//change the source neighbor relation in hex
		for (int i = 0; i < 3; i++)
		{
			int adjeId = h->es[h->vadje[vIndex][i]];
			E* adje = idEdges(adjeId);
			target->push_back_neighbor_e(adjeId);
			source->delete_neighbor_e(adjeId);
			int adjevIndex = adje->vertexIndex(source->id());
			if (adjevIndex!=-1)
			{
				adje->vs[adjevIndex] = target->id();
			}

			int adjfId = h->fs[h->vadjf[vIndex][i]];
			
			F* adjf = idFaces(adjfId);
			target->push_back_neighbor_f(adjfId);
			source->delete_neighbor_f(adjfId);
			int adjfvIndex = adjf->vertexIndex(source->id());
			if (adjfvIndex!=-1)
			{
				adjf->vs[adjfvIndex] = target->id();
			}
		}
		//add the target neighbor hex
		target->push_back_neighbor_h(h->id());
		source->delete_neighbor_h(h->id());
		return true;
	}
	
	template<typename V, typename E, typename F, typename H>
	bool topoM<V, E, F, H>::change_Hex_E(H* h, E* source, E* target)
	{
		//Is source in hex?
		int eIndex = h->edgeIndex(source->id());
		if (eIndex==-1)
		{
			return false;
		}
		//change the e in hex
		h->es[eIndex] = target->id();
		//change the source neighbor relation in hex
		for (int i = 0; i < 2; i++)
		{
			int adjfId = h->fs[h->eadjf[eIndex][i]];
			F* adjf = idFaces(adjfId);
			
			target->push_back_neighbor_f(adjfId);
			source->delete_neighbor_f(adjfId);
			int adjfeIndex = adjf->edgeIndex(source->id());
			if (adjfeIndex!=-1)
			{
				adjf->es[adjfeIndex] = target->id();
			}
		}		
		//add the target neighbor hex
		target->push_back_neighbor_h(h->id());
		source->delete_neighbor_h(h->id());
		return true;
	}
	
	template<typename V, typename E, typename F, typename H>
	bool topoM<V, E, F, H>::change_Hex_F(H* h, F* source, F* target)
	{
		//is source in hex?
		int fIndex = h->faceIndex(source->id());
		if (fIndex==-1)
		{
			return false;
		}
		//change the f in hex
		h->fs[fIndex] = target->id();
		//change the source neighbor relation in hex
		target->push_back_neighbor_h(h->id());
		source->delete_neighbor_h(h->id());
		return true;
	}

	/*change element in face*/
	template<typename V, typename E, typename F, typename H>
	bool topoM<V, E, F, H>::change_Face_V(F* f, V* source, V* target)
	{
		int vIndex = f->vertexIndex(source->id());
		if (vIndex == -1)
		{
			return false;
		}
		//change the v in hex
		f->vs[vIndex] = target->id();
		////change the source neighbor relation in hex
		//for (int i = 0; i < 3; i++)
		//{
		//	int adjeId = h->es[h->vadje[vIndex][i]];
		//	E* adje = idEdges(adjeId);
		//	target->push_back_neighbor_e(adjeId);
		//	source->delete_neighbor_e(adjeId);
		//	int adjevIndex = adje->vertexIndex(source->id());
		//	if (adjevIndex != -1)
		//	{
		//		adje->vs[adjevIndex] = target->id();
		//	}

		//	int adjfId = h->fs[h->vadjf[vIndex][i]];

		//	F* adjf = idFaces(adjfId);
		//	target->push_back_neighbor_f(adjfId);
		//	source->delete_neighbor_f(adjfId);
		//	int adjfvIndex = adjf->vertexIndex(source->id());
		//	if (adjfvIndex != -1)
		//	{
		//		adjf->vs[adjfvIndex] = target->id();
		//	}
		//}
		////add the target neighbor hex
		//target->push_back_neighbor_h(h->id());
		//source->delete_neighbor_h(h->id());
		return true;
	}

	/*change element in edge*/
	template<typename V, typename E, typename F, typename H>
	bool topoM<V, E, F, H>::change_Edge_V(E* e, V* source, V* target)
	{
		// 检查source是否是边e的一个端点
		int vIndex = e->vertexIndex(source->id());
		if (vIndex == -1)
		{
			return false;
		}
		
		// 修改边上的顶点
		e->vs[vIndex] = target->id();
		
		// 更新邻接关系
		target->push_back_neighbor_e(e->id());
		source->delete_neighbor_e(e->id());
		
		// 更新与这条边相关的面的顶点
		for (int i = 0; i < e->neighbor_fs.size(); i++)
		{
			int fid = e->neighbor_fs[i];
			F* f = idFaces(fid);
			if (f == NULL) continue;
			
			int fvIndex = f->vertexIndex(source->id());
			if (fvIndex != -1)
			{
				f->vs[fvIndex] = target->id();
				target->push_back_neighbor_f(fid);
				source->delete_neighbor_f(fid);
			}
		}
		
		return true;
	}
	

	/*revised normal*/
	template<typename V, typename E, typename F, typename H>
	void topoM<V, E, F, H>::revise_hex_face_normal(H* h)
	{
		// // 首先获取日志函数，如果在HMeshLib命名空间中存在
		// auto logFunc = [](const std::string& msg) {
		// 	// 检查是否存在//topoMeshLog函数
		// 	if (HMeshLib::gettopoMeshLog()) {
		// 		HMeshLib::topoMeshLog(msg);
		// 	}
		// 	// 检查是否存在sheet_operation日志函数
		// 	else if (HMeshLib::getSheetLog()) {
		// 		HMeshLib::log(msg);
		// 	}
		// };
		
		try {
			//logFunc("In revise_hex_face_normal function.");
			// print_hexs_map();
			// print_faces_map();
			// print_edges_map();
			// print_vertices_map();
			// Inside revise_hex_face_normal(H* h)
			//logFunc("Entering revise_hex_face_normal for hex ID: " + std::to_string(h->id()));
			
			// 先检查六面体的顶点是否都存在
			//logFunc("Checking hex vertices:");
			// for (int i = 0; i < h->vs.size(); i++) {
			// 	int vid = h->vs[i];
			// 	auto iter = m_map_vertices.find(vid);
			// 	if (iter == m_map_vertices.end()) {
			// 		logFunc("  Hex vertex " + std::to_string(vid) + " at index " + std::to_string(i) + " does not exist in vertex map!");
			// 	} else if (iter->second == nullptr) {
			// 		logFunc("  Hex vertex " + std::to_string(vid) + " at index " + std::to_string(i) + " exists in map but is nullptr!");
			// 	} else {
			// 		logFunc("  Hex vertex " + std::to_string(vid) + " at index " + std::to_string(i) + " exists and is valid.");
			// 	}
			// }
			
			// 打印与六面体相连的所有面的ID
			// logFunc("Connected faces:");
			// for (int i = 0; i < h->fs.size(); i++) {
			// 	logFunc("  Face ID: " + std::to_string(h->fs[i]) + " at index " + std::to_string(i));
			// }
			
			for (int fIndex = 0; fIndex < h->fs.size(); fIndex++)
			{
				int fid = h->fs[fIndex];
				//logFunc("  Processing face index " + std::to_string(fIndex) + ", face ID: " + std::to_string(fid));
				F* f = idFaces(fid);
				if (f == nullptr) {
					//logFunc("    Error: Face pointer is null!");
					continue; // 跳过这个面
				}
				//logFunc("    Face pointer retrieved successfully.");
				
				// 检查面的邻接六面体
				if (f->neighbor_hs.empty()) {
					//logFunc("    Warning: Face " + std::to_string(f->id()) + " has no neighbors!");
					continue;
				} else if (f->neighbor_hs[0] != h->id()) {
					//logFunc("    Skipping face " + std::to_string(f->id()) + " as h is not the first neighbor.");
					continue;
				}
				
				

				// 检查面的顶点
				//logFunc("    Face vertices:");
				// for (int i = 0; i < f->vs.size(); i++) {
				// 	int vid = f->vs[i];
				// 	auto iter = m_map_vertices.find(vid);
				// 	if (iter == m_map_vertices.end()) {
				// 		logFunc("      Face vertex " + std::to_string(vid) + " at index " + std::to_string(i) + " does not exist in vertex map!");
				// 	} else if (iter->second == nullptr) {
				// 		logFunc("      Face vertex " + std::to_string(vid) + " at index " + std::to_string(i) + " exists in map but is nullptr!");
				// 	}
				// }
				
				// 检查draw_order中的顶点
				// logFunc("    Checking draw_order vertices for face index " + std::to_string(fIndex) + ":");
				// for (int i = 0; i < 6/*draw_order[fIndex].size()*/; i++) {
				// 	int orderIndex = h->draw_order[fIndex][i];
				// 	if (orderIndex < 0 || orderIndex >= h->vs.size()) {
				// 		logFunc("      draw_order[" + std::to_string(fIndex) + "][" + std::to_string(i) + 
				// 			"] = " + std::to_string(orderIndex) + " is out of bounds for hex vertices array!");
				// 		continue;
				// 	}
				// 	int vid = h->vs[orderIndex];
				// 	auto iter = m_map_vertices.find(vid);
				// 	if (iter == m_map_vertices.end()) {
				// 		logFunc("      draw_order vertex " + std::to_string(vid) + " does not exist in vertex map!");
				// 	} else if (iter->second == nullptr) {
				// 		logFunc("      draw_order vertex " + std::to_string(vid) + " exists in map but is nullptr!");
				// 	}
				// }
				
				// 检查所有顶点是否存在
				// bool missingVertex = false;
				// for (int vIndex = 0; vIndex < 4; vIndex++) {
				// 	if (vIndex >= f->vs.size()) {
				// 		logFunc("    Error: Face vertex index " + std::to_string(vIndex) + 
				// 			" is out of bounds (face vs size: " + std::to_string(f->vs.size()) + ")");
				// 		missingVertex = true;
				// 		break;
				// 	}
					
				// 	int vID = f->vs[vIndex];
				// 	if (idVertices(vID) == nullptr) {
				// 		logFunc("    Error: Face " + std::to_string(f->id()) + " references non-existent vertex " + std::to_string(vID));
				// 		missingVertex = true;
				// 		break;
				// 	}
				// }
                
				// if (missingVertex) {
				// 	logFunc("    Skipping face " + std::to_string(f->id()) + " due to missing vertex references");
				// 	continue;
				// }
				
				std::vector<V*> orderVs;
				std::vector<V*> fVs;
				
				// logFunc("    Processing face vertices...");
				// // 检查draw_order数组的边界
				// if (fIndex >= 6/*draw_order[fIndex].size()*/) {
				// 	logFunc("    Error: draw_order index " + std::to_string(fIndex) + 
				// 		" is out of bounds (size: " + std::to_string(4/*draw_order[fIndex].size()*/) + ")");
				// 	continue;
				// }
				
				//revise normal
				for (int vIndex = 0; vIndex < 4; vIndex++)
				{
					if (vIndex >= 4/*draw_order[fIndex].size()*/) {
						// logFunc("      Error: draw_order[" + std::to_string(fIndex) + "] index " + 
						// 	std::to_string(vIndex) + " is out of bounds (size: " + 
						// 	std::to_string(4/*draw_order[fIndex].size()*/) + ")");
						missingVertex = true;
						break;
					}
					
					if (h->draw_order[fIndex][vIndex] >= h->vs.size()) {
						// logFunc("      Error: draw_order value " + std::to_string(h->draw_order[fIndex][vIndex]) + 
						// 	" is out of bounds for hex vertices (size: " + std::to_string(h->vs.size()) + ")");
						missingVertex = true;
						break;
					}
					
					int orderVid = h->vs[h->draw_order[fIndex][vIndex]];
					//logFunc("      Accessing ordered vertex ID: " + std::to_string(orderVid));
					V* orderedV = idVertices(orderVid);
					if (orderedV == nullptr) {
						//logFunc("        Error: Ordered vertex pointer is null!");
						missingVertex = true;
						break;
					} else {
						// logFunc("        Ordered vertex position: (" + std::to_string(orderedV->position()[0]) + ", " + 
						// 	std::to_string(orderedV->position()[1]) + ", " + 
						// 	std::to_string(orderedV->position()[2]) + ")");
						orderVs.push_back(orderedV);
					}
					
					if (vIndex >= f->vs.size()) {
						// logFunc("      Error: Face vertex index " + std::to_string(vIndex) + 
						// 	" is out of bounds (face vs size: " + std::to_string(f->vs.size()) + ")");
						missingVertex = true;
						break;
					}
					
					int faceVid = f->vs[vIndex];
					//logFunc("      Accessing face vertex ID: " + std::to_string(faceVid));
					V* faceV = idVertices(faceVid);
					if (faceV == nullptr) {
						//logFunc("        Error: Face vertex pointer is null!");
						missingVertex = true;
						break;
					} else {
						// logFunc("        Face vertex position: (" + std::to_string(faceV->position()[0]) + ", " + 
						// 	std::to_string(faceV->position()[1]) + ", " + 
						// 	std::to_string(faceV->position()[2]) + ")");
						fVs.push_back(faceV);
					}
				}
				
				if (missingVertex) {
					//logFunc("    Skipping face normal calculation due to missing vertices");
					continue;
				}
				
				logFunc("    Computing normals...");
				CPoint normal;
				CPoint fNormal;
				for (int vIndex = 0; vIndex < 4; vIndex++)
				{
					//logFunc("      Computing cross product for vertex index " + std::to_string(vIndex));
					try {
						CPoint v1 = orderVs[(vIndex+1)%4]->position() - orderVs[vIndex]->position();
						CPoint v2 = orderVs[(vIndex+3)%4]->position() - orderVs[vIndex]->position();
						CPoint crossProduct = v1 ^ v2;
						//logFunc("        Cross product result: (" + std::to_string(crossProduct[0]) + ", " + 
							//std::to_string(crossProduct[1]) + ", " + 
							//std::to_string(crossProduct[2]) + ")");
						normal += crossProduct;
						
						v1 = fVs[(vIndex+1)%4]->position() - fVs[vIndex]->position();
						v2 = fVs[(vIndex+3)%4]->position() - fVs[vIndex]->position();
						crossProduct = v1 ^ v2;
						fNormal += crossProduct;
					} catch (const std::exception& e) {
						//logFunc("        Exception during cross product calculation: " + std::string(e.what()));
						missingVertex = true;
						break;
					} catch (...) {
						//logFunc("        Unknown exception during cross product calculation");
						missingVertex = true;
						break;
					}
				}
				
				if (missingVertex) {
					//logFunc("    Skipping normal normalization due to calculation errors");
					continue;
				}
				
				//logFunc("    Normalizing normals...");
				normal = normal * 0.25;
				double normalNorm = normal.norm();
				if (normalNorm < 1e-10) {
					//logFunc("      Warning: Very small normal norm: " + std::to_string(normalNorm));
					continue;
				}
				normal /= normalNorm;
				//logFunc("      Normalized normal: (" + std::to_string(normal[0]) + ", " + 
					//std::to_string(normal[1]) + ", " + std::to_string(normal[2]) + ")");
							
				fNormal = fNormal * 0.25;
				double fNormalNorm = fNormal.norm();
				if (fNormalNorm < 1e-10) {
					//logFunc("      Warning: Very small face normal norm: " + std::to_string(fNormalNorm));
					continue;
				}
				fNormal /= fNormalNorm;
				//logFunc("      Normalized face normal: (" + std::to_string(fNormal[0]) + ", " + 
					//std::to_string(fNormal[1]) + ", " + std::to_string(fNormal[2]) + ")");
				
				f->normal() = normal;

				double theta = normal * fNormal / (normalNorm * fNormalNorm);
				theta = theta < -1 ? -1 : theta;
				theta = theta > 1 ? 1 : theta;
				double angle = acos(theta)/3.1415926*180;
				//logFunc("      Angle between normals: " + std::to_string(angle) + " degrees");

				if (f->neighbor_hs[0] == h->id())
				{
					if (angle > 90)
					{
						//logFunc("      Flipping face " + std::to_string(f->id()) + " normals (angle = " + std::to_string(angle) + ")");
						int tempId = f->vs[0];
						f->vs[0] = f->vs[3];
						f->vs[3] = tempId;
						tempId = f->vs[1];
						f->vs[1] = f->vs[2];
						f->vs[2] = tempId;

						int tempEId = f->es[0];
						f->es[0] = f->es[3];
						f->es[3] = tempEId;
						tempEId = f->es[1];
						f->es[1] = f->es[2];
						f->es[2] = tempEId;
						//logFunc("      Face vertices and edges reordered successfully");
					}
				}
			}
			
			//logFunc("Exiting revise_hex_face_normal for hex ID: " + std::to_string(h->id()));
		}
		catch (const std::exception& e) {
			logFunc("在修正法线时发生异常: " + std::string(e.what()));
		}
		catch (...) {
			logFunc("在修正法线时发生未知异常");
		}
	}

	template<typename V, typename E, typename F, typename H>
	void topoM<V, E, F, H>::revise_all_face_normal()
	{
		for (std::list<H*>::iterator hite = hs.begin(); hite != hs.end(); hite++)
		{
			H* h = *hite;
			if (!h->boundary()) continue;
			revise_hex_face_normal(h);
		}
	}

	/*output mesh*/
	template<typename V, typename E, typename F, typename H>
	void topoM<V, E, F, H>::output_Quad_mesh(const char* fileName)
	{
		std::fstream _os(fileName, std::fstream::out);
		if (_os.fail())
		{
			fprintf(stderr, "Error is opening file %s\n", fileName);
			return;
		}
		for (std::list<V*>::iterator vite = vs.begin(); vite != vs.end(); vite++)
		{
			V* v = *vite;
			if (v->boundary())
			{
				_os << "Vertex " << v->id() << " " << v->position() << std::endl;
			}
		}
		for (std::list<F*>::iterator fite = fs.begin(); fite != fs.end(); fite++)
		{
			F* f = *fite;
			if (f->boundary())
			{
				_os << "Face " << f->id();
				for (int fvIndex = 0; fvIndex < 4; fvIndex++)
				{
					_os << " " << f->vs[fvIndex];
				}
				_os << std::endl;
			}
		}
		_os.close();
	}

	template<typename V, typename E, typename F, typename H>
	void topoM<V, E, F, H>::_clear()
	{
		
		for (std::list<V*>::iterator vite = vs.begin();vite!=vs.end();vite++)
		{
			V* v = *vite;
			delete v;			
		}
		for (std::list<E*>::iterator eite = es.begin(); eite != es.end(); eite++)
		{
			E* e = *eite;
			delete e;
		}
		for (std::list<F*>::iterator fite = fs.begin(); fite != fs.end(); fite++)
		{
			F* f = *fite;
			delete f;
		}
		for (std::list<H*>::iterator hite = hs.begin(); hite != hs.end(); hite++)
		{
			H* h = *hite;
			delete h;
		}
		
		m_map_vertices.clear();
		m_map_edges.clear();
		m_map_faces.clear();
		m_map_hexs.clear();
	};

	// 打印map信息的辅助函数
	template<typename V, typename E, typename F, typename H>
	void topoM<V,E,F,H>::print_vertices_map() {
		std::string logMsg = "Vertices Map (size=" + std::to_string(m_map_vertices.size()) + "):";
		if (HMeshLib::gettopoMeshLog()) {
			HMeshLib::topoMeshLog(logMsg);
			
			int count = 0;
			for (auto& pair : m_map_vertices) {
				std::string msg = "  [" + std::to_string(count++) + "] key=" + std::to_string(pair.first) + 
					", value=" + (pair.second ? "valid ptr" : "NULL");
				HMeshLib::topoMeshLog(msg);
				
				if (count > 500) {
					HMeshLib::topoMeshLog("  ... and " + std::to_string(m_map_vertices.size() - 20) + " more items");
					break;
				}
			}
		}
	}
	template<typename V, typename E, typename F, typename H>
	void topoM<V,E,F,H>::print_edges_map() {
		std::string logMsg = "Edges Map (size=" + std::to_string(m_map_edges.size()) + "):";
		if (HMeshLib::gettopoMeshLog()) {
			HMeshLib::topoMeshLog(logMsg);
			
			int count = 0;
			for (auto& pair : m_map_edges) {
				std::string msg = "  [" + std::to_string(count++) + "] key=" + std::to_string(pair.first) + 
					", value=" + (pair.second ? "valid ptr" : "NULL");
				HMeshLib::topoMeshLog(msg);
				
				if (false) {
					HMeshLib::topoMeshLog("  ... and " + std::to_string(m_map_edges.size() - 20) + " more items");
					break;
				}
			}
		}
	}
	template<typename V, typename E, typename F, typename H>
	void topoM<V,E,F,H>::print_faces_map() {
		std::string logMsg = "Faces Map (size=" + std::to_string(m_map_faces.size()) + "):";
		if (HMeshLib::gettopoMeshLog()) {
			HMeshLib::topoMeshLog(logMsg);
			
			int count = 0;
			for (auto& pair : m_map_faces) {
				std::string msg = "  [" + std::to_string(count++) + "] key=" + std::to_string(pair.first) + 
					", value=" + (pair.second ? "valid ptr" : "NULL");
				HMeshLib::topoMeshLog(msg);
				
				if (false) {
					HMeshLib::topoMeshLog("  ... and " + std::to_string(m_map_faces.size() - 20) + " more items");
					break;
				}
			}
		}
	}
	template<typename V, typename E, typename F, typename H>
	void topoM<V,E,F,H>::print_hexs_map() {
		std::string logMsg = "Hexs Map (size=" + std::to_string(m_map_hexs.size()) + "):";
		if (HMeshLib::gettopoMeshLog()) {
			HMeshLib::topoMeshLog(logMsg);
			
			int count = 0;
			for (auto& pair : m_map_hexs) {
				std::string msg = "  [" + std::to_string(count++) + "] key=" + std::to_string(pair.first) + 
					", value=" + (pair.second ? "valid ptr" : "NULL");
				HMeshLib::topoMeshLog(msg);
				
				if (count > 500) {
					HMeshLib::topoMeshLog("  ... and " + std::to_string(m_map_hexs.size() - 20) + " more items");
					break;
				}
			}
		}
	}
	template<typename V, typename E, typename F, typename H>
	void topoM<V,E,F,H>::print_vertex(int id) {
		if (HMeshLib::gettopoMeshLog()) {
			auto it = m_map_vertices.find(id);
			if (it != m_map_vertices.end()) {
				V* v = it->second;
				if (v) {
					std::stringstream ss;
					ss << "Vertex " << id << ": ";
					ss << "pos=(" << v->position()[0] << "," << v->position()[1] << "," << v->position()[2] << ") ";
					ss << "neighbor_es.size=" << v->neighbor_es.size() << " ";
					ss << "neighbor_fs.size=" << v->neighbor_fs.size() << " ";
					ss << "neighbor_hs.size=" << v->neighbor_hs.size();
					HMeshLib::topoMeshLog(ss.str());
				} else {
					HMeshLib::topoMeshLog("Vertex " + std::to_string(id) + ": NULL pointer");
				}
			} else {
				HMeshLib::topoMeshLog("Vertex " + std::to_string(id) + ": Not found in map");
			}
		}
	}
	template<typename V, typename E, typename F, typename H>
	void topoM<V,E,F,H>::print_edge(int id) {
		if (HMeshLib::gettopoMeshLog()) {
			auto it = m_map_edges.find(id);
			if (it != m_map_edges.end()) {
				E* e = it->second;
				if (e) {
					std::stringstream ss;
					ss << "Edge " << id << ": ";
					ss << "vs=[";
					for (size_t i = 0; i < e->vs.size(); i++) {
						if (i > 0) ss << ",";
						ss << e->vs[i];
					}
					ss << "] ";
					ss << "neighbor_fs.size=" << e->neighbor_fs.size() << " ";
					ss << "neighbor_hs.size=" << e->neighbor_hs.size();
					HMeshLib::topoMeshLog(ss.str());
				} else {
					HMeshLib::topoMeshLog("Edge " + std::to_string(id) + ": NULL pointer");
				}
			} else {
				HMeshLib::topoMeshLog("Edge " + std::to_string(id) + ": Not found in map");
			}
		}
	}
	template<typename V, typename E, typename F, typename H>
	void topoM<V,E,F,H>::print_face(int id) {
		if (HMeshLib::gettopoMeshLog()) {
			auto it = m_map_faces.find(id);
			if (it != m_map_faces.end()) {
				F* f = it->second;
				if (f) {
					std::stringstream ss;
					ss << "Face " << id << ": ";
					ss << "vs=[";
					for (size_t i = 0; i < f->vs.size(); i++) {
						if (i > 0) ss << ",";
						ss << f->vs[i];
					}
					ss << "] ";
					ss << "es=[";
					for (size_t i = 0; i < f->es.size(); i++) {
						if (i > 0) ss << ",";
						ss << f->es[i];
					}
					ss << "] ";
					ss << "neighbor_hs.size=" << f->neighbor_hs.size();
					HMeshLib::topoMeshLog(ss.str());
				} else {
					HMeshLib::topoMeshLog("Face " + std::to_string(id) + ": NULL pointer");
				}
			} else {
				HMeshLib::topoMeshLog("Face " + std::to_string(id) + ": Not found in map");
			}
		}
	}
	template<typename V, typename E, typename F, typename H>
	void topoM<V,E,F,H>::print_hex(int id) {
		if (HMeshLib::gettopoMeshLog()) {
			auto it = m_map_hexs.find(id);
			if (it != m_map_hexs.end()) {
				H* h = it->second;
				if (h) {
					std::stringstream ss;
					ss << "Hex " << id << ": ";
					ss << "vs=[";
					for (size_t i = 0; i < h->vs.size(); i++) {
						if (i > 0) ss << ",";
						ss << h->vs[i];
					}
					ss << "] ";
					ss << "es=[";
					for (size_t i = 0; i < h->es.size() && i < 10; i++) {
						if (i > 0) ss << ",";
						ss << h->es[i];
					}
					if (h->es.size() > 10) ss << ",...";
					ss << "] ";
					ss << "fs=[";
					for (size_t i = 0; i < h->fs.size(); i++) {
						if (i > 0) ss << ",";
						ss << h->fs[i];
					}
					ss << "]";
					HMeshLib::topoMeshLog(ss.str());
				} else {
					HMeshLib::topoMeshLog("Hex " + std::to_string(id) + ": NULL pointer");
				}
			} else {
				HMeshLib::topoMeshLog("Hex " + std::to_string(id) + ": Not found in map");
			}
		}
	};


}

#endif // !TOPOLOGY_HEX_H
