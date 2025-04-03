#pragma once
#include "topoMesh.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>
#include <iomanip> // 添加该头文件以便格式化时间输出

namespace HMeshLib
{
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

		// 添加日志记录函数
		void log(const std::string& message);
		std::ofstream& get_log_stream();

		V* construct_new_vertex();
		E* construct_new_edge(std::vector<V*> vs);
		F* constrcut_new_face(std::vector<V*> vs);
		void vertex_position_policy(V* v1, V* v2, V* target);
		void edge_sharp_policy(E* e1, E* e2, E* target);
		void edge_total_angle_policy_collapse(E* e1, E* e2, E* target);
		void mark_elements_attribute_collapse(std::vector<F*> fs);

		std::vector<E*> get_one_sheet(E* e);
		int get_mesh_sheet_number();
		// 添加重置sheet属性的函数
		void reset_sheet_attributes();
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
		std::ofstream log_file;
		std::string session_id;
	};

	template<typename M>
	sheet_operation<M>::sheet_operation(M* input_mesh)
	{
		mesh = input_mesh;
		
		// 获取当前会话ID
		char* session_id_env = getenv("RL_HMESH_SESSION_ID");
		if (session_id_env) {
			session_id = std::string(session_id_env);
		} else {
			// 如果环境变量不存在，创建默认会话ID
			auto now = std::chrono::system_clock::now();
			std::time_t now_time = std::chrono::system_clock::to_time_t(now);
			std::tm now_tm;
			localtime_s(&now_tm, &now_time);
			
			char timestamp[32];
			std::strftime(timestamp, sizeof(timestamp), "%Y%m%d_%H%M%S", &now_tm);
			session_id = std::string(timestamp);
		}
		
		// 创建日志文件
		std::string log_dir = "F:/RL_HMesh/logs/" + session_id;
		std::string log_filename = log_dir + "/sheet_operation.log";
		log_file.open(log_filename, std::ios::app);
		
		if (log_file.is_open()) {
			log_file << "================ Sheet Operation Log ================" << std::endl;
			log_file << "Session ID: " << session_id << std::endl;
			
			// 修复时间点输出问题
			auto now = std::chrono::system_clock::now();
			std::time_t now_time = std::chrono::system_clock::to_time_t(now);
			std::tm now_tm;
			localtime_s(&now_tm, &now_time);
			log_file << "Log started at: " << std::put_time(&now_tm, "%Y-%m-%d %H:%M:%S") << std::endl;
		} else {
			std::cerr << "Failed to open log file: " << log_filename << std::endl;
		}
	}

	template<typename M>
	void sheet_operation<M>::log(const std::string& message) {
		if (log_file.is_open()) {
			auto now = std::chrono::system_clock::now();
			std::time_t now_time = std::chrono::system_clock::to_time_t(now);
			std::tm now_tm;
			localtime_s(&now_tm, &now_time);
			
			char timestamp[32];
			std::strftime(timestamp, sizeof(timestamp), "%H:%M:%S", &now_tm);
			
			log_file << "[" << timestamp << "] " << message << std::endl;
			log_file.flush();
		}
	}

	template<typename M>
	std::ofstream& sheet_operation<M>::get_log_stream() {
		return log_file;
	}

	template<typename M>
	sheet_operation<M>::sheet_operation()
	{

	}

	template<typename M>
	sheet_operation<M>::~sheet_operation()
	{
		if (log_file.is_open()) {
			log_file << "================ Sheet Operation Log Closed ================" << std::endl;
			log_file.close();
		}
	}

	/*collaspe*/
	template<typename M>
	std::vector<typename sheet_operation<M>::E*> sheet_operation<M>::get_one_sheet(E* e)
	{
		//log("Getting sheet from edge ID: " + std::to_string(e->id()));
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
		//log("Sheet size: " + std::to_string(sheet.size()) + " edges");
		return sheet;
	}

	// 添加重置sheet属性的函数实现
	template<typename M>
	void sheet_operation<M>::reset_sheet_attributes()
	{
		log("Resetting all edges' sheet attributes to 0");
		int reset_count = 0;
		int non_zero_count = 0;

		for (M::MEIterator eite(mesh); !eite.end(); eite++)
		{
			E* e = *eite;
			// 记录非零sheet值的边
			if (e->sheet() != 0) {
				non_zero_count++;
				// 如果有异常大的值，记录下来
				if (e->sheet() > 1000) {
					log("WARNING: Edge ID " + std::to_string(e->id()) + " had unusually large sheet ID: " + std::to_string(e->sheet()));
				}
			}
			// 重置sheet属性为0
			e->sheet() = 0;
			reset_count++;
		}

		log("Reset completed: processed " + std::to_string(reset_count) + " edges, found " + 
			std::to_string(non_zero_count) + " edges with non-zero sheet values");
	}

template<typename M>
int sheet_operation<M>::get_mesh_sheet_number()
{
	log("Starting to get mesh sheet number");
	
	// 第一步：识别所有独立的sheet，但不分配ID
	std::vector<std::vector<E*>> all_sheets;
	log("Step 1: Identifying all independent sheets without assigning IDs");
	
	for (M::MEIterator eite(mesh); !eite.end(); eite++)
	{
		E* be = *eite;
		if (be->sheet() != 0) continue;  // 跳过已处理的边
		
		std::vector<E*> sheet = get_one_sheet(be);
		log("Found sheet starting from edge ID " + std::to_string(be->id()) + " with " + std::to_string(sheet.size()) + " edges");
		
		// 临时标记这个sheet中的所有边，以避免重复处理
		for (int sheetIndex = 0; sheetIndex < sheet.size(); sheetIndex++)
		{
			E* sheete = sheet[sheetIndex];
			sheete->sheet() = -1;  // 使用-1作为临时标记
		}
		
		all_sheets.push_back(sheet);
	}
	
	log("Step 1 completed: Found " + std::to_string(all_sheets.size()) + " independent sheets");
	
	// 第二步：为所有sheet分配唯一ID
	log("Step 2: Assigning unique IDs to all sheets");
	int sheet_id = 0;
	
	// 创建一个映射来跟踪每个sheet ID包含的边数量
	std::map<int, int> sheet_id_edge_count;
	
	for (int i = 0; i < all_sheets.size(); i++)
	{
		std::vector<E*>& sheet = all_sheets[i];
		sheet_id++;
		
		log("Assigning ID " + std::to_string(sheet_id) + " to sheet " + std::to_string(i+1) + 
			" with " + std::to_string(sheet.size()) + " edges");
		
		// 记录这个sheet ID包含的边数量
		sheet_id_edge_count[sheet_id] = sheet.size();
		
		// 为这个sheet中的所有边分配ID
		for (int sheetIndex = 0; sheetIndex < sheet.size(); sheetIndex++)
		{
			E* sheete = sheet[sheetIndex];
			sheete->sheet() = sheet_id;
			
			// 记录边ID和分配给它的sheet ID
			log("Edge ID " + std::to_string(sheete->id()) + " assigned to sheet ID " + std::to_string(sheet_id));
		}
	}
	
	log("Step 2 completed: Assigned " + std::to_string(sheet_id) + " unique sheet IDs");
	
	// 第三步：检查是否有未标记的边
	log("Step 3: Checking for unmarked edges");
	int unmarked_count = 0;
	
	for (M::MEIterator eite(mesh); !eite.end(); eite++)
	{
		E* be = *eite;
		if (be->sheet() <= 0) {  // 检查未标记或临时标记的边
			unmarked_count++;
			log("WARNING: Edge ID " + std::to_string(be->id()) + " was not properly marked with a sheet ID");
			// 为未标记的边分配一个新的sheet ID
			be->sheet() = ++sheet_id;
			sheet_id_edge_count[sheet_id] = 1;
			log("Assigned new sheet ID " + std::to_string(sheet_id) + " to unmarked edge ID " + std::to_string(be->id()));
		}
	}
	
	log("Step 3 completed: Found " + std::to_string(unmarked_count) + " unmarked edges");
	
	// 第四步：收集和分析最终的sheet ID分布
	log("Step 4: Analyzing final sheet ID distribution");
	
	// 去重打印所有sheet ID
	std::set<int> sheet_ids;
	std::map<int, int> sheet_id_count;  // 记录每个sheet ID出现的次数
	
	for (M::MEIterator eite(mesh); !eite.end(); eite++)
	{
		E* be = *eite;
		sheet_ids.insert(be->sheet());
		sheet_id_count[be->sheet()]++;
	}
	
	// 打印所有sheet ID
	std::string idsStr;
	for (int id : sheet_ids) {
		idsStr += std::to_string(id) + " ";
	}
	log("All sheet IDs: " + idsStr);
	
	// 打印每个sheet ID的统计信息
	log("Sheet ID statistics:");
	for (const auto& pair : sheet_id_count) {
		int id = pair.first;
		int count = pair.second;
		int expected = sheet_id_edge_count[id];
		
		log("Sheet ID " + std::to_string(id) + ": assigned to " + 
			std::to_string(expected) + " edges, found in " + 
			std::to_string(count) + " edges");
		
		if (count != expected) {
			log("WARNING: Mismatch for sheet ID " + std::to_string(id) + 
				" - expected " + std::to_string(expected) + 
				" edges but found " + std::to_string(count));
		}
	}
	
	log("Total sheets found: " + std::to_string(sheet_id));
	log("Unique sheet IDs in final mesh: " + std::to_string(sheet_ids.size()));
	
	if (sheet_ids.size() != sheet_id) {
		log("WARNING: Mismatch between total sheets (" + std::to_string(sheet_id) + 
			") and unique sheet IDs (" + std::to_string(sheet_ids.size()) + ")");
		
		// 第五步：重新编号，使sheet ID从1开始连续
		log("Step 5: Renumbering sheet IDs to start from 1");
		
		// 创建一个映射，将当前sheet ID映射到新的连续ID
		std::map<int, int> id_mapping;
		int new_id = 1;
		
		for (int old_id : sheet_ids) {
			id_mapping[old_id] = new_id++;
		}
		
		// 应用新的ID映射
		for (M::MEIterator eite(mesh); !eite.end(); eite++) {
			E* be = *eite;
			be->sheet() = id_mapping[be->sheet()];
		}
		
		// 打印新的sheet ID分布
		std::set<int> new_sheet_ids;
		for (M::MEIterator eite(mesh); !eite.end(); eite++) {
			E* be = *eite;
			new_sheet_ids.insert(be->sheet());
		}
		
		std::string new_idsStr;
		for (int id : new_sheet_ids) {
			new_idsStr += std::to_string(id) + " ";
		}
		log("New sheet IDs after renumbering: " + new_idsStr);
		log("Total unique sheet IDs after renumbering: " + std::to_string(new_sheet_ids.size()));
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
		log("Marking element attributes after collapse, faces count: " + std::to_string(fs.size()));
		for (int fIndex = 0; fIndex < fs.size(); fIndex++)
		{
			F* f = fs[fIndex];
			//mark boundary
			if (f->neighbor_hs.size() == 1)
			{
				f->boundary() = true;
			}
			else
			{
				f->boundary() = false;
			}

			for (int fhIndex = 0; fhIndex < f->neighbor_hs.size(); fhIndex++)
			{
				H* h = mesh->idHexs(f->neighbor_hs[fhIndex]);
				mesh->revise_hex_face_normal(h);
			}
			//compute the normal
		}

		for (int fIndex = 0; fIndex < fs.size(); fIndex++)
		{
			F* f = fs[fIndex];
			//mark face vertex
			for (int fvIndex = 0; fvIndex < f->vs.size(); fvIndex++)
			{
				V* fv = mesh->idVertices(f->vs[fvIndex]);
				if (f->boundary())
				{
					fv->boundary() = true;
				}
				else
				{
					bool is_boundary = false;
					for (int vfIndex = 0; vfIndex < fv->neighbor_fs.size(); vfIndex++)
					{
						F* fvf = mesh->idFaces(fv->neighbor_fs[vfIndex]);
						if (fvf->boundary())
						{
							is_boundary = true;
							break;
						}
					}
					fv->boundary() = is_boundary;
				}
			}

			//mark face edge 
			for (int feIndex = 0; feIndex < f->es.size(); feIndex++)
			{
				E* fe = mesh->idEdges(f->es[feIndex]);
				if (f->boundary())
				{
					fe->boundary() = true;
				}
				else
				{
					bool is_boundary = false;
					for (int efIndex = 0; efIndex < fe->neighbor_fs.size(); efIndex++)
					{
						F* fef = mesh->idFaces(fe->neighbor_fs[efIndex]);
						if (fef->boundary())
						{
							is_boundary = true;
							break;
						}
					}
					fe->boundary() = is_boundary;
				}

				//mark singularity
				if (fe->boundary())
				{
					if (fe->neighbor_hs.size() != 2)
					{
						fe->singularity() = true;
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
					}
					else
					{
						fe->singularity() = false;
					}
				}
				//compute simplification energy
				if (fe->boundary())
				{
					if (fe->sharp())
					{
						int ideal_degree = round(fe->total_angle() / 90.0);
						ideal_degree = ideal_degree == 0 ? 1 : ideal_degree;
						fe->ideal_degree() = ideal_degree;
					}
					else
					{
						fe->ideal_degree() = 2;
					}
				}
				else
				{
					fe->ideal_degree() = 4;
				}
				fe->sim_energy() = (int)fe->ideal_degree() - (int)fe->neighbor_hs.size();
			}
		}
		log("Element attributes marking completed");
	}

	template<typename M>
	void sheet_operation<M>::collapse_one_sheet2(std::vector<E*> sheet)
	{
		log("Starting to collapse sheet with " + std::to_string(sheet.size()) + " edges");
		clock_t start_time = clock();
		/*get all the hexs and the include elements*/
		std::set<H*> delete_hs;//get all hex in sheet
		std::set<F*> delete_fs;//delete faces
		std::set<E*> delete_es;//delete edges
		std::set<V*> delete_vs;//delete vs
		for (int eIndex = 0; eIndex < sheet.size(); eIndex++)
		{
			E* e = sheet[eIndex];
			e->sheetId() = 1;
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
		
		clock_t time1 = clock();
		log("Get all hex in the sheet time: " + std::to_string(double(time1 - start_time)) + "ms");
		log("Total hexes to delete: " + std::to_string(delete_hs.size()));
		log("Total faces to delete: " + std::to_string(delete_fs.size()));
		log("Total edges to delete: " + std::to_string(delete_es.size()));
		log("Total vertices to delete: " + std::to_string(delete_vs.size()));
		
		/*exit if there is self-intersection*/
		bool has_self_intersection = false;
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
				has_self_intersection = true;
				log("ERROR: The current sheet has self-intersection");
				return;
			}
		}
		
		if (!has_self_intersection) {
			log("Self-intersection check passed");
		}

		clock_t time2 = clock();
		log("Self-intersection check time: " + std::to_string(double(time2 - time1)) + "ms");
		
		/*new elements*/
		std::vector<V*> newvs;
		std::vector<E*> newes;
		std::vector<F*> newfs;
		
		/*add new vertices*/
		for (int eIndex = 0; eIndex < sheet.size(); eIndex++)
		{
			E* e = sheet[eIndex];
			V* newv = construct_new_vertex();
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
		
		clock_t time3 = clock();
		log("Added " + std::to_string(newvs.size()) + " new vertices, time: " + std::to_string(double(time3 - time2)) + "ms");
		
		/*create new edges and faces*/
		for (std::set<H*>::iterator hite = delete_hs.begin(); hite != delete_hs.end(); hite++)
		{
			H* h = *hite;
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
					//hf->sheetId() = 1;
					F* parallel_f = mesh->f_parallel_f_in_hex(h, hf);
					parallel_f->sheetId() = 1;
					fs = { hf,parallel_f };
					//pfs.push_back(fs);
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
		
		log("Created " + std::to_string(newes.size()) + " new edges");
		log("Created " + std::to_string(newfs.size()) + " new faces");

		/*delete the elements*/
		clock_t time4 = clock();
		log("Creating new elements time: " + std::to_string(double(time4 - time3)) + "ms");
		
		//delete hex
		int deleted_hexes = 0;
		std::list<H*>::iterator hite = mesh->hs.begin();
		while (hite != mesh->hs.end())
		{
			H* h = *hite;
			if (h == NULL)
			{
				hite = mesh->hs.erase(hite);
				deleted_hexes++;
			}
			else if (h->is_delete())
			{
				mesh->m_map_hexs.erase(h->id());
				hite = mesh->hs.erase(hite);
				deleted_hexes++;
			}
			else
			{
				++hite;
			}
		}
		log("Deleted " + std::to_string(deleted_hexes) + " hexes");

		//delete face
		int deleted_faces = 0;
		std::list<F*>::iterator fite = mesh->fs.begin();
		while (fite != mesh->fs.end())
		{
			F* f = *fite;
			if (f == NULL)
			{
				fite = mesh->fs.erase(fite);
				deleted_faces++;
			}
			else if (f->is_delete())
			{
				mesh->m_map_faces.erase(f->id());
				fite = mesh->fs.erase(fite);
				deleted_faces++;
			}
			else
			{
				++fite;
			}
		}
		log("Deleted " + std::to_string(deleted_faces) + " faces");
		
		//delete edge
		int deleted_edges = 0;
		std::list<E*>::iterator eite = mesh->es.begin();
		while (eite != mesh->es.end())
		{
			E* e = *eite;
			if (e == NULL)
			{
				eite = mesh->es.erase(eite);
				deleted_edges++;
			}
			else if (e->is_delete())
			{
				mesh->m_map_edges.erase(e->id());
				eite = mesh->es.erase(eite);
				deleted_edges++;
			}
			else
			{
				++eite;
			}
		}
		log("Deleted " + std::to_string(deleted_edges) + " edges");
		
		//delete vertex
		int deleted_vertices = 0;
		std::list<V*>::iterator vite = mesh->vs.begin();
		while (vite != mesh->vs.end())
		{
			V* v = *vite;
			if (v == NULL)
			{
				vite = mesh->vs.erase(vite);
				deleted_vertices++;
			}
			else if (v->is_delete())
			{
				mesh->m_map_vertices.erase(v->id());
				vite = mesh->vs.erase(vite);
				deleted_vertices++;
			}
			else
			{
				++vite;
			}
		}
		log("Deleted " + std::to_string(deleted_vertices) + " vertices");

		clock_t time5 = clock();
		log("Deleting elements time: " + std::to_string(double(time5 - time4)) + "ms");

		mark_elements_attribute_collapse(newfs);
		mesh->singularities = mesh->mark_singularity();
		
		clock_t time6 = clock();
		log("Mark elements attribute time: " + std::to_string(double(time6 - time5)) + "ms");
		log("Total collapse time: " + std::to_string(double(time6 - start_time)) + "ms");
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

		double total_angle = 0;
		for (M::EHIterator ehite(mesh, e); !ehite.end(); ehite++)
		{
			H* eh = *ehite;
			std::vector<F*> ehfs = mesh->e_adj_f_in_hex(eh, e);
			std::vector<CPoint> fn;
			for (int fIndex = 0; fIndex < ehfs.size(); fIndex++)
			{
				F* f = ehfs[fIndex];

				CPoint n = f->normal();
				if (f->neighbor_hs[0] == eh->id())
				{
					fn.push_back(n);

				}
				else
				{
					fn.push_back(-n);
				}
			}
			double temp_angle = 180 - vector_angle(fn[0], fn[1]);
			//int ehIndex = eh->edgeIndex(e->id());
			eh->edge_angle(e->id()) = temp_angle;
			total_angle += temp_angle;
		}
		if (e->total_angle() != 0)
		{
			return;
		}
		if (!e->boundary())
		{
			e->total_angle() = 360.0;
		}
		else
		{
			e->total_angle() = total_angle;

			//e->total_angle() = 180;
			//e->total_angle() = total_angle;
		}
		e->ave_angle() = e->total_angle() / e->neighbor_hs.size();
	}

	template<typename M>
	void sheet_operation<M>::edge_ideal_degree(E* e)
	{
		if (e->ideal_degree() != 0)
		{
			return;
		}
		/*int ideal_degree = round(e->total_angle() / 90.0);
		ideal_degree = ideal_degree == 0 ? 1 : ideal_degree;
		e->ideal_degree() = ideal_degree;*/

		if (e->boundary())
		{
			if (e->sharp())
			{

				int ideal_degree = round(e->total_angle() / 90.0);
				ideal_degree = ideal_degree == 0 ? 1 : ideal_degree;
				e->ideal_degree() = ideal_degree;
				if (mesh->feature_ideal_degree.find(e->sharp()) == mesh->feature_ideal_degree.end())
				{
					std::vector<E*> sharps;
					sharps.push_back(e);
					mesh->feature_ideal_degree.insert(std::pair<int, int>(e->sharp(), ideal_degree));
				}

			}
			else
			{
				e->ideal_degree() = 2;
			}
		}
		else
		{
			e->ideal_degree() = 4;
		}
	}

	template<typename M>
	void sheet_operation<M>::compute_edge_energy()
	{
		log("Computing edge energy for all edges");
		int edge_count = 0;
		mesh->computeNormal();
		for (M::MEIterator eite(mesh); !eite.end(); eite++)
		{

			E* e = *eite;
			edge_count++;
			edge_angle(e);
			edge_ideal_degree(e);
			e->sim_energy() = (int)e->ideal_degree() - (int)e->neighbor_hs.size();
		}
		log("Edge energy computation completed for " + std::to_string(edge_count) + " edges");
	}

	template<typename M>
	std::vector<std::vector<typename sheet_operation<M>::E*>> sheet_operation<M>::get_sheet_parallel_edges(std::vector<E*>sheet)
	{
		log("Getting parallel edges for sheet with " + std::to_string(sheet.size()) + " edges");
		std::vector<H*> hs;
		std::vector<std::vector<E*>> parallel_es;
		for (int eIndex = 0; eIndex < sheet.size(); eIndex++)
		{
			E* e = sheet[eIndex];
			for (int ehIndex = 0; ehIndex < e->neighbor_hs.size(); ehIndex++)
			{
				//get the neighbor hex
				H* eh = mesh->idHexs(e->neighbor_hs[ehIndex]);
				if (eh->mark()) continue;
				eh->mark() = true;
				hs.push_back(eh);

				//get the parallel edges
				std::vector<F*> nfs = mesh->e_adj_f_in_hex(eh, e);
				F* nf = nfs[0];

				V* ev1 = mesh->idVertices(e->vs[0]);
				V* ev2 = mesh->idVertices(e->vs[1]);

				E* ne1 = mesh->flip_e(nf, e, ev1);
				E* ne2 = mesh->flip_e(nf, e, ev2);
				F* nf1 = mesh->flip_f(eh, nf, ne1);
				F* nf2 = mesh->flip_f(eh, nf, ne2);
				//std::cout << "eh id: " << eh->id() << std::endl;

				for (int peIndex = 0; peIndex < 4; peIndex++)
				{
					ne1 = mesh->flip_e(nf1, ne1, ev1);
					ev1 = mesh->flip_v(ne1, ev1);
					ne2 = mesh->flip_e(nf2, ne2, ev2);
					ev2 = mesh->flip_v(ne2, ev2);

					if (ne1->mark()) continue;

					ne1->mark() = true;
					ne2->mark() = true;
					std::vector<E*> pair_es = { ne1,ne2 };
					parallel_es.push_back(pair_es);

				}
			}
		}

		//remark the mark
		for (int hIndex = 0; hIndex < hs.size(); hIndex++)
		{
			H* h = hs[hIndex];
			h->mark() = false;
			for (int heIndex = 0; heIndex < h->es.size(); heIndex++)
			{
				E* he = mesh->idEdges(h->es[heIndex]);
				he->mark() = false;
			}
		}
		log("Found " + std::to_string(parallel_es.size()) + " parallel edge pairs");
		return parallel_es;
	}

	template<typename M>
	double sheet_operation<M>::predict_sheet_collapse_energy(std::vector<E*> sheet)
	{
		log("Predicting collapse energy for sheet with " + std::to_string(sheet.size()) + " edges");
		/*energy = ideal-real，energyԽСԽ��*/
		double original_total_energy = 0;
		double predict_total_energy = 0;
		double cannot_collapse = -99999;//99999 means the the sheet can't be collapsed

		//can't collapse two corners
		for (int eIndex = 0; eIndex < sheet.size(); eIndex++)
		{
			E* e = sheet[eIndex];
			V* ev1 = mesh->idVertices(e->vs[0]);
			V* ev2 = mesh->idVertices(e->vs[1]);

			if (ev1->corner() && ev2->corner())
			{
				log("Cannot collapse: sheet has two corners");
				return cannot_collapse;
			}
			if (ev1->feature_vertex() < 0 && ev2->feature_vertex() < 0)
			{
				log("Cannot collapse: sheet has two feature vertices < 0");
				return cannot_collapse;
			}
			if (ev1->feature_vertex() < 0 && ev2->feature_vertex() > 0)
			{
				int sharp_id = ev2->feature_vertex();
				int corner_id = ev1->feature_vertex();
				std::vector<int>connect_corners_id = mesh->feature_edge_corner[sharp_id];
				std::vector<int>::iterator ite = std::find(connect_corners_id.begin(), connect_corners_id.end(), corner_id);
				if (ite == connect_corners_id.end())
				{
					log("Cannot collapse: sheet has a current feature vertex is going to connected to other corner");
					return cannot_collapse;
				}
			}
			if (ev1->feature_vertex() > 0 && ev2->feature_vertex() < 0)
			{
				int sharp_id = ev1->feature_vertex();
				int corner_id = ev2->feature_vertex();
				std::vector<int>connect_corners_id = mesh->feature_edge_corner[sharp_id];
				std::vector<int>::iterator ite = std::find(connect_corners_id.begin(), connect_corners_id.end(), corner_id);
				if (ite == connect_corners_id.end())
				{
					log("Cannot collapse: sheet has a current feature vertex is going to connected to other corner");
					return cannot_collapse;
				}
			}

			if (ev1->feature_vertex() > 0 && ev2->feature_vertex() > 0)
			{
				if (ev1->feature_vertex() != ev2->feature_vertex())
				{
					log("Cannot collapse: sheet has two different feature vertex will be merge");
					return cannot_collapse;
				}
			}
		}

		//get all edges in the sheet
		std::vector<std::vector<E*>>parallel_es = get_sheet_parallel_edges(sheet);

		//evaluate the collapse or not
		int temp_less_num = 0;
		int temp_equal_num = 0;
		int temp_large_num = 0;
		int temp_boundary_num = 0;//eatimate the bounary collapse parameter

		for (int peIndex = 0; peIndex < parallel_es.size(); peIndex++)
		{
			std::vector<E*> one_pes = parallel_es[peIndex];
			E* e1 = one_pes[0];
			E* e2 = one_pes[1];
			
			if (e1->sharp() && e2->sharp())
			{
				log("Cannot collapse: sheet has two parallel sharp edges");
				return cannot_collapse;
			}
			double predict_angle = 0;
			double predict_degree = 0;

			//compute predict angle
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

			//compute predict degree
			if (e1->boundary() && e2->boundary())
			{
				predict_degree = e1->neighbor_hs.size() + e2->neighbor_hs.size() - 2;
			}
			else
			{
				predict_degree = e1->neighbor_hs.size() + e2->neighbor_hs.size() - 4;
			}
			//The degree of too small edge is prohibited
			if (e1->boundary() || e2->boundary())
			{
				if (predict_degree < 1)
				{
					log("Cannot collapse: boundary predict edge degree < 1");
					return cannot_collapse;
				}
			}
			else
			{
				if (predict_degree < 2)
				{
					log("Cannot collapse: inner predict edge degree < 2");
					return cannot_collapse;
				}
			}


			//compute the energy
			int temp_predict_energy = abs(round(predict_angle / 90.0) - predict_degree);
			predict_total_energy += temp_predict_energy;

			int temp_original_energy1 = 0;
			int temp_original_energy2 = 0;
			original_total_energy += abs(e1->sim_energy());
			original_total_energy += abs(e2->sim_energy());
			temp_original_energy1 += abs(e1->sim_energy());
			temp_original_energy2 += abs(e2->sim_energy());

			//boundary energy		
			int boundary_energy1 = 0;
			int boundary_energy2 = 0;
			if (e1->boundary())
			{
				boundary_energy1 = abs(e1->sim_energy()) - temp_predict_energy;
			}
			if (e2->boundary())
			{
				boundary_energy2 = abs(e2->sim_energy()) - temp_predict_energy;
			}
			//estimate the boundary collapse energy
			if ((e1->boundary() || e2->boundary()) && !(e1->boundary() && e2->boundary()))
			{
				if (temp_boundary_num <= 0)
				{
					if (boundary_energy1 > 0 || boundary_energy2 > 0)
					{
						temp_boundary_num = 1;
					}
					else if (boundary_energy1 == 0 && boundary_energy2 == 0)
					{
						temp_boundary_num = temp_boundary_num;
					}
					else
					{
						temp_boundary_num = -1;
					}
				}
			}

			//estimate the collapse energy
			if (temp_original_energy1 > temp_predict_energy || temp_original_energy2 > temp_predict_energy)
			{
				temp_large_num++;
			}
			else if (temp_original_energy1 == temp_predict_energy || temp_original_energy2 == temp_predict_energy)
			{
				temp_equal_num++;
			}
			else if (temp_original_energy1 < temp_predict_energy && temp_original_energy2 < temp_predict_energy)
			{
				temp_less_num++;
			}
		}
		// 最终决策日志
		if (temp_boundary_num < 0)
		{
			log("Cannot collapse: sheet creates boundary edges with larger energy");
			return -999;
		}
		else if (temp_large_num > 0)
		{
			double energy_diff = original_total_energy - predict_total_energy;
			log("Can collapse: positive energy gain " + std::to_string(energy_diff));
			return energy_diff;
		}
		else if (temp_less_num == 0)
		{
			double energy_diff = original_total_energy - predict_total_energy;
			log("Can collapse: no negative energy " + std::to_string(energy_diff));
			return energy_diff;
		}
		else
		{
			log("Cannot collapse: sheet doesn't have positive energy");
			return -999;
		}
	}
}
