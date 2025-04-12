#pragma once
#ifndef MY_BASE_COMPLEX_MESH_H
#define MY_BASE_COMPLEX_MESH_H
#include<core/Topology/TopoIterator.h>
#include<core/Topology/TopoOperation.h>
#include<unordered_set>
#include <core/Parser/parser.h>
namespace HMeshLib
{
	class TV;
	class TE;
	class TF;
	class TH;

	class singularity;

	class singularity
	{
	public:
		singularity() {
			m_id = 0;
			m_loop = false;
			m_degree = 0;
			m_boundary = false;
		};
		~singularity() {};
		std::vector<TE*> es;
		std::vector<TV*> vs;
		int& id() { return m_id; };
		bool& loop() { return m_loop; };
		int& degree() { return m_degree; };
		bool& boundary() { return m_boundary; };
	private:
		int m_id;
		bool m_loop;
		int m_degree;
		bool m_boundary;
	};

	class TV :public topoV
	{
	public:
		TV() {
			m_outside = false;
			m_base_complex_node = false;
			m_saved = 0;
			m_coupled = false;
			m_this_sheet = false;
			m_corner = false;
			m_feature_curve = 0;
			m_boundary_id = 0;
		};
		~TV() {};

		bool& outside() { return m_outside; };
		bool& base_complex_node() { return m_base_complex_node; };//锟�?�凤拷锟斤拷base_complex 锟侥节�?�拷
		int& saved() { return m_saved; };
		bool& coupled() { return m_coupled; };//锟斤拷sheet锟斤拷时锟斤拷锟角凤拷锟�?碉拷couple
		bool& this_sheet() { return m_this_sheet; };//sheet锟津化的憋拷牵锟斤拷锟绞憋拷锟斤拷锟轿猼rue
		bool& corner() { return m_corner; };//锟�?�凤拷锟角�?��?�拷
		int& feature_curve() { return m_feature_curve; };//锟斤拷锟斤拷锟斤拷锟斤拷锟�??碉拷id
		int& boundary_id() { return m_boundary_id; };//通锟斤拷锟斤拷锟斤拷锟�??伙拷锟�?��?�拷锟斤拷id
		int& degree() { return m_degree; };
		int& mark() { return m_mark; };
		/*!	Vertex spherical harmonic map image coordinates */
		CPoint& u() { return m_u; };
		/*!	Vertex Laplacian*/
		CPoint& L() { return m_L; };
		int& feature_vertex() { return m_feature_vertex; };
		bool& is_delete() { return m_is_delete; };
		bool& sheet_inflate() { return m_sheet_inflate; };
		int& newv() { return m_newv; };
		int& sheet() { return m_sheet; };
		std::vector<TF*> stc_fs;
		bool& not_output() { return m_not_output; };
		CPoint& pre_position() {
			return m_pre_position;
		};
		void _from_string()
		{
			CParser parser(m_string);
			for (std::list<CToken*>::iterator iter = parser.tokens().begin(); iter != parser.tokens().end(); ++iter)
			{
				CToken* token = *iter;
				if (token->m_key == "feature_vertex")
				{
					std::string feature_string = strutil::trim(token->m_value, "()");//锟揭碉拷锟斤拷GP锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷锟�	
					//m_feature_vertex = (feature_string =="1");
					sscanf_s(feature_string.c_str(), "%d ", &m_feature_vertex);
				}
			}
		}
		void _to_string()
		{
			CParser parser(m_string);
			parser._removeToken("feature_vertex");

			parser._toString(m_string);

			std::string line;
			std::stringstream iss(line);
			iss << "feature_vertex=(" << m_feature_vertex << ")";
			/*iss << " normal=(" << m_normal<<")";*/
			if (m_string.length() > 0)
			{
				m_string += " ";
			}
			m_string += iss.str();
		}
	protected:

		bool m_outside;
		bool m_base_complex_node;
		int m_saved;
		bool m_coupled;
		bool m_this_sheet;
		bool m_corner;
		int m_feature_curve;
		int m_boundary_id;
		int m_degree;
		int m_mark;
		CPoint m_u;
		CPoint m_L;
		int m_feature_vertex;
		bool m_is_delete;
		bool m_sheet_inflate;
		int m_newv;
		int m_sheet;
		bool m_not_output;
		CPoint m_pre_position;
	};

	class TE :public topoE
	{
	public:
		TE()
		{
			m_base_complex_e_id = 0;
			m_singularity_id = 0;
			m_saved = false;
			m_this_sheet = false;
			m_sharp = 0;
			m_stc = 0;
		};
		~TE() {};

		bool& outside() { return m_outside; };
		int& base_complex_e_id() { return m_base_complex_e_id; };
		int& singularity_id() { return m_singularity_id; };
		bool& saved() { return m_saved; };
		bool& this_sheet() { return m_this_sheet; };
		int& sharp() { return m_sharp; };//锟斤拷锟斤拷锟斤�?
		int& stc() { return m_stc; };
		int& mark() { return m_mark; };
		/*! Edge weight */
		double& weight() { return m_weight; };
		bool& is_delete() { return m_is_delete; };
		int& newe() { return m_new_e; };
		bool& sheet_inflate() { return m_sheet_inflate; };
		double& sim_energy() { return m_sim_energy; };
		int& ideal_degree() { return m_ideal_degree; };
		void setNewv(TV* inputV) { m_newv = inputV; };
		TV* newv() { return m_newv; };
		double& ave_angle() { return m_ave_angle; };
		double& total_angle() { return m_total_angle; };
		int& sheet() { return m_sheet; };

		std::map<int, TV*> sheet_newv;//int: sheet id. TV new vertex of sheet
		void _from_string()
		{
			CParser parser(m_string);

			for (std::list<CToken*>::iterator iter = parser.tokens().begin(); iter != parser.tokens().end(); ++iter)
			{
				CToken* token = *iter;
				if (token->m_key == "sharp")
				{
					std::string value = strutil::trim(token->m_value, "()");
					std::stringstream sStream;
					sStream << value;
					sStream >> m_sharp;
				}
				if (token->m_key == "sheet")
				{
					m_sheet = 1;
				}
			}
		}
		void _to_string()
		{
			CParser parser(m_string);
			parser._removeToken("zero_length");
			parser._removeToken("sharp");
			parser._toString(m_string);

			if (m_sharp)
			{
				m_string = m_string + "sharp=(" + (std::to_string(m_sharp)) + ")";
			}
		}
	private:
		bool m_outside;
		int m_base_complex_e_id;
		int m_singularity_id;
		bool m_saved;
		bool m_this_sheet;
		int m_sharp;
		int m_stc;
		int m_mark;
		double m_weight;
		bool m_is_delete;
		int m_new_e;
		bool m_sheet_inflate;
		double m_sim_energy;
		int m_ideal_degree;
		TV* m_newv;
		double m_ave_angle;
		double m_total_angle;
		int m_sheet;

	};

	class TF :public topoF
	{
	public:
		TF()
		{
			m_base_complex_f_id = 0;
			m_base_complex_segmentation_f_id = 0;
			m_saved = false;
			m_this_sheet = false;
			m_boundary_id = 0;
			m_traced = false;
			m_delete = false;
		};
		~TF() {};
		bool& inverse() { return m_inverse; };
		int& base_complex_f_id() { return m_base_complex_f_id; };
		int& base_complex_segmentation_f_id() { return m_base_complex_segmentation_f_id; };
		bool& saved() { return m_saved; };
		bool& this_sheet() { return m_this_sheet; };
		int& boundary_id() { return m_boundary_id; };
		bool& traced() { return m_traced; };
		//TV* newv;
		bool& is_delete() { return m_delete; };
		int& mark() { return m_mark; };
		bool& sheet_inflate() { return m_sheet_inflate; };
		int& newf() { return m_newf; };
		void setNewv(TV* inputV) { m_newv = inputV; };
		TV* newv() { return m_newv; };
		void _from_string()
		{
			CParser parser(m_string);
			//std::cout << "m string: " << m_string << std::endl;
			for (std::list<CToken*>::iterator iter = parser.tokens().begin(); iter != parser.tokens().end(); ++iter)
			{
				CToken* token = *iter;
				if (token->m_key == "inflate")
				{
					m_sheet_inflate = true;
				}
			}
		};
	protected:
		bool m_inverse;//draw face param
		int m_base_complex_f_id;
		int m_base_complex_segmentation_f_id;
		bool m_saved;
		bool m_this_sheet;
		int m_boundary_id;
		bool m_traced;
		bool m_delete;
		int m_mark;
		bool m_sheet_inflate;
		int m_newf;
		TV* m_newv;
	};

	class TH :public topoH
	{
	public:
		TH()
		{
			m_base_complex_h_id = 0;
			m_this_sheet = false;
			m_stc = 0;
			m_delete = false;
		};
		~TH() {};


		bool& outside() { return m_outside; };
		int& base_complex_h_id() { return m_base_complex_h_id; };
		bool& saved() { return m_saved; };
		bool& this_sheet() { return m_this_sheet; };
		int& stc() { return m_stc; };
		bool& is_delete() { return m_delete; };
		int& mark() { return m_mark; };
		bool& sheet_inflate() { return m_sheet_inflate; };
		int& part() { return m_part; };
		int& sheet() { return m_sheet; };
		bool& not_output() { return m_not_output; };
		std::vector<TF*> stc_fs;
		TE* newe;

		void _from_string()
		{

		};
	protected:
		bool m_outside;
		int m_base_complex_h_id;
		bool m_saved;
		bool m_this_sheet;
		int m_stc;
		bool m_delete;
		int m_mark;
		bool m_sheet_inflate;
		int m_part;
		int m_sheet;
		bool m_not_output;
	};


	template<typename V, typename E, typename F, typename H>
	class TopoMesh :public topoM< V, E, F, H>
	{

	public:
		typedef V V;
		typedef E E;
		typedef F F;
		typedef H H;

		typedef TopoMeshVertexIterator< V, E, F, H>  MVIterator;
		typedef TopoMeshEdgeIterator< V, E, F, H>  MEIterator;
		typedef TopoMeshFaceIterator< V, E, F, H>  MFIterator;
		typedef TopoMeshHexIterator< V, E, F, H>  MHIterator;
		typedef TopoVertexVertexIterator< V, E, F, H> VVIterator;
		typedef TopoVertexEdgeIterator< V, E, F, H> VEIterator;
		typedef TopoVertexFaceIterator< V, E, F, H> VFIterator;
		typedef TopoVertexHexIterator< V, E, F, H> VHIterator;
		typedef TopoEdgeFaceIterator< V, E, F, H> EFIterator;
		typedef TopoEdgeHexIterator< V, E, F, H> EHIterator;
		typedef TopoFaceHexIterator< V, E, F, H> FHIterator;
		typedef TopoHexVertexIterator< V, E, F, H> HVIterator;
		typedef TopoHexEdgeIterator< V, E, F, H> HEIterator;
		typedef TopoHexFaceIterator< V, E, F, H> HFIterator;
		typedef TopoVertexVertexInFaceIterator< V, E, F, H> VVInFaceIterator;
		typedef TopoVertexVertexInHexIterator< V, E, F, H> VVInHexIterator;
		typedef TopoVertexEdgeInHexIterator<V, E, F, H> VEInHexIterator;
		typedef TopoVertexFaceInHexIterator<V, E, F, H> VFInHexIterator;
		typedef TopoEdgeFaceInHexIterator<V, E, F, H> EFInHexIterator;
		typedef TopoEdgeParalelEdgeIterator< V, E, F, H> EPEInHexItarator;
		typedef TopoEdgeAdjEdgeInHexIterator< V, E, F, H> EAEInHexIterator;
		typedef TopoFaceAdjFaceInHexIterator< V, E, F, H> FAFInHexIterator;
	public:
		void _cut(CPlane& p);
		~TopoMesh() 
		{
			// 锟斤拷锟斤拷 m_Faces_Above 锟叫碉拷指锟斤拷
			for (auto face : m_Faces_Above) {
				delete face;
			}
			m_Faces_Above.clear();

			// 锟斤拷锟斤拷 m_Faces_Below 锟叫碉拷指锟斤拷
			for (auto face : m_Faces_Below) {
				delete face;
			}
			m_Faces_Below.clear();

			// 锟斤拷锟斤拷 m_Edges_Above 锟叫碉拷指锟斤拷
			for (auto edge : m_Edges_Above) {
				delete edge;
			}
			m_Edges_Above.clear();

			// 锟斤拷锟斤拷 m_Edges_Below 锟叫碉拷指锟斤拷
			for (auto edge : m_Edges_Below) {
				delete edge;
			}
			m_Edges_Below.clear();

			// 锟斤拷锟斤拷 classify_features 锟叫碉拷指锟斤拷
			for (auto& edge_vec : classify_features) {
				for (auto edge : edge_vec) {
					delete edge;
				}
				edge_vec.clear();
			}
			classify_features.clear();

		};
		/*value*/
		std::vector<F*> m_Faces_Above;
		std::vector<F*> m_Faces_Below;
		std::vector<E*> m_Edges_Above;
		std::vector<E*> m_Edges_Below;

		/*features*/
		std::vector<std::vector<E*>> classify_features;
		std::unordered_map<int, int> feature_ideal_degree;
		std::unordered_map<int, std::vector<int>> feature_edge_corner;
	};
	template<typename V, typename E, typename F, typename H>
	void TopoMesh<V, E, F, H>::_cut(CPlane& p)
	{
		m_Faces_Above.clear();
		m_Faces_Below.clear();
		m_Edges_Below.clear();
		m_Edges_Above.clear();

		for (std::list<V*>::iterator vite = vs.begin(); vite != vs.end(); vite++)
		{
			V* pV = *vite;
			CPoint pos = pV->position();
			pV->outside() = (p.side(pos) >= 0);
		}



		for (std::list<H*>::iterator hite = hs.begin(); hite != hs.end(); hite++)
		{
			H* pH = *hite;
			V* pV[8];
			for (int i = 0; i < 8; i++)
			{
				pV[i] = idVertices(pH->vs[i]);
			}
			pH->outside() = true;
			for (int i = 0; i < 8; i++)
			{
				if (pV[i]->outside())
				{
					pH->outside() = false;
					break;
				}
			}
		}
		for (std::list<E*>::iterator eite = es.begin(); eite != es.end(); eite++)
		{
			E* pE = *eite;
			/*if (pE->boundary())
			{
				m_Edges_Above.push_back(pE);
			}*/
			int outside = 0;
			int inside = 0;
			for (int i = 0; i < pE->neighbor_hs.size(); i++)
			{
				H* h = idHexs(pE->neighbor_hs[i]);
				if (h->outside())
				{
					outside++;
				}
				else
				{
					inside++;
				}
			}

			if (outside > 0 && inside > 0)
			{
				m_Edges_Above.push_back(pE);
			}
			if (outside > 0 && pE->boundary())
			{
				m_Edges_Above.push_back(pE);
			}
			if (outside == 0 && pE->boundary())
			{
				m_Edges_Below.push_back(pE);
			}
		}


		for (MFIterator fite(this); !fite.end(); fite++)
		{
			F* pF = *fite;
			if (pF->neighbor_hs.size() == 1)
			{
				H* tempH = idHexs(pF->neighbor_hs[0]);
				if (!tempH->outside())
				{
					m_Faces_Below.push_back(pF);
				}
				else
				{
					m_Faces_Above.push_back(pF);
				}
			}
			else
			{
				H* pH[2];
				for (int i = 0; i < pF->neighbor_hs.size(); i++)
				{
					pH[i] = idHexs(pF->neighbor_hs[i]);
				}
				if (pH[0]->outside() && !pH[1]->outside())
				{
					m_Faces_Above.push_back(pF);
					//m_Faces_Below.push_back(pF);
					pF->inverse() = false;
				}
				else if (!pH[0]->outside() && pH[1]->outside())
				{
					//m_Faces_Below.push_back(pF);
					m_Faces_Above.push_back(pF);
					pF->inverse() = true;
				}

			}
		}
	}

	typedef TopoMesh<TV, TE, TF, TH> TMesh;

	template<typename V, typename E, typename F, typename H, typename M>
	class TopoSheet :public sheet<V, E, F, H, M>
	{

	};
	typedef TopoSheet<TV, TE, TF, TH, TMesh> TSheet;

	template<typename V, typename E, typename F, typename H, typename M>
	class TopoColumn :public column<V, E, F, H, M>
	{

	};
	typedef TopoColumn<TV, TE, TF, TH, TMesh> TColumn;

	template<typename V, typename E, typename F, typename H, typename M>
	class TopoOperator :public topoOperation<V, E, F, H, M>
	{};

	typedef TopoOperator<TV, TE, TF, TH, TMesh> TOperator;

}


#endif // !MY_BASE_COMPLEX_MESH_H
