#ifndef TOPOLOGY_MESH_H
#define TOPOLOGY_MESH_H
#include<vector>
#include<map>
#include<unordered_map>
#include "../Parser/StrUtil.h"
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
			return m_map_vertices[vid];
			/*std::map<int, V*>::iterator  iter = m_map_vertices.find(vid);
			if (iter!= m_map_vertices.end())
			{
				return iter->second;
			}
			else 
				return NULL;*/
		};
		E* idEdges(int eid) { 
			return m_map_edges[eid];
			/*std::map<int, V*>::iterator  iter = m_map_edges.find(vid);
			if (iter != m_map_edges.end())
			{
				return iter->second;
			}
			else
				return NULL;*/
		};
		F* idFaces(int fid) { 
			return m_map_faces[fid];
			/*std::map<int, V*>::iterator  iter = m_map_faces.find(vid);
			if (iter != m_map_faces.end())
			{
				return iter->second;
			}
			else
				return NULL;*/
		};
		H* idHexs(int hid) { 
			return m_map_hexs[hid];
			/*std::map<int, V*>::iterator  iter = m_map_hexs.find(vid);
			if (iter != m_map_hexs.end())
			{
				return iter->second;
			}
			else
				return NULL; */
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
	std::vector<E*> topoM< V, E, F, H>::f_adj_e_in_hex(H* hex, F* face)
	{
		std::vector<F*> results;
		int fIndex = hex->faceIndex(face->id());
		if (fIndex == -1)
		{
			return NULL;
		}
		for (int i = 0; i < 4; i++)
		{
			E* resultE = idEdges(hex->fs[hex->fadje[fIndex][i]]);
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
			
			//�ж��Ƿ���ڹ������
			V* facev = idVertices(fv[0]);
			
			bool construct_new_face = true;
			for (int adjFIndex = 0; adjFIndex < facev->neighbor_fs.size(); adjFIndex++)
			{
				F* adjF = idFaces(facev->neighbor_fs[adjFIndex]);

				if (adjF->is_v_equal(fv))
				{				
					construct_new_face = false;
					h->fs.push_back(adjF->id());
					hexFs[i] = adjF;//���Ѿ����ڵ�face����fs��
					adjF->neighbor_hs.push_back(h->id());
					
					//�����ǰ����ڣ�����Ҫ����Ӧ�߷���hexEs��
					if (i<4)
					{
						int edgeOrder[3][2] = { {0,1},{3,0},{1,2} };
						//���ߴ浽hexEs��
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
				//��ֵ��������Ϣ;
				for (int fvIndex = 0; fvIndex < 4; fvIndex++)
				{
					newf->vs.push_back(fv[fvIndex]);
					V* pV = this->idVertices(fv[fvIndex]);
					pV->neighbor_fs.push_back(newf->id());
				}
				//��ֵ�����Ĺ�ϵ
				newf->neighbor_hs.push_back(h->id());
				hexFs[i] = newf;
				h->fs.push_back(newf->id());
							
				if (i<4)//ǰ�ĸ��棬ÿ���湹��������	
				{
					//������
					//�ҵ��ߵ�������				
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
													
							//������ͱߵĹ�ϵ
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
		

		//�������а����ߵĹ�ϵ		
		int feOrder[6][4] = { {0,2,3,1},{3, 5 ,6 ,4},{6,8,9,7},{9,11,0,10}, {2,11,8,5},{1,4,7,10} };
		for (int hexFIndex = 0; hexFIndex < 6; hexFIndex++)
		{		
			F* hexf = hexFs[hexFIndex];	
			if (hexf->es.size()==4)//˵�������½�������
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
		//�����ߺ�������Ĺ�ϵ
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

			//�ж��Ƿ���ڹ������
			V* facev = idVertices(fv[0]);

			bool construct_new_face = true;
			for (int adjFIndex = 0; adjFIndex < facev->neighbor_fs.size(); adjFIndex++)
			{
				F* adjF = idFaces(facev->neighbor_fs[adjFIndex]);

				if (adjF->is_v_equal(fv))
				{
					construct_new_face = false;
					h->fs.push_back(adjF->id());
					hexFs[i] = adjF;//���Ѿ����ڵ�face����fs��
					adjF->neighbor_hs.push_back(h->id());

					//�����ǰ����ڣ�����Ҫ����Ӧ�߷���hexEs��
					if (i < 4)
					{
						int edgeOrder[3][2] = { {0,1},{3,0},{1,2} };
						//���ߴ浽hexEs��	
						
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

						////���Ѿ����ڵ����ϵı߷���hexEs��

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
				//��ֵ��������Ϣ;
				for (int fvIndex = 0; fvIndex < 4; fvIndex++)
				{
					newf->vs.push_back(fv[fvIndex]);
					V* pV = this->idVertices(fv[fvIndex]);
					pV->neighbor_fs.push_back(newf->id());
				}
				//��ֵ�����Ĺ�ϵ
				newf->neighbor_hs.push_back(h->id());
				hexFs[i] = newf;
				h->fs.push_back(newf->id());

				if (i < 4)//ǰ�ĸ��棬ÿ���湹��������	
				{
					//������
					//�ҵ��ߵ�������				
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

							//������ͱߵĹ�ϵ
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
		//�������а����ߵĹ�ϵ		
		int feOrder[6][4] = { {0,2,3,1},{3, 5 ,6 ,4},{6,8,9,7},{9,11,0,10}, {2,11,8,5},{1,4,7,10} };
		for (int hexFIndex = 0; hexFIndex < 6; hexFIndex++)
		{
			F* hexf = hexFs[hexFIndex];
			if (hexf->es.size() == 4)//˵�������½�������
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
		
		//�����ߺ�������Ĺ�ϵ
		for (int hexEIndex = 0; hexEIndex < 12; hexEIndex++)
		{
			hexEs[hexEIndex]->neighbor_hs.push_back(h->id());
			h->es.push_back(hexEs[hexEIndex]->id());
			
		}
		std::cout << "-----------------------------" << std::endl;

	}

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
		//mark boundary
		mark_boundary();
		singularities = mark_singularity();
		/*compute the face normal*/
		computeNormal();
	}

	template<typename V, typename E, typename F, typename H>
	void topoM<V, E, F, H>::write_Qhex(const char* output)
	{
		for (std::list<V*>::iterator vite = vs.begin(); vite != vs.end(); vite++)
		{
			V* v = *vite;
			v->_to_string();
		}
		for (std::list<E*>::iterator eite = es.begin(); eite != es.end(); eite++)
		{
			E* e = *eite;
			e->_to_string();
		}
		for (std::list<H*>::iterator hite = hs.begin(); hite != hs.end(); hite++)
		{
			H* h = *hite;
			h->_to_string();
		}

		std::fstream _os(output, std::fstream::out);

		if (_os.fail())
		{
			fprintf(stderr, "Error in opening file %s\n", output);
			return;
		}
		//output vertices
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
		}
		//output hex
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
		}

		//output edge
		for (std::list<E*>::iterator eite = es.begin(); eite != es.end(); eite++)
		{
			E* e = *eite;
			if (e->string().size() > 0)
			{
				_os << "Edge " << e->vs[0] << " " << e->vs[1] << " ";
				_os << "{" << e->string() << "}" << std::endl;
			}
		}

		_os.close();
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

	}
	

	/*revised normal*/
	template<typename V, typename E, typename F, typename H>
	void topoM<V, E, F, H>::revise_hex_face_normal(H* h)
	{
		for (int fIndex = 0; fIndex < h->fs.size(); fIndex++)
		{
			F* f = idFaces(h->fs[fIndex]);
			if (f->neighbor_hs[0] != h->id()) continue;
			//if (f->neighbor_hs.size()==1)
			{
				std::vector<V*> orderVs;
				std::vector<V*> fVs;
				//revise normal
				for (int vIndex = 0; vIndex < 4; vIndex++)
				{
					orderVs.push_back(idVertices(h->vs[h->draw_order[fIndex][vIndex]]));
					
					fVs.push_back(idVertices(f->vs[vIndex]));
				}
				
				CPoint normal;
				CPoint fNormal;
				for (int vIndex = 0; vIndex < 4; vIndex++)
				{
					normal+= (orderVs[(vIndex+1)%4]->position() - orderVs[vIndex]->position()) ^ (orderVs[(vIndex + 3) % 4]->position() - orderVs[vIndex]->position());
					fNormal+= (fVs[(vIndex + 1) % 4]->position() - fVs[vIndex]->position()) ^ (fVs[(vIndex + 3) % 4]->position() - fVs[vIndex]->position());
				}			
				normal = normal * 0.25;
				normal /= normal.norm();			
				fNormal = fNormal * 0.25;
				fNormal /= fNormal.norm();
				f->normal() = normal;

				double theta = normal * fNormal / (normal.norm() * fNormal.norm());
				theta = theta < -1 ? -1 : theta;
				theta = theta > 1 ? 1 : theta;
				double angle = acos(theta)/3.1415926*180;

				if (f->neighbor_hs[0]==h->id())
				{
					if (angle > 90)
					{
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

						//f->normal() = normal;
					}
				}
				
				
			}
		}
	}

	template<typename V, typename E, typename F, typename H>
	void topoM<V, E, F, H>::revise_all_face_normal()
	{
		for (std::list<H*>::iterator hite = hs.begin(); hite != hs.end(); hite++)
		{
			H* h = *hite;
			//if (!h->boundary()) continue;
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
	}

}

#endif // !TOPOLOGY_HEX_H