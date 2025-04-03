#ifndef TOPOLOGICAL_OPERATION_H
#define TOPOLOGICAL_OPERATION_H
#include<time.h>
#include<unordered_set>
namespace HMeshLib
{
	template<typename V, typename E, typename F, typename H, typename M>
	class sheet
	{
	public:
		sheet() { m_id = -1; };
		~sheet() {};

		int& id() { return m_id; };
		std::vector<E*> sheets;
		std::vector<H*> hs;
		std::vector<std::vector<F*>> couple_faces;

	protected:
		int m_id;
	};

	template<typename V, typename E, typename F, typename H, typename M>
	class column
	{
	public:
		column() { m_id = 0; };
		~column() {};
		std::vector<F*>& faces() { return m_faces; };
		std::vector<H*>& hexs() { return m_hexs; };
		std::vector<E*>& es1() { return m_es1; };
		std::vector<E*>& es2() { return m_es2; };
		int& id() { return m_id; };
	private:
		std::vector<F*> m_faces;
		std::vector<H*> m_hexs;
		std::vector<E*> m_es1;
		std::vector<E*> m_es2;
	protected:
		int m_id;
	};

	template<typename V, typename E, typename F, typename H, typename M>
	class section
	{
	public:
		section() { m_sectionId = -1; };
		~section() {};
		std::vector<F*>& faces() { return m_faces; };
		std::vector<E*>& fixed_edges() { return m_fixed_edges; };
		std::vector<std::unordered_set<V*>>& layers_vertieces() { return m_layers_vertices; };
		int& sectionId() { return m_sectionId; };
		void init_sectionId();
	private:
		std::vector<F*> m_faces;
		std::vector<E*> m_fixed_edges;
		int m_sectionId;
		std::vector<std::unordered_set<V*>> m_layers_vertices;

	};
	
	template<typename V, typename E, typename F, typename H, typename M>
	void section<V,E,F,H,M>::init_sectionId()
	{
		for (int i = 0; i < m_faces.size(); i++)
		{
			F* f = m_faces[i];
			f->sectionId() = m_sectionId;
		}

	}

	template<typename V, typename E, typename F, typename H, typename M>
	class both_end_section 
	{
	public:
		both_end_section() { m_sectionId = -1; };
		~both_end_section() {};
		std::vector<F*>& faces() { return m_faces; };
		std::vector<E*>& fixed_edges1() { return m_fixed_edges1; };
		std::vector<E*>& fixed_edges2() { return m_fixed_edges2; };
		std::vector<std::unordered_set<V*>>& layers_vertieces() { return m_layers_vertices; };
		int& sectionId() { return m_sectionId; };
		void init_sectionId();
	private:
		std::vector<F*> m_faces;
		std::vector<E*> m_fixed_edges1;
		std::vector<E*> m_fixed_edges2;
		int m_sectionId;
		std::vector<std::unordered_set<V*>> m_layers_vertices;
	};


	template<typename V, typename E, typename F, typename H, typename M>
	class topoOperation
	{
		typedef section<V, E, F, H, M> section;
		typedef both_end_section<V, E, F, H, M> both_end_section;
		typedef sheet<V, E, F, H, M> sheet;
		typedef column<V, E, F, H, M> column;
	public:
		topoOperation() {};
		~topoOperation() {};
		void singularity_insert(M*mesh,section* s, int insert_layer_num);
		void both_end_singularity_insert(M* mesh, both_end_section* s, int insert_layer_num);
		//void singularity_insert(M* mesh, section* s);
		int num_of_layers(M* mesh, section* s);
		int both_end_num_of_layers(M* mesh, both_end_section* s);

		/*sheet opeation*/
		void get_sheet(M* mesh, E* e,sheet* sh);
		void mark_sheet_information(M* mesh, sheet* sh);
		void sheet_collapse(M* mesh, sheet*s);
		void get_sheet_faces();
		void inflate_faecs_to_sheet();


		/*column operation*/
		void get_column();
		void column_collapse(M* mesh, column* c);
		void get_column_faces();
		void inflate_faces_to_column();

		/*constrcut elementes*/
		V* new_vertex(M* mesh);
		E* new_edge(M* mesh);
		F* new_face(M* mesh);
		H* new_hex(M* mesh);

	//	/*id*/
	//	int& sheet_num() { return sheet_num; };
	//	int& column_num() { return column_num; };

	//private:
	//	/*id*/
	//	int sheet_num = 0;
	//	int column_num = 0;

	
	};

	template<typename V, typename E, typename F, typename H, typename M>
	void topoOperation<V, E, F, H, M>::both_end_singularity_insert(M* mesh, both_end_section* s, int insert_layer_num)
	{
		if (insert_layer_num <= 1)
		{
			std::cout << "insert_layer_num must bigger than 1" << std::endl;
			return;
		}
		mesh->mark_singularity();
		//���section�ϵ���
		for (int fIndex = 0; fIndex < s->faces().size(); fIndex++)
		{
			F* f = s->faces()[fIndex];
			f->singularity_insert() = true;
		}
		//��ǹ̶���
		for (int fixEIndex = 0; fixEIndex < s->fixed_edges1().size(); fixEIndex++)
		{
			E* e = s->fixed_edges1()[fixEIndex];
			for (int i = 0; i < 2; i++)
			{
				mesh->idVertices(e->vs[i])->fixed() = true;
			}
			e->fixed() = true;
		}
		for (int fixEIndex = 0; fixEIndex < s->fixed_edges2().size(); fixEIndex++)
		{
			E* e = s->fixed_edges2()[fixEIndex];
			for (int i = 0; i < 2; i++)
			{
				mesh->idVertices(e->vs[i])->fixed() = true;
			}
			e->fixed() = true;
		}
		int num_of_layer = both_end_num_of_layers(mesh, s);//��Ĳ���
		if (insert_layer_num > num_of_layer||insert_layer_num*2+1>num_of_layer)
		{
			std::cout << "the insert layer number exceeds the range of num of layer " << std::endl;
			return;
		}

		//ȡ�����㣬��
		std::vector<E*> oldEs;
		std::vector<V*> oldVs;
		std::unordered_set<H*> sidesHs;//�������hexs

		//ÿ�㴴���µĵ�
		int newvNum = 1;
		for (int layerNum = 1; layerNum < s->layers_vertieces().size()-1; layerNum++)
		{
			std::unordered_set<V*> layerVs = s->layers_vertieces()[layerNum];
			if (layerNum < insert_layer_num)
			{
				newvNum = layerNum * 2;
			}
			if (layerNum > (num_of_layer - insert_layer_num))
			{
				newvNum = newvNum - 2;
			}
			for (std::unordered_set<V*>::iterator vite = layerVs.begin(); vite != layerVs.end(); vite++)
			{
				V* v = *vite;
				//v->_from_string();
				oldVs.push_back(v);
				for (int i = 0; i < newvNum; i++)
				{
					V* newv = new V();
					mesh->maxVid()++;
					newv->id() = mesh->maxVid();
					mesh->vs.push_back(newv);
					mesh->m_map_vertices.insert(std::pair<int, V*>(newv->id(), newv));
					newv->position() = v->position();
					v->newv().push_back(newv);
					newv->oldv() = v;

					//��ֵ����
					newv->tail() = v->tail();
					newv->pBoundary() = v->pBoundary();
					newv->bctype() = v->bctype();
					newv->slit() = v->slit();
					newv->boundary() = v->boundary();
					newv->tail_farfield() = v->tail_farfield();
					newv->slit_surface() = v->slit_surface();

				}
			}
		}

		//�ڱߣ����ϴ����µ���Ϣ
		for (int fIndex = 0; fIndex < s->faces().size(); fIndex++)
		{
			F* f = s->faces()[fIndex];
			//��ֵsidesHs
			for (int hIndex = 0; hIndex < f->neighbor_hs.size(); hIndex++)
			{
				sidesHs.insert(mesh->idHexs(f->neighbor_hs[hIndex]));
			}

			std::vector<V*> fvs;
			std::vector<E*> fes;
			//�����µ�
			for (int fvIndex = 0; fvIndex < f->vs.size(); fvIndex++)
			{
				V* v = mesh->idVertices(f->vs[fvIndex]);
				fvs.push_back(v);
			}
			//�����±�
			for (int feIndex = 0; feIndex < f->es.size(); feIndex++)
			{
				E* e = mesh->idEdges(f->es[feIndex]);
				e->singularity_insert() = true;
				fes.push_back(e);
				//ÿ�����ϴ���һ���µı�
				for (int eIndex = 0; eIndex < f->es.size(); eIndex++)
				{
					E* e = mesh->idEdges(f->es[eIndex]);
					if (!e->fixed())
					{
						if (e->newe().size() == 0)
						{
							oldEs.push_back(e);
							for (int i = 0; i < f->neighbor_hs.size() - 1; i++)
							{
								E* newe = new E();
								mesh->maxEid()++;
								newe->id() = mesh->maxEid();
								mesh->es.push_back(newe);
								mesh->m_map_edges.insert(std::pair<int, E*>(newe->id(), newe));

								V* ev1 = mesh->idVertices(e->vs[0]);
								newe->vs.push_back(ev1->id());

								V* ev2 = mesh->idVertices(e->vs[1]);
								newe->vs.push_back(ev2->id());

								e->newe().push_back(newe);
								e->singularity_insert() = true;
							}

						}
					}
				}

			}

			//��������	
			for (int i = 0; i < f->neighbor_hs.size() - 1; i++)
			{
				//����
				F* newf = new F();
				newf->singularity_insert() = true;
				mesh->maxFid()++;
				newf->id() = mesh->maxFid();
				mesh->fs.push_back(newf);
				mesh->m_map_faces.insert(std::pair<int, F*>(newf->id(), newf));
				for (int vIndex = 0; vIndex < 4; vIndex++)
				{
					newf->vs.push_back(fvs[vIndex]->id());

					if (fvs[vIndex]->fixed())
					{
						fvs[vIndex]->neighbor_fs.push_back(newf->id());
					}

					newf->es.push_back(fes[vIndex]->id());
					if (fes[vIndex]->fixed())
					{
						fes[vIndex]->neighbor_fs.push_back(newf->id());
					}
				}
				f->newf().push_back(newf);
				newf->oldf() = f;
			}
			//�滻��
			std::vector<int> adjHs = f->neighbor_hs;
			for (int adjhIndex = 0; adjhIndex < (adjHs.size() - 1); adjhIndex++)
			{
				H* h = mesh->idHexs(adjHs[adjhIndex]);
				mesh->change_Hex_F(h, f, f->newf()[adjhIndex]);
			}
		}

		//�滻��
		for (int eIndex = 0; eIndex < oldEs.size(); eIndex++)
		{
			E* e = oldEs[eIndex];
			if (e->fixed())
			{
				continue;
			}
			std::vector<H*> adjHs;
			//�ҵ�singularity_insert����
			F* traceF = NULL;
			for (int efIndex = 0; efIndex < e->neighbor_fs.size(); efIndex++)
			{
				F* adjf = mesh->idFaces(e->neighbor_fs[efIndex]);
				if (adjf->singularity_insert())
				{
					traceF = adjf->newf()[0];
					break;
				}
			}
			H* traceH = mesh->idHexs(traceF->neighbor_hs[0]);
			adjHs.push_back(traceH);

			//flipe�ҵ�һ�����е�hex
			while (true)
			{
				traceF = mesh->flip_f(traceH, traceF, e);
				if (traceF->singularity_insert() || traceF->neighbor_hs.size() == 1)
				{
					break;
				}
				traceH = mesh->flip_H(traceH, traceF);
				adjHs.push_back(traceH);
			}
			//�滻���еı�
			for (int adjhIndex = 0; adjhIndex < adjHs.size(); adjhIndex++)
			{
				H* adjH = adjHs[adjhIndex];
				E* targetE = e->newe()[0];
				mesh->change_Hex_E(adjH, e, targetE);
			}

		}

		//�滻��
		std::vector<V*> newvs;
		for (int vIndex = 0; vIndex < oldVs.size(); vIndex++)
		{
			V* v = oldVs[vIndex];
			if (v->fixed())
			{
				continue;
			}
			std::unordered_set<H*> adjHs;

			F* traceF = NULL;
			for (int adjfIndex = 0; adjfIndex < v->neighbor_fs.size(); adjfIndex++)
			{
				F* f = mesh->idFaces(v->neighbor_fs[adjfIndex]);
				if (f->singularity_insert())
				{
					if (f->newf().size() != 0)
					{
						traceF = f->newf()[0];
						break;
					}
				}
			}

			if (traceF != NULL)
			{
				E* traceE = NULL;
				for (int feIndex = 0; feIndex < 4; feIndex++)
				{
					E* e = mesh->idEdges(traceF->es[feIndex]);
					int evIndex = e->vertexIndex(v->id());
					if (evIndex != -1)
					{
						traceE = e;
						break;
					}
				}
				//trace�ڽӵ��壬�滻��
				H* traceH = mesh->idHexs(traceF->neighbor_hs[0]);
				F* startF = traceF;
				E* startE = traceE;
				H* startH = traceH;

				bool is_circle = false;
				adjHs.insert(traceH);
				//����trace
				while (true)
				{
					traceE = mesh->flip_e(traceF, traceE, v);
					if (traceE == startE)
					{
						is_circle = true;
						break;

					}
					bool is_break = false;
					while (true)
					{
						traceF = mesh->flip_f(traceH, traceF, traceE);
						if (traceF->singularity_insert())
						{
							break;
						}
						else if (traceF->neighbor_hs.size() == 1)
						{
							is_break = true;
							break;
						}
						traceH = mesh->flip_H(traceH, traceF);
						adjHs.insert(traceH);
					}
					if (is_break)
					{
						break;
					}

				}

				//�����Ȧ������trace
				if (!is_circle)
				{

					while (true)
					{
						bool is_break = false;
						while (true)
						{
							traceF = mesh->flip_f(traceH, traceF, traceE);
							if (traceF->singularity_insert())
							{
								break;
							}
							else if (traceF->neighbor_hs.size() == 1)
							{
								is_break = true;
								break;
							}
							traceH = mesh->flip_H(traceH, traceF);
							adjHs.insert(traceH);
						}
						if (is_break)
						{
							break;
						}
						traceE = mesh->flip_e(traceF, traceE, v);
					}

				}

				for (std::unordered_set<H*>::iterator hite = adjHs.begin(); hite != adjHs.end(); hite++)
				{
					H* h = *hite;
					mesh->change_Hex_V(h, v, v->newv()[0]);
					newvs.push_back(v->newv()[0]);
				}
			}
		}

		//��������ķ���
		mesh->Facenormal();
		for (std::unordered_set<H*>::iterator ite = sidesHs.begin(); ite != sidesHs.end(); ite++)
		{
			H* h = *ite;
			mesh->reviese_hex_face_normal(h);
		}

		//Ų���������
		std::vector<V*> move_vs = oldVs;
		move_vs.insert(move_vs.end(), newvs.begin(), newvs.end());
		for (int vIndex = 0; vIndex < move_vs.size(); vIndex++)
		{
			V* v = move_vs[vIndex];
			if (v->fixed())
			{
				continue;
			}
			//������̵ıߵľ���
			double min_dis = 0;
			for (int adjeIndex = 0; adjeIndex < v->neighbor_es.size(); adjeIndex++)
			{
				E* adje = mesh->idEdges(v->neighbor_es[adjeIndex]);
				CPoint p1 = mesh->idVertices(adje->vs[0])->position();
				CPoint p2 = mesh->idVertices(adje->vs[1])->position();
				CPoint p3 = p1 - p2;
				double dis = sqrt(p3[0] * p3[0] + p3[1] * p3[1] + p3[2] * p3[2]);
				if (adjeIndex == 0)
				{
					min_dis = dis;
				}
				else
				{
					min_dis = min_dis > dis ? dis : min_dis;
				}
			}
			min_dis /= 2;
			CPoint v_normal;
			for (int adjfIndex = 0; adjfIndex < v->neighbor_fs.size(); adjfIndex++)
			{
				F* adjf = mesh->idFaces(v->neighbor_fs[adjfIndex]);
				v_normal += adjf->normal();
			}
			v_normal /= v->neighbor_fs.size();
			v->position() += (-v_normal) * min_dis;
		}
		
		//Ų����ֵ�ĵ�
		for (int vIndex = 0; vIndex < oldVs.size(); vIndex++)
		{
			V* v = oldVs[vIndex];
			if (v->fixed())
			{
				continue;
			}
			CPoint p1 = v->newv()[0]->position();
			CPoint p2 = v->position();
			CPoint vec = p1 - p2;
			CPoint step = vec / v->newv().size();
			for (int newvIndex = 1; newvIndex < v->newv().size(); newvIndex++)
			{
				V* newv = v->newv()[newvIndex];
				newv->position() = p2 + step * newvIndex;
			}
		}
		
		//���ڲ����hex
		std::vector<H*> newhs;
		for (int fixEIndex = 0; fixEIndex < s->fixed_edges1().size(); fixEIndex++)
		{
			//�ҵ��̶�������һ��ıߣ�����hex
			E* e = s->fixed_edges1()[fixEIndex];
			F* f = NULL;
			for (int efIndex = 0; efIndex < e->neighbor_fs.size(); efIndex++)
			{
				F* adjf = mesh->idFaces(e->neighbor_fs[efIndex]);
				if (adjf->singularity_insert())
				{
					f = adjf;
				}
			}

			//�ҵ����������İ˸���
			int hvs[8];
			for (int vIndex = 0; vIndex < 4; vIndex++)
			{
				V* v = mesh->idVertices(f->vs[vIndex]);
				if (v->fixed())
				{
					V* adjv = NULL;
					for (M::VVInFaceIterator vite(mesh, f, v); !vite.end(); vite++)
					{
						V* tempv = *vite;
						if (!tempv->fixed())
						{
							adjv = tempv;
							break;
						}
					}
					if (adjv->newv().size() == 0)
					{
						hvs[vIndex + 4] = adjv->oldv()->id();
						
					}
					else
					{
						hvs[vIndex + 4] = adjv->newv()[0]->id();
						
					}


				}
				else
				{
					if (v->newv().size() == 0)
					{
						hvs[vIndex + 4] = v->oldv()->newv()[1]->id();
						
					}
					else
					{
						hvs[vIndex + 4] = v->newv()[1]->id();
						
					}

				}
				hvs[vIndex] = v->id();
				
			}

			H* newh = new H();
			mesh->maxHid()++;
			mesh->m_map_hexs.insert(std::pair<int, H*>(mesh->maxHid(), newh));
			mesh->hs.push_back(newh);
			mesh->construct_hex(newh, mesh->maxHid(), hvs);
			newhs.push_back(newh);

			//���µ�hex�����ӱ��
			bool hex_slit = true;
			for (int hvIndex = 0; hvIndex < 8; hvIndex++)
			{
				V* v = mesh->idVertices(hvs[hvIndex]);
				if (!v->slit())
				{
					hex_slit = false;
					break;
				}
			}
			newh->fix() = hex_slit;
			
			
			

			//����ϲ�
			int hex_layer_num = 1;//������Ĳ���������������ȵ�Ĳ�������һ��
			E* traceE = mesh->e_parallel_e_in_face(f, e);
			F* traceF = f;

			while (true)
			{
				std::vector<F*> adj_insert_f;
				for (int efIndex = 0; efIndex < traceE->neighbor_fs.size(); efIndex++)
				{
					F* adjf = mesh->idFaces(traceE->neighbor_fs[efIndex]);
					if (adjf->singularity_insert())
					{
						adj_insert_f.push_back(adjf);
					}
				}
				if ((num_of_layer- insert_layer_num) == hex_layer_num)
				{
					break;
				}
				/*if (adj_insert_f.size() == 1)
				{
					break;
				}*/
				traceF = adj_insert_f[0]->id() == traceF->id() ? adj_insert_f[1] : adj_insert_f[0];
				hex_layer_num++;

				if (hex_layer_num < insert_layer_num)
				{
					//��ǰ�����������������
					int h_num = (hex_layer_num - 1) * 2 + 1;
					int half_num = h_num / 2;

					//����ǰ�벿�ֵ�hex
					for (int hexNum = 1; hexNum <= half_num; hexNum++)
					{
						//����hex
						for (int vIndex = 0; vIndex < 4; vIndex++)
						{
							V* v = mesh->idVertices(traceF->vs[vIndex]);
							if (hexNum == 1)
							{
								hvs[vIndex] = v->id();
								if (v->newv().size() == 0)
								{
									hvs[vIndex + 4] = v->oldv()->newv()[hexNum]->id();
								}
								else
								{
									hvs[vIndex + 4] = v->newv()[hexNum]->id();
								}
							}
							else
							{
								if (v->newv().size() == 0)
								{
									hvs[vIndex] = v->oldv()->newv()[hexNum - 1]->id();
									hvs[vIndex + 4] = v->oldv()->newv()[hexNum]->id();
								}
								else
								{
									hvs[vIndex] = v->newv()[hexNum - 1]->id();
									hvs[vIndex + 4] = v->newv()[hexNum]->id();
								}
							}
						}
						H* newh1 = new H();
						mesh->maxHid()++;
						mesh->m_map_hexs.insert(std::pair<int, H*>(mesh->maxHid(), newh1));
						mesh->hs.push_back(newh1);
						mesh->construct_hex(newh1, mesh->maxHid(), hvs);
						newhs.push_back(newh1);
						
						//���µ�hex�����ӱ��
						hex_slit = true;
						for (int hvIndex = 0; hvIndex < 8; hvIndex++)
						{
							V* v = mesh->idVertices(hvs[hvIndex]);
							if (!v->slit())
							{
								hex_slit = false;
								break;
							}
						}
						newh1->fix() = hex_slit;
						
					}

					//��������㲿�ֵ�hex
					for (int vIndex = 0; vIndex < 4; vIndex++)
					{
						V* v = mesh->idVertices(traceF->vs[vIndex]);
						if (v->newv().size() == 0)
						{
							if (v->oldv()->newv().size() == h_num - 1)//����ĵ�
							{
								hvs[vIndex] = v->oldv()->newv()[half_num]->id();
								//�ҵ���һ��ĵ㣬index��half_num+2;
								for (M::VVInFaceIterator vvite(mesh, traceF, v); !vvite.end(); vvite++)
								{
									V* adjv = *vvite;
									if (adjv->newv().size() != 0 && adjv->newv().size() != h_num - 1)
									{
										hvs[vIndex + 4] = adjv->newv()[half_num + 2]->id();
										break;
									}
									else if (adjv->newv().size() == 0 && adjv->oldv()->newv().size() != h_num - 1)
									{
										hvs[vIndex + 4] = adjv->oldv()->newv()[half_num + 2]->id();
										break;
									}
								}
							}
							else//����ĵ�
							{

								hvs[vIndex] = v->oldv()->newv()[half_num]->id();
								hvs[vIndex + 4] = v->oldv()->newv()[half_num + 1]->id();

							}
						}
						else
						{
							if (v->newv().size() == h_num - 1)//����ĵ�
							{
								hvs[vIndex] = v->newv()[half_num]->id();
								//�ҵ���һ��ĵ㣬index��half_num+2;
								for (M::VVInFaceIterator vvite(mesh, traceF, v); !vvite.end(); vvite++)
								{
									V* adjv = *vvite;
									if (adjv->newv().size() != 0 && adjv->newv().size() != h_num - 1)
									{
										hvs[vIndex + 4] = adjv->newv()[half_num + 2]->id();
										break;
									}
									else if (adjv->newv().size() == 0 && adjv->oldv()->newv().size() != h_num - 1)
									{
										hvs[vIndex + 4] = adjv->oldv()->newv()[half_num + 2]->id();
										break;
									}
								}
							}
							else
							{
								hvs[vIndex] = v->newv()[half_num]->id();
								hvs[vIndex + 4] = v->newv()[half_num + 1]->id();
							}
						}
					}
					H* singularityH = new H();
					mesh->maxHid()++;
					mesh->m_map_hexs.insert(std::pair<int, H*>(mesh->maxHid(), singularityH));
					mesh->hs.push_back(singularityH);
					mesh->construct_hex(singularityH, mesh->maxHid(), hvs);
					newhs.push_back(singularityH);

					//���µ�hex�����ӱ��
					hex_slit = true;
					for (int hvIndex = 0; hvIndex < 8; hvIndex++)
					{
						V* v = mesh->idVertices(hvs[hvIndex]);
						if (!v->slit())
						{
							hex_slit = false;
							break;
						}
					}
					singularityH->fix() = hex_slit;
					

					//������벿�ֵ�hex
					for (int hexNum = half_num + 2; hexNum <= h_num; hexNum++)
					{
						//����hex
						for (int vIndex = 0; vIndex < 4; vIndex++)
						{
							V* v = mesh->idVertices(traceF->vs[vIndex]);
							
							//�������һ��
							if (hexNum == h_num)
							{
								if (v->newv().size() == 0)
								{

									if (v->oldv()->newv().size() == h_num - 1)//����һ��
									{
										hvs[vIndex] = v->oldv()->newv()[hexNum - 2]->id();
										hvs[vIndex + 4] = v->oldv()->id();
									}
									else//����һ��
									{

										hvs[vIndex] = v->oldv()->newv()[hexNum]->id();
										hvs[vIndex + 4] = v->oldv()->id();
									}
								}
								else
								{
									if (v->newv().size() == h_num - 1)//����һ��
									{
										hvs[vIndex] = v->newv()[hexNum - 2]->id();
										hvs[vIndex + 4] = v->newv()[0]->id();
									}
									else
									{
										hvs[vIndex] = v->newv()[hexNum]->id();
										hvs[vIndex + 4] = v->newv()[0]->id();
									}
								}

							}
							else
							{
								if (v->newv().size() == 0)
								{
									if (v->oldv()->newv().size() == h_num - 1)//����һ��
									{
										hvs[vIndex] = v->oldv()->newv()[hexNum - 2]->id();
										hvs[vIndex + 4] = v->oldv()->newv()[hexNum - 1]->id();
									}
									else
									{
										hvs[vIndex] = v->oldv()->newv()[hexNum]->id();
										hvs[vIndex + 4] = v->oldv()->newv()[hexNum + 1]->id();
									}
								}
								else
								{
									if (v->newv().size() == h_num - 1)//����һ��
									{
										hvs[vIndex] = v->newv()[hexNum - 2]->id();
										hvs[vIndex + 4] = v->newv()[hexNum - 1]->id();
									}
									else
									{
										hvs[vIndex] = v->newv()[hexNum]->id();
										hvs[vIndex + 4] = v->newv()[hexNum + 1]->id();
									}
								}
							}
						}
						H* newh1 = new H();
						mesh->maxHid()++;
						
						mesh->m_map_hexs.insert(std::pair<int, H*>(mesh->maxHid(), newh1));
						mesh->hs.push_back(newh1);
						mesh->construct_hex(newh1, mesh->maxHid(), hvs);
						newhs.push_back(newh1);
						//���µ�hex�����ӱ��
						hex_slit = true;
						for (int hvIndex = 0; hvIndex < 8; hvIndex++)
						{
							V* v = mesh->idVertices(hvs[hvIndex]);
							if (!v->slit())
							{
								hex_slit = false;
								break;
							}
						}
						newh1->fix() = hex_slit;
					}
				}
				else
				{
					//��ǰ�����������������
					int h_num = (insert_layer_num - 1) * 2;
					for (int hexNum = 1; hexNum <= h_num; hexNum++)
					{
						for (int vIndex = 0; vIndex < 4; vIndex++)
						{
							V* v = mesh->idVertices(traceF->vs[vIndex]);
							if (hexNum == 1)
							{
								hvs[vIndex] = v->id();
								if (v->newv().size() == 0)
								{
									hvs[vIndex + 4] = v->oldv()->newv()[hexNum]->id();
								}
								else
								{
									hvs[vIndex + 4] = v->newv()[hexNum]->id();
								}
							}
							else if (hexNum == h_num)
							{
								if (v->newv().size() == 0)
								{
									hvs[vIndex] = v->oldv()->newv()[hexNum - 1]->id();
									hvs[vIndex + 4] = v->oldv()->id();
								}
								else
								{
									hvs[vIndex] = v->newv()[hexNum - 1]->id();
									hvs[vIndex + 4] = v->newv()[0]->id();
								}
							}
							else
							{
								if (v->newv().size() == 0)
								{
									hvs[vIndex] = v->oldv()->newv()[hexNum - 1]->id();
									hvs[vIndex + 4] = v->oldv()->newv()[hexNum]->id();;
								}
								else
								{
									hvs[vIndex] = v->newv()[hexNum - 1]->id();
									hvs[vIndex + 4] = v->newv()[hexNum]->id();
								}
							}
						}
						H* newh1 = new H();
						mesh->maxHid()++;
						mesh->m_map_hexs.insert(std::pair<int, H*>(mesh->maxHid(), newh1));
						mesh->hs.push_back(newh1);
						mesh->construct_hex(newh1, mesh->maxHid(), hvs);
						newhs.push_back(newh1);
						//���µ�hex�����ӱ��
						hex_slit = true;
						for (int hvIndex = 0; hvIndex < 8; hvIndex++)
						{
							V* v = mesh->idVertices(hvs[hvIndex]);
							if (!v->slit())
							{
								hex_slit = false;
								break;
							}
						}
						newh1->fix() = hex_slit;
					}
				}

				traceE = mesh->e_parallel_e_in_face(traceF, traceE);

			}

		}
		
		//���fixed_edge2�ϵĵ㣬��䵽�����
		for (int fixEIndex = 0; fixEIndex < s->fixed_edges2().size(); fixEIndex++)
		{
			//�ҵ��̶�������һ��ıߣ�����hex
			E* e = s->fixed_edges2()[fixEIndex];
			F* f = NULL;
			for (int efIndex = 0; efIndex < e->neighbor_fs.size(); efIndex++)
			{
				F* adjf = mesh->idFaces(e->neighbor_fs[efIndex]);
				if (adjf->singularity_insert())
				{
					f = adjf;
				}
			}

			//�ҵ����������İ˸���
			int hvs[8];
			for (int vIndex = 0; vIndex < 4; vIndex++)
			{
				V* v = mesh->idVertices(f->vs[vIndex]);
				if (v->fixed())
				{
					V* adjv = NULL;
					for (M::VVInFaceIterator vite(mesh, f, v); !vite.end(); vite++)
					{
						V* tempv = *vite;
						if (!tempv->fixed())
						{
							adjv = tempv;
							break;
						}
					}
					if (adjv->newv().size() == 0)
					{
						hvs[vIndex + 4] = adjv->oldv()->id();
					}
					else
					{
						hvs[vIndex + 4] = adjv->newv()[0]->id();
					}


				}
				else
				{
					if (v->newv().size() == 0)
					{
						hvs[vIndex + 4] = v->oldv()->newv()[1]->id();
					}
					else
					{
						hvs[vIndex + 4] = v->newv()[1]->id();
					}

				}
				hvs[vIndex] = v->id();
			}

			H* newh = new H();
			mesh->maxHid()++;
			mesh->m_map_hexs.insert(std::pair<int, H*>(mesh->maxHid(), newh));
			mesh->hs.push_back(newh);
			mesh->construct_hex(newh, mesh->maxHid(), hvs);
			newhs.push_back(newh);
			//���µ�hex�����ӱ��
			bool hex_slit = true;
			for (int hvIndex = 0; hvIndex < 8; hvIndex++)
			{
				V* v = mesh->idVertices(hvs[hvIndex]);
				if (!v->slit())
				{
					hex_slit = false;
					break;
				}
			}
			newh->fix() = hex_slit;

			//����ϲ�
			int hex_layer_num = 1;//������Ĳ���������������ȵ�Ĳ�������һ��
			E* traceE = mesh->e_parallel_e_in_face(f, e);
			F* traceF = f;

			while (true)
			{
				std::vector<F*> adj_insert_f;
				for (int efIndex = 0; efIndex < traceE->neighbor_fs.size(); efIndex++)
				{
					F* adjf = mesh->idFaces(traceE->neighbor_fs[efIndex]);
					if (adjf->singularity_insert())
					{
						adj_insert_f.push_back(adjf);
					}
				}

				/*if (adj_insert_f.size() == 1)
				{
					break;
				}*/
				traceF = adj_insert_f[0]->id() == traceF->id() ? adj_insert_f[1] : adj_insert_f[0];
				hex_layer_num++;

				if (hex_layer_num < insert_layer_num)
				{
					//��ǰ�����������������
					int h_num = (hex_layer_num - 1) * 2 + 1;
					int half_num = h_num / 2;

					//����ǰ�벿�ֵ�hex
					for (int hexNum = 1; hexNum <= half_num; hexNum++)
					{
						//����hex
						for (int vIndex = 0; vIndex < 4; vIndex++)
						{
							V* v = mesh->idVertices(traceF->vs[vIndex]);
							if (hexNum == 1)
							{
								hvs[vIndex] = v->id();
								if (v->newv().size() == 0)
								{
									hvs[vIndex + 4] = v->oldv()->newv()[hexNum]->id();
								}
								else
								{
									hvs[vIndex + 4] = v->newv()[hexNum]->id();
								}
							}
							else
							{
								if (v->newv().size() == 0)
								{
									hvs[vIndex] = v->oldv()->newv()[hexNum - 1]->id();
									hvs[vIndex + 4] = v->oldv()->newv()[hexNum]->id();
								}
								else
								{
									hvs[vIndex] = v->newv()[hexNum - 1]->id();
									hvs[vIndex + 4] = v->newv()[hexNum]->id();
								}
							}
						}
						H* newh1 = new H();
						mesh->maxHid()++;
						mesh->m_map_hexs.insert(std::pair<int, H*>(mesh->maxHid(), newh1));
						mesh->hs.push_back(newh1);
						mesh->construct_hex(newh1, mesh->maxHid(), hvs);
						newhs.push_back(newh1);
						//���µ�hex�����ӱ��
						hex_slit = true;
						for (int hvIndex = 0; hvIndex < 8; hvIndex++)
						{
							V* v = mesh->idVertices(hvs[hvIndex]);
							if (!v->slit())
							{
								hex_slit = false;
								break;
							}
						}
						newh1->fix() = hex_slit;
					}

					//��������㲿�ֵ�hex
					for (int vIndex = 0; vIndex < 4; vIndex++)
					{
						V* v = mesh->idVertices(traceF->vs[vIndex]);
						if (v->newv().size() == 0)
						{
							if (v->oldv()->newv().size() == h_num - 1)//����ĵ�
							{
								hvs[vIndex] = v->oldv()->newv()[half_num]->id();
								//�ҵ���һ��ĵ㣬index��half_num+2;
								for (M::VVInFaceIterator vvite(mesh, traceF, v); !vvite.end(); vvite++)
								{
									V* adjv = *vvite;
									if (adjv->newv().size() != 0 && adjv->newv().size() != h_num - 1)
									{
										hvs[vIndex + 4] = adjv->newv()[half_num + 2]->id();
										break;
									}
									else if (adjv->newv().size() == 0 && adjv->oldv()->newv().size() != h_num - 1)
									{
										hvs[vIndex + 4] = adjv->oldv()->newv()[half_num + 2]->id();
										break;
									}
								}
							}
							else//����ĵ�
							{

								hvs[vIndex] = v->oldv()->newv()[half_num]->id();
								hvs[vIndex + 4] = v->oldv()->newv()[half_num + 1]->id();

							}
						}
						else
						{
							if (v->newv().size() == h_num - 1)//����ĵ�
							{
								hvs[vIndex] = v->newv()[half_num]->id();
								//�ҵ���һ��ĵ㣬index��half_num+2;
								for (M::VVInFaceIterator vvite(mesh, traceF, v); !vvite.end(); vvite++)
								{
									V* adjv = *vvite;
									if (adjv->newv().size() != 0 && adjv->newv().size() != h_num - 1)
									{
										hvs[vIndex + 4] = adjv->newv()[half_num + 2]->id();
										break;
									}
									else if (adjv->newv().size() == 0 && adjv->oldv()->newv().size() != h_num - 1)
									{
										hvs[vIndex + 4] = adjv->oldv()->newv()[half_num + 2]->id();
										break;
									}
								}
							}
							else
							{
								hvs[vIndex] = v->newv()[half_num]->id();
								hvs[vIndex + 4] = v->newv()[half_num + 1]->id();
							}
						}
					}
					H* singularityH = new H();
					mesh->maxHid()++;
					mesh->m_map_hexs.insert(std::pair<int, H*>(mesh->maxHid(), singularityH));
					mesh->hs.push_back(singularityH);
					mesh->construct_hex(singularityH, mesh->maxHid(), hvs);
					newhs.push_back(singularityH);
					//���µ�hex�����ӱ��
					hex_slit = true;
					for (int hvIndex = 0; hvIndex < 8; hvIndex++)
					{
						V* v = mesh->idVertices(hvs[hvIndex]);
						if (!v->slit())
						{
							hex_slit = false;
							break;
						}
					}
					singularityH->fix() = hex_slit;


					//������벿�ֵ�hex
					for (int hexNum = half_num + 2; hexNum <= h_num; hexNum++)
					{
						//����hex
						for (int vIndex = 0; vIndex < 4; vIndex++)
						{
							V* v = mesh->idVertices(traceF->vs[vIndex]);
							
							//�������һ��
							if (hexNum == h_num)
							{
								if (v->newv().size() == 0)
								{

									if (v->oldv()->newv().size() == h_num - 1)//����һ��
									{
										hvs[vIndex] = v->oldv()->newv()[hexNum - 2]->id();
										hvs[vIndex + 4] = v->oldv()->id();
									}
									else//����һ��
									{

										hvs[vIndex] = v->oldv()->newv()[hexNum]->id();
										hvs[vIndex + 4] = v->oldv()->id();
									}
								}
								else
								{
									if (v->newv().size() == h_num - 1)//����һ��
									{
										hvs[vIndex] = v->newv()[hexNum - 2]->id();
										hvs[vIndex + 4] = v->newv()[0]->id();
									}
									else
									{
										hvs[vIndex] = v->newv()[hexNum]->id();
										hvs[vIndex + 4] = v->newv()[0]->id();
									}
								}

							}
							else
							{
								if (v->newv().size() == 0)
								{
									if (v->oldv()->newv().size() == h_num - 1)//����һ��
									{
										hvs[vIndex] = v->oldv()->newv()[hexNum - 2]->id();
										hvs[vIndex + 4] = v->oldv()->newv()[hexNum - 1]->id();
									}
									else
									{
										hvs[vIndex] = v->oldv()->newv()[hexNum]->id();
										hvs[vIndex + 4] = v->oldv()->newv()[hexNum + 1]->id();
									}
								}
								else
								{
									if (v->newv().size() == h_num - 1)//����һ��
									{
										hvs[vIndex] = v->newv()[hexNum - 2]->id();
										hvs[vIndex + 4] = v->newv()[hexNum - 1]->id();
									}
									else
									{
										hvs[vIndex] = v->newv()[hexNum]->id();
										hvs[vIndex + 4] = v->newv()[hexNum + 1]->id();
									}
								}
							}
						}
						H* newh1 = new H();
						mesh->maxHid()++;
						
						mesh->m_map_hexs.insert(std::pair<int, H*>(mesh->maxHid(), newh1));
						mesh->hs.push_back(newh1);
						mesh->construct_hex(newh1, mesh->maxHid(), hvs);
						newhs.push_back(newh1);

						//���µ�hex�����ӱ��
						hex_slit = true;
						for (int hvIndex = 0; hvIndex < 8; hvIndex++)
						{
							V* v = mesh->idVertices(hvs[hvIndex]);
							if (!v->slit())
							{
								hex_slit = false;
								break;
							}
						}
						newh1->fix() = hex_slit;

					}
				}
				else
				{
					break;
				}
			}
		}
		//�޸�����
		for (int hIndex = 0; hIndex < newhs.size(); hIndex++)
		{
			H* h = newhs[hIndex];
			mesh->reviese_hex_face_normal(h);
		}

	}
	
	template<typename V, typename E, typename F, typename H, typename M>
	void topoOperation<V, E, F, H, M>::singularity_insert(M* mesh, section* s, int insert_layer_num)
	{
		if (insert_layer_num <=1)
		{
			std::cout << "insert_layer_num must bigger than 1" << std::endl;
			return;
		}
		clock_t start = clock();

		mesh->mark_singularity();
		//���section�ϵ���
		for (int fIndex = 0; fIndex < s->faces().size(); fIndex++)
		{
			F* f = s->faces()[fIndex];
			f->singularity_insert() = true;
		}
		//��ǹ̶���
		for (int fixEIndex = 0; fixEIndex < s->fixed_edges().size(); fixEIndex++)
		{
			E* e = s->fixed_edges()[fixEIndex];
			for (int i = 0; i < 2; i++)
			{
				mesh->idVertices(e->vs[i])->fixed() = true;
			}
			e->fixed() = true;
		}

		int num_of_layer = num_of_layers(mesh,s);//��Ĳ���
		if (insert_layer_num> num_of_layer)
		{
			std::cout << "the insert layer number exceeds the range of num of layer " << std::endl;
			return;
		}

		//ȡ�����㣬��
		std::vector<E*> oldEs;
		std::vector<V*> oldVs;
		std::unordered_set<H*> sidesHs;//�������hexs

		//ÿ�㴴���µĵ�
		int newvNum = 1;
		for (int layerNum = 1; layerNum < s->layers_vertieces().size(); layerNum++)
		{
			std::unordered_set<V*> layerVs = s->layers_vertieces()[layerNum];
			if (layerNum<insert_layer_num)
			{
				newvNum = layerNum * 2;
			}	
			for (std::unordered_set<V*>::iterator vite = layerVs.begin();vite!=layerVs.end();vite++ )
			{
				V* v = *vite;
				oldVs.push_back(v);
				for (int i = 0; i < newvNum; i++)
				{
					V* newv = new V();
					mesh->maxVid()++;
					newv->id() = mesh->maxVid();
					mesh->vs.push_back(newv);
					mesh->m_map_vertices.insert(std::pair<int, V*>(newv->id(), newv));
					newv->position() = v->position();
					v->newv().push_back(newv);
					newv->oldv() = v;
				}
			}
		}

		//�ڱߣ����ϴ����µ���Ϣ
		for (int fIndex = 0; fIndex < s->faces().size(); fIndex++)
		{
			F* f = s->faces()[fIndex];
			//��ֵsidesHs
			for (int hIndex = 0; hIndex < f->neighbor_hs.size(); hIndex++)
			{
				sidesHs.insert(mesh->idHexs(f->neighbor_hs[hIndex]));
			}

			std::vector<V*> fvs;
			std::vector<E*> fes;
			//�����µ�
			for (int fvIndex = 0; fvIndex < f->vs.size(); fvIndex++)
			{
				V* v = mesh->idVertices(f->vs[fvIndex]);
				fvs.push_back(v);
			}
			//�����±�
			for (int feIndex = 0; feIndex < f->es.size(); feIndex++)
			{
				E* e = mesh->idEdges(f->es[feIndex]);
				e->singularity_insert() = true;
				fes.push_back(e);
				//ÿ�����ϴ���һ���µı�
				for (int eIndex = 0; eIndex < f->es.size(); eIndex++)
				{
					E* e = mesh->idEdges(f->es[eIndex]);
					if (!e->fixed())
					{
						if (e->newe().size() == 0)
						{
							oldEs.push_back(e);
							for (int i = 0; i < f->neighbor_hs.size() - 1; i++)
							{
								E* newe = new E();
								mesh->maxEid()++;
								newe->id() = mesh->maxEid();
								mesh->es.push_back(newe);
								mesh->m_map_edges.insert(std::pair<int, E*>(newe->id(), newe));

								V* ev1 = mesh->idVertices(e->vs[0]);
								newe->vs.push_back(ev1->id());

								V* ev2 = mesh->idVertices(e->vs[1]);
								newe->vs.push_back(ev2->id());

								e->newe().push_back(newe);
								e->singularity_insert() = true;
							}

						}
					}
				}

			}

			//��������	
			for (int i = 0; i < f->neighbor_hs.size() - 1; i++)
			{
				//����
				F* newf = new F();
				newf->singularity_insert() = true;
				mesh->maxFid()++;
				newf->id() = mesh->maxFid();
				mesh->fs.push_back(newf);
				mesh->m_map_faces.insert(std::pair<int, F*>(newf->id(), newf));
				for (int vIndex = 0; vIndex < 4; vIndex++)
				{
					newf->vs.push_back(fvs[vIndex]->id());

					if (fvs[vIndex]->fixed())
					{
						fvs[vIndex]->neighbor_fs.push_back(newf->id());
					}

					newf->es.push_back(fes[vIndex]->id());
					if (fes[vIndex]->fixed())
					{
						fes[vIndex]->neighbor_fs.push_back(newf->id());
					}
				}
				f->newf().push_back(newf);
				newf->oldf() = f;
			}
			//�滻��
			std::vector<int> adjHs = f->neighbor_hs;
			for (int adjhIndex = 0; adjhIndex < (adjHs.size() - 1); adjhIndex++)
			{
				H* h = mesh->idHexs(adjHs[adjhIndex]);
				mesh->change_Hex_F(h, f, f->newf()[adjhIndex]);
			}
		}

		//�滻��
		for (int eIndex = 0; eIndex < oldEs.size(); eIndex++)
		{
			E* e = oldEs[eIndex];
			if (e->fixed())
			{
				continue;
			}
			std::vector<H*> adjHs;
			//�ҵ�singularity_insert����
			F* traceF = NULL;
			for (int efIndex = 0; efIndex < e->neighbor_fs.size(); efIndex++)
			{
				F* adjf = mesh->idFaces(e->neighbor_fs[efIndex]);
				if (adjf->singularity_insert())
				{
					traceF = adjf->newf()[0];
					break;
				}
			}
			H* traceH = mesh->idHexs(traceF->neighbor_hs[0]);
			adjHs.push_back(traceH);

			//flipe�ҵ�һ�����е�hex
			while (true)
			{
				traceF = mesh->flip_f(traceH, traceF, e);
				if (traceF->singularity_insert() || traceF->neighbor_hs.size() == 1)
				{
					break;
				}
				traceH = mesh->flip_H(traceH, traceF);
				adjHs.push_back(traceH);
			}
			//�滻���еı�
			for (int adjhIndex = 0; adjhIndex < adjHs.size(); adjhIndex++)
			{
				H* adjH = adjHs[adjhIndex];
				E* targetE = e->newe()[0];
				mesh->change_Hex_E(adjH, e, targetE);
			}

		}

		//�滻��
		std::vector<V*> newvs;
		for (int vIndex = 0; vIndex < oldVs.size(); vIndex++)
		{
			V* v = oldVs[vIndex];
			if (v->fixed())
			{
				continue;
			}
			std::unordered_set<H*> adjHs;

			F* traceF = NULL;
			for (int adjfIndex = 0; adjfIndex < v->neighbor_fs.size(); adjfIndex++)
			{
				F* f = mesh->idFaces(v->neighbor_fs[adjfIndex]);
				if (f->singularity_insert())
				{
					if (f->newf().size() != 0)
					{
						traceF = f->newf()[0];
						break;
					}
				}
			}

			if (traceF != NULL)
			{
				E* traceE = NULL;
				for (int feIndex = 0; feIndex < 4; feIndex++)
				{
					E* e = mesh->idEdges(traceF->es[feIndex]);
					int evIndex = e->vertexIndex(v->id());
					if (evIndex != -1)
					{
						traceE = e;
						break;
					}
				}
				//trace�ڽӵ��壬�滻��
				H* traceH = mesh->idHexs(traceF->neighbor_hs[0]);
				F* startF = traceF;
				E* startE = traceE;
				H* startH = traceH;

				bool is_circle = false;
				adjHs.insert(traceH);
				//����trace
				while (true)
				{
					traceE = mesh->flip_e(traceF, traceE, v);
					if (traceE == startE)
					{
						is_circle = true;
						break;

					}
					bool is_break = false;
					while (true)
					{
						traceF = mesh->flip_f(traceH, traceF, traceE);
						if (traceF->singularity_insert())
						{
							break;
						}
						else if (traceF->neighbor_hs.size() == 1)
						{
							is_break = true;
							break;
						}
						traceH = mesh->flip_H(traceH, traceF);
						adjHs.insert(traceH);
					}
					if (is_break)
					{
						break;
					}

				}

				//�����Ȧ������trace
				if (!is_circle)
				{

					while (true)
					{
						bool is_break = false;
						while (true)
						{
							traceF = mesh->flip_f(traceH, traceF, traceE);
							if (traceF->singularity_insert())
							{
								break;
							}
							else if (traceF->neighbor_hs.size() == 1)
							{
								is_break = true;
								break;
							}
							traceH = mesh->flip_H(traceH, traceF);
							adjHs.insert(traceH);
						}
						if (is_break)
						{
							break;
						}
						traceE = mesh->flip_e(traceF, traceE, v);
					}

				}

				for (std::unordered_set<H*>::iterator hite = adjHs.begin(); hite != adjHs.end(); hite++)
				{
					H* h = *hite;
					mesh->change_Hex_V(h, v, v->newv()[0]);
					newvs.push_back(v->newv()[0]);
				}
			}
		}

		//��������ķ���
		mesh->Facenormal();
		for (std::unordered_set<H*>::iterator ite = sidesHs.begin(); ite != sidesHs.end(); ite++)
		{
			H* h = *ite;
			mesh->reviese_hex_face_normal(h);
		}

		//Ų���������
		std::vector<V*> move_vs = oldVs;
		move_vs.insert(move_vs.end(), newvs.begin(), newvs.end());
		for (int vIndex = 0; vIndex < move_vs.size(); vIndex++)
		{
			V* v = move_vs[vIndex];
			if (v->fixed())
			{
				continue;
			}
			//������̵ıߵľ���
			double min_dis = 0;
			for (int adjeIndex = 0; adjeIndex < v->neighbor_es.size(); adjeIndex++)
			{
				E* adje = mesh->idEdges(v->neighbor_es[adjeIndex]);
				CPoint p1 = mesh->idVertices(adje->vs[0])->position();
				CPoint p2 = mesh->idVertices(adje->vs[1])->position();
				CPoint p3 = p1 - p2;
				double dis = sqrt(p3[0] * p3[0] + p3[1] * p3[1] + p3[2] * p3[2]);
				if (adjeIndex == 0)
				{
					min_dis = dis;
				}
				else
				{
					min_dis = min_dis > dis ? dis : min_dis;
				}
			}
			min_dis /= 2;
			CPoint v_normal;
			for (int adjfIndex = 0; adjfIndex < v->neighbor_fs.size(); adjfIndex++)
			{
				F* adjf = mesh->idFaces(v->neighbor_fs[adjfIndex]);
				v_normal += adjf->normal();
			}
			v_normal /= v->neighbor_fs.size();
			v->position() += (-v_normal) * min_dis;
		}
		//Ų����ֵ�ĵ�
		for (int vIndex = 0; vIndex < oldVs.size(); vIndex++)
		{
			V* v = oldVs[vIndex];
			if (v->fixed())
			{
				continue;
			}
			CPoint p1 = v->newv()[0]->position();
			CPoint p2 = v->position();
			CPoint vec = p1 - p2;
			CPoint step = vec / v->newv().size();
			for (int newvIndex = 1; newvIndex < v->newv().size(); newvIndex++)
			{
				V* newv = v->newv()[newvIndex];
				newv->position() = p2 + step * newvIndex;
			}
		}


		//���ڲ����hex
		std::vector<H*> newhs;
		for (int fixEIndex = 0; fixEIndex < s->fixed_edges().size(); fixEIndex++)
		{
			//�ҵ��̶�������һ��ıߣ�����hex
			E* e = s->fixed_edges()[fixEIndex];
			F* f = NULL;
			for (int efIndex = 0; efIndex < e->neighbor_fs.size(); efIndex++)
			{
				F* adjf = mesh->idFaces(e->neighbor_fs[efIndex]);
				if (adjf->singularity_insert())
				{
					f = adjf;
				}
			}
			
			//�ҵ����������İ˸���
			int hvs[8];
			for (int vIndex = 0; vIndex < 4; vIndex++)
			{
				V* v = mesh->idVertices(f->vs[vIndex]);
				if (v->fixed())
				{
					V* adjv = NULL;
					for (M::VVInFaceIterator vite(mesh, f, v); !vite.end(); vite++)
					{
						V* tempv = *vite;
						if (!tempv->fixed())
						{
							adjv = tempv;
							break;
						}
					}
					if (adjv->newv().size() == 0)
					{
						hvs[vIndex + 4] = adjv->oldv()->id();
					}
					else
					{
						hvs[vIndex + 4] = adjv->newv()[0]->id();
					}


				}
				else
				{
					if (v->newv().size() == 0)
					{
						hvs[vIndex + 4] = v->oldv()->newv()[1]->id();
					}
					else
					{
						hvs[vIndex + 4] = v->newv()[1]->id();
					}

				}
				hvs[vIndex] = v->id();
			}

			H* newh = new H();
			mesh->maxHid()++;
			mesh->m_map_hexs.insert(std::pair<int, H*>(mesh->maxHid(), newh));
			mesh->hs.push_back(newh);
			mesh->construct_hex(newh, mesh->maxHid(), hvs);
			newhs.push_back(newh);

			//����ϲ�
			int hex_layer_num = 1;//������Ĳ���������������ȵ�Ĳ�������һ��
			E* traceE = mesh->e_parallel_e_in_face(f, e);
			F* traceF = f;

			while (true)
			{				
				std::vector<F*> adj_insert_f;
				for (int efIndex = 0; efIndex < traceE->neighbor_fs.size(); efIndex++)
				{
					F* adjf = mesh->idFaces(traceE->neighbor_fs[efIndex]);
					if (adjf->singularity_insert())
					{
						adj_insert_f.push_back(adjf);
					}
				}
				if (adj_insert_f.size() == 1)
				{
					break;
				}
				traceF = adj_insert_f[0]->id() == traceF->id() ? adj_insert_f[1] : adj_insert_f[0];
				hex_layer_num++;

				if (hex_layer_num <insert_layer_num)
				{
					//��ǰ�����������������
					int h_num = (hex_layer_num - 1) * 2 + 1;
					int half_num = h_num / 2;
					
					//����ǰ�벿�ֵ�hex
					for (int hexNum = 1; hexNum <= half_num; hexNum++)
					{
						//����hex
						for (int vIndex = 0; vIndex < 4; vIndex++)
						{
							V* v = mesh->idVertices(traceF->vs[vIndex]);
							if (hexNum == 1)
							{
								hvs[vIndex] = v->id();	
								if (v->newv().size() == 0)
								{
									hvs[vIndex + 4] = v->oldv()->newv()[hexNum]->id();
								}
								else
								{
									hvs[vIndex + 4] = v->newv()[hexNum]->id();
								}
							}
							else
							{
								if (v->newv().size() == 0)
								{
									hvs[vIndex] = v->oldv()->newv()[hexNum - 1]->id();
									hvs[vIndex + 4] = v->oldv()->newv()[hexNum]->id();
								}
								else
								{
									hvs[vIndex] = v->newv()[hexNum - 1]->id();
									hvs[vIndex + 4] = v->newv()[hexNum]->id();
								}						
							}			
						}
						H* newh1 = new H();
						mesh->maxHid()++;
						mesh->m_map_hexs.insert(std::pair<int, H*>(mesh->maxHid(), newh1));
						mesh->hs.push_back(newh1);
						mesh->construct_hex(newh1, mesh->maxHid(), hvs);
						newhs.push_back(newh1);
					}
					
					//��������㲿�ֵ�hex
					for (int vIndex = 0; vIndex < 4; vIndex++)
					{
						V* v = mesh->idVertices(traceF->vs[vIndex]);
						if (v->newv().size()==0)
						{
							if (v->oldv()->newv().size()== h_num-1)//����ĵ�
							{
								hvs[vIndex] = v->oldv()->newv()[half_num]->id();
								//�ҵ���һ��ĵ㣬index��half_num+2;
								for (M::VVInFaceIterator vvite(mesh,traceF,v);!vvite.end();vvite++)
								{
									V* adjv = *vvite;
									if (adjv->newv().size() != 0 && adjv->newv().size() != h_num - 1)
									{
										hvs[vIndex + 4] = adjv->newv()[half_num + 2]->id();
										break;
									}
									else if(adjv->newv().size()==0 && adjv->oldv()->newv().size()!=h_num-1)
									{
										hvs[vIndex + 4] = adjv->oldv()->newv()[half_num + 2]->id();
										break;
									}				
								}
							}
							else//����ĵ�
							{
								
								hvs[vIndex] = v->oldv()->newv()[half_num]->id();
								hvs[vIndex+4] = v->oldv()->newv()[half_num + 1]->id();
								
							}
						}
						else
						{
							if (v->newv().size() == h_num - 1)//����ĵ�
							{
								hvs[vIndex] = v->newv()[half_num]->id();
								//�ҵ���һ��ĵ㣬index��half_num+2;
								for (M::VVInFaceIterator vvite(mesh, traceF, v); !vvite.end(); vvite++)
								{
									V* adjv = *vvite;
									if (adjv->newv().size() != 0 && adjv->newv().size() != h_num - 1)
									{
										hvs[vIndex + 4] = adjv->newv()[half_num + 2]->id();
										break;
									}
									else if (adjv->newv().size() == 0 && adjv->oldv()->newv().size() != h_num - 1)
									{
										hvs[vIndex + 4] = adjv->oldv()->newv()[half_num + 2]->id();
										break;
									}
								}
							}
							else
							{
								hvs[vIndex] = v->newv()[half_num]->id();
								hvs[vIndex+4] = v->newv()[half_num + 1]->id();
							}
						}						
					}				
					H* singularityH = new H();
					mesh->maxHid()++;
					mesh->m_map_hexs.insert(std::pair<int, H*>(mesh->maxHid(), singularityH));
					mesh->hs.push_back(singularityH);
					mesh->construct_hex(singularityH, mesh->maxHid(), hvs);
					newhs.push_back(singularityH);
					
					//������벿�ֵ�hex
					for (int hexNum = half_num + 2; hexNum <= h_num; hexNum++)
					{
						//����hex
						for (int vIndex = 0; vIndex < 4; vIndex++)
						{
							V* v = mesh->idVertices(traceF->vs[vIndex]);
							
							//�������һ��
							if (hexNum==h_num)
							{
								if (v->newv().size() == 0)
								{
									
									if (v->oldv()->newv().size()==h_num-1)//����һ��
									{
										hvs[vIndex] = v->oldv()->newv()[hexNum - 2]->id();
										hvs[vIndex + 4] = v->oldv()->id();
									}
									else//����һ��
									{
										
										hvs[vIndex] = v->oldv()->newv()[hexNum]->id();
										hvs[vIndex + 4] = v->oldv()->id();
									}				
								}
								else
								{
									if (v->newv().size() == h_num - 1)//����һ��
									{
										hvs[vIndex] = v->newv()[hexNum - 2]->id();
										hvs[vIndex + 4] = v->newv()[0]->id();
									}
									else
									{
										hvs[vIndex] = v->newv()[hexNum]->id();
										hvs[vIndex + 4] = v->newv()[0]->id();
									}
								}
								
							}
							else
							{
								if (v->newv().size() == 0)
								{
									if (v->oldv()->newv().size() == h_num - 1)//����һ��
									{
										hvs[vIndex] = v->oldv()->newv()[hexNum - 2]->id();
										hvs[vIndex + 4] = v->oldv()->newv()[hexNum - 1]->id();
									}
									else
									{
										hvs[vIndex] = v->oldv()->newv()[hexNum]->id();
										hvs[vIndex + 4] = v->oldv()->newv()[hexNum+1]->id();
									}
								}
								else
								{
									if (v->newv().size() == h_num - 1)//����һ��
									{
										hvs[vIndex] = v->newv()[hexNum - 2]->id();
										hvs[vIndex + 4] = v->newv()[hexNum - 1]->id();
									}
									else
									{
										hvs[vIndex] = v->newv()[hexNum]->id();
										hvs[vIndex + 4] = v->newv()[hexNum+1]->id();
									}
								}
							}
						}
						H* newh1 = new H();
						mesh->maxHid()++;
						
						mesh->m_map_hexs.insert(std::pair<int, H*>(mesh->maxHid(), newh1));
						mesh->hs.push_back(newh1);
						mesh->construct_hex(newh1, mesh->maxHid(), hvs);
						newhs.push_back(newh1);
					}
				}
				else
				{
					//��ǰ�����������������
					int h_num = (insert_layer_num - 1) * 2;
					for (int hexNum = 1; hexNum <= h_num; hexNum++)
					{
						for (int vIndex = 0; vIndex < 4; vIndex++)
						{
							V* v = mesh->idVertices(traceF->vs[vIndex]);
							if (hexNum == 1)
							{
								hvs[vIndex] = v->id();
								if (v->newv().size()==0)
								{
									hvs[vIndex + 4] = v->oldv()->newv()[hexNum]->id();
								}
								else
								{
									hvs[vIndex + 4] = v->newv()[hexNum]->id();
								}
								
							}
							else if (hexNum == h_num)
							{
								if (v->newv().size() == 0)
								{			
									hvs[vIndex] = v->oldv()->newv()[hexNum-1]->id();
									hvs[vIndex + 4] = v->oldv()->id();
								}
								else
								{
									hvs[vIndex] = v->newv()[hexNum-1]->id();
									hvs[vIndex + 4] = v->newv()[0]->id();
								}
							}
							else
							{
								if (v->newv().size() == 0)
								{
									hvs[vIndex] = v->oldv()->newv()[hexNum - 1]->id();
									hvs[vIndex + 4] = v->oldv()->newv()[hexNum]->id();;
								}
								else
								{
									hvs[vIndex] = v->newv()[hexNum - 1]->id();
									hvs[vIndex + 4] = v->newv()[hexNum]->id();
								}
							}
						}
						H* newh1 = new H();
						mesh->maxHid()++;
						mesh->m_map_hexs.insert(std::pair<int, H*>(mesh->maxHid(), newh1));
						mesh->hs.push_back(newh1);
						mesh->construct_hex(newh1, mesh->maxHid(), hvs);
						newhs.push_back(newh1);
					}
				}
				
				traceE = mesh->e_parallel_e_in_face(traceF, traceE);

			}

		}

		//�޸�����
		for (int hIndex = 0; hIndex < newhs.size(); hIndex++)
		{
			H* h = newhs[hIndex];
			mesh->reviese_hex_face_normal(h);
		}

		clock_t end = clock();
		std::cout << "singulairty insert time : " << end - start << std::endl;

	}

	template<typename V, typename E, typename F, typename H, typename M>
	int topoOperation<V, E, F, H, M>::num_of_layers(M* mesh, section* s)
	{
		int num_of_layer = 1;
		for (int fixEIndex = 0; fixEIndex < s->fixed_edges().size(); fixEIndex++)
		{		
			//�ҵ��̶�������һ��ıߣ�����hex
			E* e = s->fixed_edges()[fixEIndex];
			
			if (s->layers_vertieces().size()<=0)
			{
				std::unordered_set<V*> temp_v_set;
				s->layers_vertieces().push_back(temp_v_set);
			}
			for (int evIndex = 0; evIndex < e->vs.size(); evIndex++)
			{
				V* v = mesh->idVertices(e->vs[evIndex]);
				s->layers_vertieces()[0].insert(v);
			}

			F* f = NULL;
			for (int efIndex = 0; efIndex < e->neighbor_fs.size(); efIndex++)
			{
				F* adjf = mesh->idFaces(e->neighbor_fs[efIndex]);
				if (adjf->singularity_insert())
				{
					f = adjf;
				}
			}
			E* traceE = mesh->e_parallel_e_in_face(f, e);

			if (s->layers_vertieces().size() <= 1)
			{
				std::unordered_set<V*> temp_v_set;
				s->layers_vertieces().push_back(temp_v_set);
			}
			for (int evIndex = 0; evIndex < traceE->vs.size(); evIndex++)
			{
				V* v = mesh->idVertices(traceE->vs[evIndex]);
				s->layers_vertieces()[1].insert(v);
			}
			F* traceF = f;
			
			int temp_layer_num = 2;
			while (true)
			{
				std::vector<F*> adj_insert_f;
				for (int efIndex = 0; efIndex < traceE->neighbor_fs.size(); efIndex++)
				{
					F* adjf = mesh->idFaces(traceE->neighbor_fs[efIndex]);
					if (adjf->singularity_insert())
					{
						adj_insert_f.push_back(adjf);
					}
				}
				if (adj_insert_f.size() == 1)
				{
					break;
				}
				traceF = adj_insert_f[0]->id() == traceF->id() ? adj_insert_f[1] : adj_insert_f[0];
				traceE = mesh->e_parallel_e_in_face(traceF, traceE);

				if (s->layers_vertieces().size() <= temp_layer_num)
				{
					std::unordered_set<V*> temp_v_set;
					s->layers_vertieces().push_back(temp_v_set);
				}
				for (int evIndex = 0; evIndex < traceE->vs.size(); evIndex++)
				{
					V* v = mesh->idVertices(traceE->vs[evIndex]);
					s->layers_vertieces()[temp_layer_num].insert(v);
				}
				temp_layer_num++;
			}
			num_of_layer = num_of_layer > temp_layer_num ? num_of_layer : temp_layer_num;

		}
		return num_of_layer;
	}

	template<typename V, typename E, typename F, typename H, typename M>
	int topoOperation<V, E, F, H, M>::both_end_num_of_layers(M* mesh, both_end_section* s)
	{
		int num_of_layer = 1;
		for (int fixEIndex = 0; fixEIndex < s->fixed_edges1().size(); fixEIndex++)
		{
			//�ҵ��̶�������һ��ıߣ�����hex
			E* e = s->fixed_edges1()[fixEIndex];

			if (s->layers_vertieces().size() <= 0)
			{
				std::unordered_set<V*> temp_v_set;
				s->layers_vertieces().push_back(temp_v_set);
			}
			for (int evIndex = 0; evIndex < e->vs.size(); evIndex++)
			{
				V* v = mesh->idVertices(e->vs[evIndex]);
				s->layers_vertieces()[0].insert(v);
			}

			F* f = NULL;
			for (int efIndex = 0; efIndex < e->neighbor_fs.size(); efIndex++)
			{
				F* adjf = mesh->idFaces(e->neighbor_fs[efIndex]);
				if (adjf->singularity_insert())
				{
					f = adjf;
				}
			}
			E* traceE = mesh->e_parallel_e_in_face(f, e);

			if (s->layers_vertieces().size() <= 1)
			{
				std::unordered_set<V*> temp_v_set;
				s->layers_vertieces().push_back(temp_v_set);
			}
			for (int evIndex = 0; evIndex < traceE->vs.size(); evIndex++)
			{
				V* v = mesh->idVertices(traceE->vs[evIndex]);
				s->layers_vertieces()[1].insert(v);
			}
			F* traceF = f;

			int temp_layer_num = 2;
			while (true)
			{
				std::vector<F*> adj_insert_f;
				for (int efIndex = 0; efIndex < traceE->neighbor_fs.size(); efIndex++)
				{
					F* adjf = mesh->idFaces(traceE->neighbor_fs[efIndex]);
					if (adjf->singularity_insert())
					{
						adj_insert_f.push_back(adjf);
					}
				}
				std::vector<E*>::iterator eite = find(s->fixed_edges2().begin(), s->fixed_edges2().end(), traceE);
				if (eite!=s->fixed_edges2().end())
				{
					break;
				}
				traceF = adj_insert_f[0]->id() == traceF->id() ? adj_insert_f[1] : adj_insert_f[0];
				traceE = mesh->e_parallel_e_in_face(traceF, traceE);

				if (s->layers_vertieces().size() <= temp_layer_num)
				{
					std::unordered_set<V*> temp_v_set;
					s->layers_vertieces().push_back(temp_v_set);
				}
				for (int evIndex = 0; evIndex < traceE->vs.size(); evIndex++)
				{
					V* v = mesh->idVertices(traceE->vs[evIndex]);
					s->layers_vertieces()[temp_layer_num].insert(v);
				}
				temp_layer_num++;
			}
			num_of_layer = num_of_layer > temp_layer_num ? num_of_layer : temp_layer_num;
		}
		return num_of_layer;
	}

	template<typename V, typename E, typename F, typename H, typename M>
	void topoOperation<V, E, F, H, M>:: get_sheet(M* mesh, E* e,sheet* sh)
	{
		if (e->sheetId())
		{
			return;
		}
		std::queue<E*> eQueue;
		int sheetId = sh->id();
		eQueue.push(e);
		e->saved() = true;
		while (!eQueue.empty())
		{
			E* pe = eQueue.front();
			eQueue.pop();
			pe->sheetId() = sheetId;
			
			sh->sheets.push_back(pe);
			//save the e and hex in sh, and push element in eQueue;
			for (int ehIndex = 0; ehIndex < pe->neighbor_hs.size(); ehIndex++)
			{
				H* nh = mesh->idHexs(pe->neighbor_hs[ehIndex]);
				int heIndex = nh->edgeIndex(pe->id());
				bool is_continue = true;
				for (int peIndex = 0; peIndex < 3; peIndex++)
				{
					E* parallel_e = mesh->idEdges(nh->es[nh->eparalle[heIndex][peIndex]]);
					if (parallel_e->sheetId()!=sh->id())
					{
						is_continue = false;
						break;
					}
				}
				if (is_continue) continue;
				if (nh->sheetId()!=sheetId)
				{
					nh->sheetId() = sheetId;
					sh->hs.push_back(nh);
				}
				
				heIndex = nh->edgeIndex(pe->id());
				for (int peIndex = 0; peIndex < 3; peIndex++)
				{
					E* parallel_e = mesh->idEdges(nh->es[nh->eparalle[heIndex][peIndex]]);
					if ((parallel_e->sheetId()!= sheetId)&&(!parallel_e->saved()))
					{
						eQueue.push(parallel_e);
						parallel_e->saved() = true;
					}
				}

			}
		}

	}

	template<typename V, typename E, typename F, typename H, typename M>
	void topoOperation<V, E, F, H, M>::mark_sheet_information(M* mesh, sheet* sh)
	{

	}

	template<typename V, typename E, typename F, typename H, typename M>
	void topoOperation<V, E, F, H, M>::sheet_collapse(M* mesh, sheet* s)
	{
		//pretreatment: find the self-corss hex and normal hex
		std::vector<H*> self_cross_hs;
		std::vector<H*> normal_hs;
		for (int hIndex = 0; hIndex < s->hs.size(); hIndex++)
		{
			TH* h = s->hs[hIndex];
			TF* f = mesh->idFaces(h->fs[0]);
			TE* e = mesh->idEdges(f->es[0]);
			TV* v = mesh->idVertices(e->vs[0]);
			std::vector<TE*> es;
			es.push_back(e);
			e = mesh->flip_e(f, e, v);
			es.push_back(e);
			f = mesh->flip_f(h, f, e);
			e = mesh->flip_e(f, e, v);
			es.push_back(e);
			int num = 0;
			for (int eIndex = 0; eIndex < es.size(); eIndex++)
			{
				TE* pe = es[eIndex];
				if (pe->sheetId() == s->id())
				{
					num++;
				}
			}
			if (num==2)
			{
				self_cross_hs.push_back(h);

			}
			else
			{
				normal_hs.push_back(h);
			}
		}
		std::cout <<"self_cross_hs: " << self_cross_hs.size() << std::endl;
		std::cout <<"normal_hs: " << normal_hs.size() << std::endl;

		//mark sheet id in sheet hex
		std::unordered_set<V*> oldVs;// need to delete vertex
		std::unordered_set<E*> oldEs;//need to delete edge
		std::unordered_set<F*> oldFs;//need to delte face
		for (int hIndex = 0; hIndex < s->hs.size(); hIndex++)
		{
			H* h = s->hs[hIndex];
			h->sheetId() = s->id();
			for (int hvIndex = 0; hvIndex < h->vs.size(); hvIndex++)
			{
				V* hv = mesh->idVertices(h->vs[hvIndex]);
				oldVs.insert(hv);
				//hv->is_delete() = true;
			}
			for (int heIndex = 0; heIndex < h->es.size(); heIndex++)
			{
				E* he = mesh->idEdges(h->es[heIndex]);
				oldEs.insert(he);
				he->is_delete() = true;
			}
			for (int hfIndex = 0; hfIndex < h->fs.size(); hfIndex++)
			{
				F* hf = mesh->idFaces(h->fs[hfIndex]);
				oldFs.insert(hf);
				hf->is_delete() = true;
			}
		}
		
		//construct new elements in cross_self hs
		for (int hIndex = 0; hIndex < self_cross_hs.size(); hIndex++)
		{
			H* h = self_cross_hs[hIndex];
			std::vector<E*> parallel_e;
			//find the edges not belong to sheet;
			E* e = NULL;
			for (int heIndex = 0; heIndex < h->es.size(); heIndex++)
			{
				E* pe = mesh->idEdges(h->es[heIndex]);
				if (pe->sheetId()!=s->id())
				{
					e = pe;
					break;
				}
			}
			parallel_e = mesh->e_parallel_e_in_hex(h, e);
			parallel_e.push_back(e);
			//find the parallel faces belong to sheet
			F* f = NULL;
			for (int efIndex = 0; efIndex < e->neighbor_fs.size(); efIndex++)
			{
				F* pf = mesh->idFaces(e->neighbor_fs[efIndex]);
				if (h->faceIndex(pf->id())!=-1)
				{
					f = pf;
					break;
				}
			}
			V* v = mesh->idVertices(e->vs[0]);
			e = mesh->flip_e(f, e, v);
			f = mesh->flip_f(h, f, e);
			std::vector<F*> parallel_f;//faces in sheets
			parallel_f.push_back(f);
			parallel_f.push_back(mesh->f_parallel_f_in_hex(h, f));
			
			//construct new edges
			E* newe = new_edge(mesh);
			//get two new verticces
			for (int fIndex = 0; fIndex < parallel_f.size(); fIndex++)
			{
				F* pf = parallel_f[fIndex];
				V* newv = NULL;
				if (pf->newv!=NULL)
				{
					newv = pf->newv;
					newe->vs.push_back(newv->id());
					newv->neighbor_es.push_back(newe->id());
					continue;
				}			
				newv = new_vertex(mesh);
				//construct new vetex and neighbor information
				for (int fvIndex = 0; fvIndex < f->vs.size(); fvIndex++)
				{
					
					V* v = mesh->idVertices(f->vs[fvIndex]);
					//assign the position information
					newv->position() += v->position() * 0.4;
					//assign the new vertex neighbor information
					//assign the neighbor edge of new vertex
					for (int veIndex = 0; veIndex < v->neighbor_es.size(); veIndex++)
					{
						E* ve = mesh->idEdges(v->neighbor_es[veIndex]);
						if (!ve->is_delete())
						{
							newv->neighbor_es.push_back(ve->id());
							int vevIndex = ve->vertexIndex(v->id());
							ve->vs[vevIndex] = newv->id();
						}
					}
					//assign the neighbor face of new vertex
					for (int vfIndex = 0; vfIndex < v->neighbor_fs.size(); vfIndex++)
					{
						F* vf = mesh->idFaces(v->neighbor_fs[vfIndex]);
						/*bool is_neighbor_information = true;
						for (int feIndex = 0; feIndex < vf->es.size(); feIndex++)
						{
							E* fe = mesh->idEdges(vf->es[feIndex]);
							if (fe->sheetId() == s->id());
							{
								is_neighbor_information = false;
								break;
							}
						}*/
						if (!vf->is_delete())
						{
							newv->neighbor_fs.push_back(vf->id());
							int vfvIndex = vf->vertexIndex(v->id());
							vf->vs[vfvIndex] = newv->id();
						}
					}
					//assign the neighbor hex of new vertex
					for (int vhIndex = 0; vhIndex < v->neighbor_hs.size(); vhIndex++)
					{
						H* vh = mesh->idHexs(v->neighbor_hs[vhIndex]);
						if (vh->sheetId()!=s->id())
						{
							newv->neighbor_hs.push_back(vh->id());
							int vhvIndex = vh->vertexIndex(v->id());
							vh->vs[vhvIndex] = newv->id();
						}
					}

				}
				//assign the new vertex to old edges
				for (int feIndex = 0; feIndex < f->es.size(); feIndex++)
				{
					E* e = mesh->idEdges(f->es[feIndex]);
					e->newv = newv;
				}
				newe->vs.push_back(newv->id());
				newv->neighbor_es.push_back(newe->id());
			}
			h->newe = newe;
			//assign the neighbor information of new edge		
			for (int peIndex = 0; peIndex < parallel_e.size(); peIndex++)
			{
				E* pe = parallel_e[peIndex];
				//assign the neighbor face of new edge
				for (int pefIndex = 0; pefIndex < pe->neighbor_fs.size(); pefIndex++)
				{
					F* pef = mesh->idFaces(pe->neighbor_fs[pefIndex]);
					/*bool is_neighbor_f = true;
					for (int pefeIndex = 0; pefeIndex < pef->es.size(); pefeIndex++)
					{
						E* pefe = mesh->idEdges(pef->es[pefeIndex]);
						if (pefe->sheetId()==s->id())
						{
							is_neighbor_f = false;
						}
					}*/
					if (!pef->is_delete())
					{
						newe->neighbor_fs.push_back(pef->id());
						int pefeIndex = pef->edgeIndex(pe->id());
						pef->es[pefeIndex] = newe->id();
					}
				}
				//assign the neighbor hex of new edge
				for (int pehIndex = 0; pehIndex < pe->neighbor_hs.size(); pehIndex++)
				{
					H* peh = mesh->idHexs(pe->neighbor_hs[pehIndex]);
					if (peh->sheetId()!=s->id())
					{
						newe->neighbor_hs.push_back(peh->id());
						int peheIndex = peh->edgeIndex(pe->id());
						peh->es[peheIndex] = newe->id();
					}
				}

			}


		}
		
		//construct new elements in normal hs
		for (int hIndex = 0; hIndex < normal_hs.size(); hIndex++)
		{		
			H* h =normal_hs[hIndex];
			F* f = NULL;//the face without sheet edge		

			//find the f without sheet edge
			for (int hfIndex = 0; hfIndex < h->fs.size(); hfIndex++)
			{
				F* pf = mesh->idFaces(h->fs[hfIndex]);
				int num = 0;
				for (int feIndex = 0; feIndex < pf->es.size(); feIndex++)
				{
					E* pe = mesh->idEdges(pf->es[feIndex]);
					if (pe->sheetId()==s->id())
					{
						num++;			
					}
				}
				if (num==0)
				{
					f = pf;
					break;
				}
			}
			std::vector<F*> newf_parallel_f;//for create new faces, find the parallel old faces don't belong to the sheets, for assign the neighbor information of new faces
			newf_parallel_f.push_back(f);
			newf_parallel_f.push_back(mesh->f_parallel_f_in_hex(h,f));

			//find the four sides in order and new edge parallel edges
			std::vector<V*> newvs;
			std::vector<E*> sheet_es;//the edges in sheets in order
			std::vector<std::vector<E*>> newe_parallel_e;//for create new edges, find the new edges parallel the old edges don't belong to sheets, for assign the neighbor information of new edges
			for (int fvIndex = 0; fvIndex < f->vs.size(); fvIndex++)
			{
				V* v = mesh->idVertices(f->vs[fvIndex]);
				int nextvId = f->vs[(fvIndex + 1) % (f->vs.size())];
				std::vector<E*>parallel_e;// the parallel edges in two face, which don't contain sheets
				for (TMesh::VEInHexIterator veite(mesh, v, h); !veite.end(); veite++)
				{		
					E* e = *veite;
					if (e->sheetId() == s->id())
					{
						sheet_es.push_back(e);
						e->mark() = true;
						
					}	
					if (e->sheetId()!=s->id()&&e->vertexIndex(nextvId) != -1)
					{
						parallel_e.push_back(e);
						F* temp_f = mesh->flip_f(h, f, e);
						parallel_e.push_back(mesh->e_parallel_e_in_face(temp_f, e));

					}
				}
				if (parallel_e.size() != 2)
				{
					std::cout << "----------------------------" << std::endl;
					std::cout << "û�ҵ�ƽ�б�" << std::endl;
					std::cout << "parallel_e.size()��"<< parallel_e.size() << std::endl;
					std::cout << "v id " << v->id() << std::endl;
					v->mark() = true;
					std::cout << "next v id " << nextvId << std::endl;
					mesh->idVertices(nextvId)->mark() = true;
					std::cout << "f id " << f->id() << std::endl;
					std::cout << "h id" << h->id() << std::endl;
					for (TMesh::VEInHexIterator veite(mesh, v, h); !veite.end(); veite++)
					{
						E* e = *veite;
						std::cout << "ve->sheetId() " << e->sheetId() << std::endl;
						std::cout << "ve->vs: " << e->vs[0] << " " << e->vs[1] << std::endl;
						mesh->idVertices(e->vs[0])->mark() = true;
						mesh->idVertices(e->vs[1])->mark() = true;
						std::cout << "ve->vertexIndex(nextvId) " << e->vertexIndex(nextvId) << std::endl;
						
					}
					std::cout << "----------------------------" << std::endl;
					return;
				}			
				newe_parallel_e.push_back(parallel_e);

			}

			/*std::cout << " h id :" << h->id() << std::endl;
			std::cout << "newe_parallel_e.size(): "<< newe_parallel_e.size() << std::endl;
			std::cout << " oldes.size() :" << oldes.size() << std::endl;*/

			std::cout << "start new vertices" << std::endl;
			//construct new vertices
			for (int oldeIndex = 0; oldeIndex < 4; oldeIndex++)
			{
				E* e = sheet_es[oldeIndex];
				V* newv = NULL;
				if (e->newv!=NULL)
				{
					newv = e->newv;
					newvs.push_back(newv);
				}
				else
				{
					V* newv = new_vertex(mesh);
					for (int evIndex = 0; evIndex < e->vs.size(); evIndex++)
					{
						V* ev = mesh->idVertices(e->vs[evIndex]);
						//assign the neighbor edge of new vertex
						for (int eveIndex = 0; eveIndex < ev->neighbor_es.size(); eveIndex++)
						{
							E* eve = mesh->idEdges(ev->neighbor_es[eveIndex]);
							if (!eve->is_delete())
							{
								std::cout << "eve->id: " << eve->id() << std::endl;
								std::cout << "eve vs: " << eve->vs[0] << " " << eve->vs[1] << std::endl;
								std::cout << "ev->id(): " << ev->id() << std::endl;
								std::cout << "-------------------------------" << std::endl;
								newv->neighbor_es.push_back(eve->id());
								int evevIndex = eve->vertexIndex(ev->id());
								eve->vs[eveIndex] = newv->id();
							}
						}
						//assign the neighbor face of new vertex
						//for (int evfIndex = 0; evfIndex < ev->neighbor_fs.size(); evfIndex++)
						//{
						//	F* evf = mesh->idFaces(ev->neighbor_fs[evfIndex]);
						//	/*bool is_neighbor_f = true;
						//	for (int evfeIndex = 0; evfeIndex < evf->es.size(); evfeIndex++)
						//	{
						//		E* evfe = mesh->idEdges(evf->es[evfeIndex]);
						//		if (evfe->sheetId()==s->id())
						//		{
						//			is_neighbor_f = false;
						//			break;
						//		}
						//	}*/
						//	if (!evf->is_delete())
						//	{
						//		newv->neighbor_fs.push_back(evf->id());
						//		int evfvIndex = evf->vertexIndex(ev->id());
						//		evf->vs[evfvIndex] = newv->id();

						//	}
						//}
						////assign the neighbor hex of new vertex
						//for (int evhIndex = 0; evhIndex < ev->neighbor_hs.size(); evhIndex++)
						//{
						//	H* evh = mesh->idHexs(ev->neighbor_hs[evhIndex]);
						//	if (evh->sheetId()!=s->id())
						//	{
						//		newv->neighbor_hs.push_back(evh->id());
						//		int evhvIndex = evh->vertexIndex(ev->id());
						//		evh->vs[evhvIndex] = newv->id();
						//	}
						//}

					}
					newv->position() = (mesh->idVertices(e->vs[0])->position() + mesh->idVertices(e->vs[1])->position()) * 0.5;
					newvs.push_back(newv);
					e->newv = newv;
				}
				//break;
			}
			
			//std::cout << "start new edges" << std::endl;
			////construct new edges 
			//std::vector<E*> newes;
			//for (int vIndex = 0; vIndex < newvs.size(); vIndex++)
			//{
			//	std::vector<V*> evs = { newvs[vIndex],newvs[(vIndex + 1) % newvs.size()] };
			//	std::vector<int> evsId = {evs[0]->id(), evs[1]->id()};
			//	V* v = evs[0];
			//	bool construct_edge = true;
			//	for (int veIndex = 0; veIndex < v->neighbor_es.size(); veIndex++)
			//	{
			//		E* ne = mesh->idEdges(v->neighbor_es[veIndex]);
			//		if (ne->is_v_equal(evsId))
			//		{
			//			newes.push_back(ne);
			//			construct_edge = false;
			//			break;
			//		}
			//	}
			//	//construct edge
			//	if (construct_edge)
			//	{
			//		std::cout <<"construct_edge" <<construct_edge << std::endl;
			//		E* newe = new_edge(mesh);
			//		newes.push_back(newe);
			//		newe->vs = evsId;
			//		for (int evIndex = 0; evIndex < evs.size(); evIndex++)
			//		{
			//			V* ev = evs[evIndex];
			//			ev->neighbor_es.push_back(newe->id());
			//		}
			//		
			//		//assign the neighbor information of hex
			//		//std::cout << "parallel_e.size(): " << newe_parallel_e[vIndex].size() << std::endl;
			//		std::vector<E*> parallel_e = newe_parallel_e[vIndex];
			//		
			//		for (int peIndex = 0; peIndex < parallel_e.size(); peIndex++)
			//		{
			//			E* pe = parallel_e[peIndex];
			//			//assign the neighbor face of new edge
			//			for (int pefIndex = 0; pefIndex < pe->neighbor_fs.size(); pefIndex++)
			//			{
			//				F* pef = mesh->idFaces(pe->neighbor_fs[pefIndex]);
			//				/*bool is_neighbor_f = true;
			//				for (int pefeIndex = 0; pefeIndex < pef->es.size(); pefeIndex++)
			//				{
			//					E* pefe = mesh->idEdges(pef->es[pefeIndex]);
			//					if (pefe->sheetId()==s->id())
			//					{
			//						is_neighbor_f = false;
			//						break;
			//					}
			//				}*/
			//				if (!pef->is_delete())
			//				{
			//					newe->neighbor_fs.push_back(pef->id());
			//					int pefeIndex = pef->edgeIndex(pe->id());
			//					pef->es[pefeIndex] = newe->id();
			//				}
			//			}
			//			
			//			//assign the neighbor hex of new edge
			//			for (int pehIndex = 0; pehIndex < pe->neighbor_hs.size(); pehIndex++)
			//			{
			//				H* peh = mesh->idHexs(pe->neighbor_hs[pehIndex]);
			//				if (peh->sheetId()!=s->id())
			//				{
			//					newe->neighbor_hs.push_back(peh->id());
			//					int peheIndex = peh->edgeIndex(pe->id());
			//					peh->es[peheIndex] = newe->id();
			//				}
			//			}
			//		}
			//	}
			//}

			//construct new faces
			/*std::cout << "start new faces" << std::endl;
			F* newf = new_face(mesh);
			for (int fvIndex = 0; fvIndex < newvs.size(); fvIndex++)
			{
				V* fv = newvs[fvIndex];
				newf->vs.push_back(fv->id());
				fv->neighbor_fs.push_back(newf->id());
			}
			for (int feIndex = 0; feIndex < newes.size(); feIndex++)
			{
				E* fe = newes[feIndex];
				newf->es.push_back(fe->id());
				fe->neighbor_fs.push_back(newf->id());
			}*/
			

			//assign the neighbor information for new face
			/*for (int pfIndex = 0; pfIndex < newf_parallel_f.size(); pfIndex++)
			{
				F* pf = newf_parallel_f[pfIndex];
				for (int pfhIndex = 0; pfhIndex < pf->neighbor_hs.size(); pfhIndex++)
				{
					H* pfh = mesh->idHexs(pf->neighbor_hs[pfhIndex]);
					if (pfh==NULL)
					{
						std::cout << "pf->id():" << pf->id() << std::endl;
						std::cout << "pf->neighbor_hs.size(): " << pf->neighbor_hs.size() << std::endl;
						std::cout << "pfhIndex:" << pfhIndex << std::endl;
						std::cout << "pf->neighbor_hs[pfhIndex]: " << pf->neighbor_hs[pfhIndex] << std::endl;
						std::cout << "h->id(): " << h->id() << std::endl;
						
					}
					if (pfh->sheetId()!=s->id())
					{
						newf->neighbor_hs.push_back(pfh->id());
						int pfIndex = pfh->faceIndex(pf->id());
						if (pfIndex==-1)
						{
							std::cout << "pfIndex = -1" << std::endl;
						}
						pfh->fs[pfIndex] = newf->id();
					}
				}
			}*/
		}

		//std::cout << "start delete elements" << std::endl;
		////delete elements
		//for (int hIndex = 0; hIndex < self_cross_hs.size(); hIndex++)
		//{
		//	H* h = self_cross_hs[hIndex];
		//	
		//	mesh->delete_hex(h);
		//}
		//for (std::unordered_set<F*>::iterator fite = oldFs.begin();fite!=oldFs.end();fite++)
		//{
		//	F* oldf = *fite;
		//	mesh->delete_face(oldf);
		//}
		//for (std::unordered_set<E*>::iterator eite = oldEs.begin();eite!=oldEs.end();eite++)
		//{
		//	E* olde = *eite;
		//	mesh->delete_edge(olde);
		//}
		//for (std::unordered_set<V*>::iterator vite = oldVs.begin(); vite != oldVs.end(); vite++)
		//{
		//	V* oldv = *vite;
		//	mesh->delete_vertex(oldv);
		//}
		std::cout << "---------------------------" << std::endl;

	}



	template<typename V, typename E, typename F, typename H, typename M>
	V* topoOperation<V, E, F, H, M>::new_vertex(M* mesh)
	{
		V* v = new V;
		mesh->maxVid()++;
		v->id() = mesh->maxVid();
		mesh->vs.push_back(v);
		mesh->m_map_vertices.insert(std::pair<int, V*>(v->id(), v));
		return v;
	}
	template<typename V, typename E, typename F, typename H, typename M>
	E* topoOperation<V, E, F, H, M>::new_edge(M* mesh)
	{
		E* e = new E;
		mesh->maxEid()++;
		e->id() = mesh->maxEid();
		mesh->es.push_back(e);
		mesh->m_map_edges.insert(std::pair<int, E*>(e->id(), e));
		return e;
	}
	template<typename V, typename E, typename F, typename H, typename M>
	F* topoOperation<V, E, F, H, M>::new_face(M* mesh)
	{
		F* f = new F;
		mesh->maxFid()++;
		f->id() = mesh->maxFid();
		mesh->fs.push_back(f);
		mesh->m_map_faces.insert(std::pair<int, F*>(f->id(), f));
		return f;
	}
	template<typename V, typename E, typename F, typename H, typename M>
	H* topoOperation<V, E, F, H, M>::new_hex(M* mesh)
	{
		H* h = new H;
		mesh->maxHid()++;
		h->id() = mesh->maxFid();
		mesh->hs.push_back(h);
		mesh->m_map_hexs.insert(std::pair<int, H*>(h->id(), h));
		return h;
	}

}


#endif // !TOPOLOGICAL_OPERATION_H
