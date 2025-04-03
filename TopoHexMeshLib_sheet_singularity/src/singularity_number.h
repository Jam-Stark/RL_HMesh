#pragma once
#include"topoMesh.h"
#include<queue>

namespace HMeshLib
{
	template<typename M>
	class get_singularity_number
	{
		typedef typename M::V V;
		typedef typename M::E E;
		typedef typename M::F F;
		typedef typename M::H H;
	public:
		get_singularity_number(M* input_mesh);
		~get_singularity_number();
		void get_connect_singualarity_edge();
		singularity* trace_singularity(E* e, int singularity_id);
		void generate_singularity_number(M* input_mesh);

		/*param*/
		int singualarity_id;
		std::vector<singularity*> singularities;
	private:
		M* mesh;
	};
	template<typename M>
	get_singularity_number<M>::get_singularity_number(M* input_mesh)
	{
		mesh = input_mesh;
		singualarity_id = 0;
		get_connect_singualarity_edge();
		
	}
	template<typename M>
	get_singularity_number<M>::~get_singularity_number()
	{
	}

	template<typename M>
	void get_singularity_number<M>::get_connect_singualarity_edge()
	{
		for (M::MEIterator eite(mesh); !eite.end(); eite++)
		{
			E* e = *eite;
			if (e->singularity() && !e->singularity_id())
			{
				singualarity_id++;
				singularity* s = trace_singularity(e, singualarity_id);
			}
		}

	}

	template<typename M>
	singularity* get_singularity_number<M>::trace_singularity(E* e, int singularity_id)
	{
		if (e->singularity_id())
		{
			return NULL;
		}
		singularity* s = new singularity();
		singularities.push_back(s);
		s->degree() = e->neighbor_hs.size();
		s->id() = singualarity_id;

		std::queue<E*> eQueue;
		eQueue.push(e);
		while (!eQueue.empty())
		{
			E* pE = eQueue.front();
			eQueue.pop();
			pE->singularity_id() = singularity_id;
			s->es.push_back(pE);
			int v_saved_count = 0;
			for (int evIndex = 0; evIndex < pE->vs.size(); evIndex++)
			{
				V* v = mesh->idVertices(pE->vs[evIndex]);
				if (v->saved() == singularity_id)
				{
					v_saved_count++;
					if (v_saved_count == 2)
					{
						s->loop() = true;
					}
					continue;
				}
				v->saved() = singualarity_id;
				s->vs.push_back(v);

				int num = 0;
				E* nextE = NULL;
				for (int veIndex = 0; veIndex < v->neighbor_es.size(); veIndex++)
				{
					E* adjE = mesh->idEdges(v->neighbor_es[veIndex]);
					if (adjE->id() != pE->id() && adjE->singularity())
					{
						num++;
						nextE = adjE;
					}
				}
				if (num == 1 && nextE->neighbor_hs.size() == pE->neighbor_hs.size() && nextE->boundary() == pE->boundary())
				{
					eQueue.push(nextE);
					nextE->saved() = true;

				}
			}
		}
		return s;
	}

	template<typename M>
	void get_singularity_number<M>::generate_singularity_number(M* input_mesh)
	{
		mesh = input_mesh;
		singualarity_id = 0;
		get_connect_singualarity_edge();
	}


}