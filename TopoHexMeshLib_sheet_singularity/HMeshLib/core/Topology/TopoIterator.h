/*iterator*/
#ifndef TOPO_ITERATOR_H
#define TOPO_ITERATOR_H

#include"TopoE.h"
#include"TopoV.h"
#include"TopoF.h"
#include"TopoM.h"
#include"TopoH.h"
#include<string>
#include<list>
namespace HMeshLib
{
	/*
	* go though all vertices
	*/
	template<typename V, typename E, typename F, typename H>
	class TopoMeshVertexIterator
	{
	public:
		TopoMeshVertexIterator(topoM< V,  E,  F,  H>* InputtopoMesh)
		{
			topoMesh = InputtopoMesh;
			m_iter = topoMesh->vs.begin();
		}
		~TopoMeshVertexIterator() {};
		V* value() { return *m_iter; };
		void operator++() { m_iter++; };
		void operator++(int) { m_iter ++; };
		bool end() { return m_iter == topoMesh->vs.end(); };
		V* operator*() { return value(); };
	protected:
		topoM< V,  E,  F,  H>* topoMesh;
		typename std::list<V*>::iterator m_iter;
	};

	/*
	* go though all edges
	*/
	template<typename V, typename E, typename F, typename H>
	class TopoMeshEdgeIterator
	{
	public:
		TopoMeshEdgeIterator(topoM< V,  E,  F,  H>* InputtopoMesh)
		{
			topoMesh = InputtopoMesh;
			m_iter = topoMesh->es.begin();
		}
		~TopoMeshEdgeIterator() {};
		E* value() { return *m_iter; };
		void operator++() { m_iter++; };
		void operator++(int) { m_iter++; };
		bool end() { return m_iter == topoMesh->es.end(); };
		E* operator*() { return value(); };
	protected:
		topoM< V,  E,  F,  H>* topoMesh;
		typename std::list<E*>::iterator m_iter;
	};

	/*
	* go though all faces
	*/
	template<typename V, typename E, typename F, typename H>
	class TopoMeshFaceIterator
	{
	public:
		TopoMeshFaceIterator(topoM< V,  E,  F,  H>* InputtopoMesh)
		{
			topoMesh = InputtopoMesh;
			m_iter = topoMesh->fs.begin();
		}
		~TopoMeshFaceIterator() {};
		F* value() { return *m_iter; };
		void operator++() { m_iter ++; };
		void operator++(int) { m_iter ++; };
		bool end() { return m_iter == topoMesh->fs.end(); };
		F* operator*() { return value(); };

	protected:
		topoM< V,  E,  F,  H>* topoMesh;
		typename std::list<F*>::iterator m_iter;
	};
	
	/*	
	* go though all hexs
	*/
	template<typename V, typename E, typename F, typename H>
	class TopoMeshHexIterator
	{
	public:
		TopoMeshHexIterator(topoM< V,  E,  F,  H>* InputtopoMesh)
		{
			topoMesh = InputtopoMesh;
			m_iter = topoMesh->hs.begin();
		}
		~TopoMeshHexIterator() {};
		H* value() { return *m_iter; };
		void operator++() { m_iter++; };
		void operator++(int) { m_iter++; };
		bool end() { return m_iter == topoMesh->hs.end(); };
		H* operator*() { return value(); };
	protected:
		topoM< V,  E,  F,  H>* topoMesh;
		typename std::list<H*>::iterator m_iter;
	};

	/*	
	* go through all the adjcent vertex of a vertex
	*/
	template<typename V, typename E, typename F, typename H>
	class TopoVertexVertexIterator 
	{
	public:
		TopoVertexVertexIterator(topoM< V,  E,  F,  H>* inputM,V* inputV)
		{
			m = inputM;
			v = inputV;
			m_iter = v->neighbor_es.begin();
		}
		~TopoVertexVertexIterator() {};
		V* value() {
			int adjeId = *m_iter;
			E* adje = m->idEdges(adjeId);
			V* adjv = m->flip_v(adje, v);
			return adjv;
		}
		void operator++() { m_iter++; };
		void operator++(int) { m_iter++; };
		bool end() { return m_iter == v->neighbor_es.end(); };
		V* operator*() { return value(); };
	protected:
		V* v;
		topoM< V, E, F, H>* m;
		std::vector<int>::iterator m_iter;
	};

	/*go through all the adjcent edges of a vertex*/
	template<typename V, typename E, typename F, typename H>
	class TopoVertexEdgeIterator
	{
	public:
		TopoVertexEdgeIterator(topoM< V, E, F, H>* inputM, V* inputV)
		{
			m = inputM;
			v = inputV;
			m_iter = v->neighbor_es.begin();
		};
		~TopoVertexEdgeIterator() {};
		E* value() {
			int adjEId = *m_iter;
			return m->idEdges(adjEId);
		}
		void operator++() { m_iter++; };
		void operator++(int) { m_iter++; };
		bool end() { return m_iter == v->neighbor_es.end(); };
		E* operator*() { return value(); };

	protected:
		V* v;
		topoM< V, E, F, H>* m;
		std::vector<int>::iterator m_iter;
	};

	/*go through ajdcent face of a vertex*/
	template<typename V, typename E, typename F, typename H>
	class TopoVertexFaceIterator
	{
	public:
		TopoVertexFaceIterator(topoM< V, E, F, H>* inputM, V* inputV)
		{
			m = inputM;
			v = inputV;
			m_iter = v->neighbor_fs.begin();
		};
		~TopoVertexFaceIterator() {};
		F* value() {
			int adjFId = *m_iter;
			return m->idFaces(adjFId);
			
		}
		void operator++() { m_iter++; };
		void operator++(int) { m_iter++; };
		bool end() { return m_iter == v->neighbor_fs.end(); };
		F* operator*() { return value(); };

	protected:
		V* v;
		topoM< V, E, F, H>* m;
		std::vector<int>::iterator m_iter;
	};

	/*go through ajdcent hex of a vertex*/
	template<typename V, typename E, typename F, typename H>
	class TopoVertexHexIterator
	{
	public:
		TopoVertexHexIterator(topoM< V, E, F, H>* inputM, V* inputV)
		{
			m = inputM;
			v = inputV;
			m_iter = v->neighbor_hs.begin();
		};
		~TopoVertexHexIterator() {};
		H* value() {
			int adjHId = *m_iter;
			return m->idHexs(adjHId);
		}
		void operator++() { m_iter++; };
		void operator++(int) { m_iter++; };
		bool end() { return m_iter == v->neighbor_hs.end(); };
		H* operator*() { return value(); };

	protected:
		V* v;
		topoM< V, E, F, H>* m;
		std::vector<int>::iterator m_iter;
	};

	/*go through all adjcent face of a edge*/
	template<typename V, typename E, typename F, typename H>
	class TopoEdgeFaceIterator
	{
	public:
		TopoEdgeFaceIterator(topoM< V, E, F, H>* inputM, E* inputE)
		{
			m = inputM;
			e = inputE;
			m_iter = e->neighbor_fs.begin();
		};
		~TopoEdgeFaceIterator() {};
		F* value() {
			int adjFId = *m_iter;
			return m->idFaces(adjFId);
		}
		void operator++() { m_iter++; };
		void operator++(int) { m_iter++; };
		bool end() { return m_iter == e->neighbor_fs.end(); };
		F* operator*() { return value(); };

	protected:
		E* e;
		topoM< V, E, F, H>* m;
		std::vector<int>::iterator m_iter;
	};

	/*go through all adjcent hex of a edge*/
	template<typename V, typename E, typename F, typename H>
	class TopoEdgeHexIterator
	{
	public:
		TopoEdgeHexIterator(topoM< V, E, F, H>* inputM, E* inputE)
		{
			m = inputM;
			e = inputE;
			m_iter = e->neighbor_hs.begin();
		};
		~TopoEdgeHexIterator() {};
		H* value() {
			int adjHId = *m_iter;
			return m->idHexs(adjHId);
		}
		void operator++() { m_iter++; };
		void operator++(int) { m_iter++; };
		bool end() { return m_iter == e->neighbor_hs.end(); };
		H* operator*() { return value(); };

	protected:
		E* e;
		topoM< V, E, F, H>* m;
		std::vector<int>::iterator m_iter;
	};

	/*go through all adjcent hex of a edge*/
	template<typename V, typename E, typename F, typename H>
	class TopoFaceHexIterator
	{
	public:
		TopoFaceHexIterator(topoM< V, E, F, H>* inputM, F* inputF)
		{
			m = inputM;
			f = inputF;
			m_iter = f->neighbor_hs.begin();
		};
		~TopoFaceHexIterator() {};
		H* value() {
			int adjHId = *m_iter;
			return m->idHexs(adjHId);
		}
		void operator++() { m_iter++; };
		void operator++(int) { m_iter++; };
		bool end() { return m_iter == f->neighbor_hs.end(); };
		H* operator*() { return value(); };

	protected:
		F* f;
		topoM< V, E, F, H>* m;
		std::vector<int>::iterator m_iter;
	};
	
	/*go through all vertices of hex*/
	template<typename V, typename E, typename F, typename H>
	class TopoHexVertexIterator
	{
	public:
		TopoHexVertexIterator(topoM< V, E, F, H>* inputm, H* inputh)
		{
			h = inputh;
			m = inputm;
			m_iter = h->vs.begin();
		};
		~TopoHexVertexIterator() {};
		V* value() {
			int eId = *m_iter;
			return m->idVertices(eId);
		}
		void operator++() { m_iter++; };
		void operator++(int) { m_iter++; };
		bool end() { return m_iter == h->vs.end(); };
		V* operator*() { return value(); };

	protected:
		H* h;
		topoM< V, E, F, H>* m;
		std::vector<int>::iterator m_iter;
	};

	/*go through all edges of hex*/
	template<typename V, typename E, typename F, typename H>
	class TopoHexEdgeIterator
	{
	public:
		TopoHexEdgeIterator(topoM< V, E, F, H>* inputm,H* inputh)
		{
			h = inputh;
			m = inputm;
			m_iter = h->es.begin();
		};
		~TopoHexEdgeIterator() {};
		E* value() {
			int eId = *m_iter;
			return m->idEdges(eId);
		}
		void operator++() { m_iter++; };
		void operator++(int) { m_iter++; };
		bool end() { return m_iter == h->es.end(); };
		E* operator*() { return value(); };

	protected:
		H* h;
		topoM< V, E, F, H>* m;
		std::vector<int>::iterator m_iter;
	};

	/*go through all faces of hex*/
	template<typename V, typename E, typename F, typename H>
	class TopoHexFaceIterator
	{
	public:
		TopoHexFaceIterator(topoM< V, E, F, H>* inputm, H* inputh)
		{
			h = inputh;
			m = inputm;
			m_iter = h->fs.begin();
		};
		~TopoHexFaceIterator() {};
		F* value() {
			int fId = *m_iter;
			return m->idFaces(fId);
		}
		void operator++() { m_iter++; };
		void operator++(int) { m_iter++; };
		bool end() { return m_iter == h->fs.end(); };
		F* operator*() { return value(); };

	protected:
		H* h;
		topoM< V, E, F, H>* m;
		std::vector<int>::iterator m_iter;
	};

	/*go through all adjcent vertex of a vertex in a face*/
	template<typename V, typename E, typename F, typename H>
	class TopoVertexVertexInFaceIterator
	{
	public:
		TopoVertexVertexInFaceIterator(topoM< V, E, F, H>* inputM, F* inputF, V* inputV)
		{
			m = inputM;
			f = inputF;
			v = inputV;
			vIndex = f->vertexIndex(v->id());
			assert(vIndex != -1);
			m_iter = 0;
		}
		~TopoVertexVertexInFaceIterator() {};
		V* value() {
			
			int adjVId = f->vs[f->vadjv[vIndex][m_iter]];
			V* adjv = m->idVertices(adjVId);
			return adjv;
		}
		void operator++() { m_iter++; };
		void operator++(int) {
			m_iter++; 
		};
		bool end() { return m_iter == 2; };
		V* operator*() { return value(); };
	protected:
		V* v;
		F* f;
		topoM<V, E, F, H>* m;
		int m_iter;
		int vIndex;

	};
	
	/*go through adjcent vertex of a vertex in a hex*/
	template<typename V, typename E, typename F, typename H>
	class TopoVertexVertexInHexIterator
	{
	public:
		TopoVertexVertexInHexIterator(topoM< V, E, F, H>* inputM, V* inputV, H* inputH)
		{
			v = inputV;
			m = inputM;
			h = inputH;
			vIndex = h->vertexIndex(v->id());
			m_iter = 0;
			assert(vIndex != -1);
		};
		~TopoVertexVertexInHexIterator() {};
		V* value()
		{
			int adjVid = h->vs[h->vadjv[vIndex][m_iter]];
			V* adjV = m->idVertices(adjVid);
			return adjV;
		}
		void operator++() { m_iter++; };
		void operator++(int) { m_iter++; };
		bool end() { return m_iter == 3; };
		V* operator*() { return value(); };
	protected:
		V* v;
		H* h;
		topoM< V, E, F, H>* m;
		int vIndex;
		int m_iter;
	};

	/*go through adjcent edges of a vertex in a hex*/
	template<typename V, typename E, typename F, typename H>
	class TopoVertexEdgeInHexIterator
	{
	public:
		TopoVertexEdgeInHexIterator(topoM< V, E, F, H>* inputM, V* inputV, H* inputH)
		{
			v = inputV;
			m = inputM;
			h = inputH;
			vIndex = h->vertexIndex(v->id());
			assert(vIndex != -1);
			m_iter = 0;
		};
		~TopoVertexEdgeInHexIterator() {};
		E* value()
		{
			int adjEid = h->es[h->vadje[vIndex][m_iter]];
			E* adjE = m->idEdges(adjEid);
			return adjE;
		}
		void operator++() { m_iter++; };
		void operator++(int) { m_iter++; };
		bool end() { return m_iter == 3; };
		E* operator*() { return value(); };
	protected:
		V* v;
		H* h;
		topoM< V, E, F, H>* m;
		int vIndex;
		int m_iter;
	};
	
	/*go through adjcent faces of a vertex in a hex*/
	template<typename V, typename E, typename F, typename H>
	class TopoVertexFaceInHexIterator
	{
	public:
		TopoVertexFaceInHexIterator(topoM< V, E, F, H>* inputM, V* inputV, H* inputH)
		{
			v = inputV;
			m = inputM;
			h = inputH;
			vIndex = h->vertexIndex(v->id);
			assert(vIndex != -1);
			m_iter = 0;
		};
		~TopoVertexFaceInHexIterator() {};
		F* value()
		{
			int adjFid = h->es[h->vadjf[vIndex][m_iter]];
			F* adjF = m->idFaces(adjFid);
			return adjF;
		}
		void operator++() { m_iter++; };
		void operator++(int) { m_iter++; };
		bool end() { return m_iter == 3; };
		F* operator*() { return value(); };
	protected:
		V* v;
		H* h;
		topoM< V, E, F, H>* m;
		int vIndex;
		int m_iter;
	};

	/*go through parallel edges of a edge in a hex*/
	template<typename V, typename E, typename F, typename H>
	class TopoEdgeParalelEdgeIterator 
	{
	public:
		TopoEdgeParalelEdgeIterator(topoM< V, E, F, H>* inputM, E* inputE,H* inputH)
		{
			m = inputM;
			e = inputE;
			h = inputH;
			eIndex = h->edgeIndex(e->id());
			assert(eIndex != -1);
			m_iter = 0;
		};
		~TopoEdgeParalelEdgeIterator() {};
		E* value() {
			
			int paraEId = h->es[h->eparalle[eIndex][m_iter]];
			return m->idEdges(paraEId);
		}
		void operator++() { m_iter++; };
		void operator++(int) { m_iter++; };
		bool end() { return m_iter == 3; };
		E* operator*() { return value(); };
	protected:
		E* e;
		H* h;
		topoM< V, E, F, H>* m;
		int eIndex;
		int m_iter = 0;
	};

	/*go through adjcent edge of a edge in a hex*/
	template<typename V, typename E, typename F, typename H>
	class TopoEdgeAdjEdgeInHexIterator
	{
	public:
		TopoEdgeAdjEdgeInHexIterator(topoM< V, E, F, H>* inputM, E* inputE, H* inputH)
		{
			m = inputM;
			e = inputE;
			h = inputH;
			eIndex = h->edgeIndex(e->id);
			assert(eIndex != -1);
			m_iter = 0;
		};
		~TopoEdgeAdjEdgeInHexIterator() {};
		E* value() {

			int adjEId = h->es[h->eadje[eIndex][m_iter]];
			return m->idEdges(adjEId);
		}
		void operator++() { m_iter++; };
		void operator++(int) { m_iter++; };
		bool end() { return m_iter == 4; };
		E* operator*() { return value(); };
	protected:
		E* e;
		H* h;
		topoM< V, E, F, H>* m;
		int eIndex;
		int m_iter = 0;
	};
	
	/*go through adjcent face of a edge in a hex*/
	template<typename V, typename E, typename F, typename H>
	class TopoEdgeFaceInHexIterator
	{
	public:
		TopoEdgeFaceInHexIterator(topoM< V, E, F, H>* inputM, H* inputH, E* inputE)
		{
			m = inputM;
			e = inputE;
			h = inputH;
			eIndex = h->edgeIndex(e->id);
			assert(eIndex != -1);
			m_iter = 0;
		};
		~TopoEdgeFaceInHexIterator() {};
		E* value() {
			int adjFId = h->es[h->eadjf[eIndex][m_iter]];
			return m->idFaces(adjFId);
		}
		void operator++() { m_iter++; };
		void operator++(int) { m_iter++; };
		bool end() { return m_iter == 2; };
		E* operator*() { return value(); };
	protected:
		E* e;
		H* h;
		topoM< V, E, F, H>* m;
		int eIndex;
		int m_iter = 0;
	};

	/*go through adjcent face of a face in a hex*/
	template<typename V, typename E, typename F, typename H>
	class TopoFaceAdjFaceInHexIterator
	{
	public:
		TopoFaceAdjFaceInHexIterator(topoM< V, E, F, H>* inputM, H* inputH, F* inputF )
		{
			m_iter = 0;
			m = inputM;
			f = inputF;
			h = inputH;
			fIndex = h->faceIndex(f->id());
			assert(fIndex != -1);
			m_iter = 0;
		};
		~TopoFaceAdjFaceInHexIterator() {};
		F* value() {

			int adjFid = h->fs[h->fadjf[fIndex][m_iter]];
			return m->idFaces(adjFid);
		}
		void operator++() { m_iter++; };
		void operator++(int) { m_iter++; };
		bool end() { return m_iter == 4; };
		F* operator*() { return value(); };
	protected:
		F* f;
		H* h;
		topoM< V, E, F, H>* m;
		int fIndex;
		int m_iter;
	};

	/*go through adjcent face of a edge in order*/
	template<typename V, typename E, typename F, typename H>
	class TopoEdgeAdjFaceInOrderIterator
	{
	public:
		TopoEdgeAdjFaceInOrderIterator(topoM< V, E, F, H>* inputM, E* inputE)
		{
			miter = 0;
			m = inputM;
			e = inputE;
			startF = NULL;
			for (int i = 0; i < e->neighbor_fs.size(); i++)
			{
				F* f = m->idFaces(e->neighbor_fs[i]);
				if (f->boundary())
				{
					startF = f;
					break;
				}
				startF = f;
			}
			itF = startF;
			itH = m->idHexs(itF->neighbor_hs[0]);
			faceNum = e->neighbor_fs.size();
		};
		~TopoEdgeAdjFaceInOrderIterator() {};
		F* value() { return itF; };
		void operator++() 
		{
			miter++;
			itF = m->flip_f(itH, itF, e);
			if (miter!=faceNum-1)
			{
				itH = m->flip_H(itH, itF);
			}		
		};
		void operator++(int) {
			miter++;
			itF = m->flip_f(itH, itF, e);
			if (miter != faceNum - 1)
			{
				itH = m->flip_H(itH, itF);
			}
		};
		bool end() 
		{ 
			bool result = false;
			if (miter>=faceNum)
			{
				result = true;
			}
			/*if (miter!=0)
			{
				if (itF->boundary())
				{
					
				}
				if (itF==startF)
				{
					result = true;
				}
			}*/
			return result;
		};
		F* operator*() { return value(); };
	protected:
		topoM< V, E, F, H>* m;
		E* e;
		F* startF;
		F* itF;	
		H* itH;
		int miter;
		int faceNum;
	};
}

#endif // TOPO_ITERATOR_H

	