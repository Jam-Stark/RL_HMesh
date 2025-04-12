/*!
*      \file titerators.h
*      \brief iterators for Tetrahedron Mesh Classes
*
*	   \author David Gu
*      \date 10/02/2011
*
*/

#ifndef _HMESHLIB_TITERATORS_H_
#define _HMESHLIB_TITERATORS_H_

#include <list>
#include <set>
#include "vertex.h"
#include "hvertex.h"
#include "edge.h"
#include "hedge.h"
#include "halfedge.h"
#include "face.h"
#include "halfface.h"
#include "hex.h"
#include "BasehMesh.h"

namespace HMeshLib{

/*!
\brief HMeshVertex Iterator

go through all the vertices in the TetMesh
*/
template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
class HMeshVertexIterator
{
public:
	HMeshVertexIterator( CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> *  pMesh )
	{ 
		m_pMesh = pMesh; 
		//m_id = 0;
		m_iter = pMesh->vertices().begin();
	};

	~HMeshVertexIterator(){};

	CVertex * value() { return *m_iter; };

	/*! iterate ++ */
	void operator++()
	{ m_iter ++; };
	/*! ++iterator */
	void operator++(int) { m_iter ++; };
	/*! whether the iterator reaches the end */
	bool end(){ return m_iter == m_pMesh->vertices().end(); };
	/*! get the current vertex */
	CVertex * operator*() { return value(); };

protected:
	/*! HexMesh pointer */
	CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> *        m_pMesh;
	/*! Current vertex list iterator */
	//int             m_id;
	typename std::list<CVertex*>::iterator m_iter;
};

/*!
\brief VertexVertex Iterator

go through all the adjacent vertices in of a vertex
*/
template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
class HVertexVertexIterator
{
public:
	/*! Constructor */
	HVertexVertexIterator( CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> *  pMesh, CVertex *  pV )
	{ 
		m_pMesh = pMesh; 
		m_vertex = pV; 
		std::list<CEdge*> * pEdges = pMesh->VertexEdgeList( pV );
		std::set<CVertex*> vertices;

		for(std::list<CEdge*>::iterator eiter = (*pEdges).begin(); eiter != (*pEdges).end(); eiter ++ )
		 {
		   CEdge * e = *eiter;
		   CVertex * v1 = pMesh->EdgeVertex1( e );
		   CVertex * v2 = pMesh->EdgeVertex2( e );
		   if( v1 != pV ) 
			 m_pVertices.insert( v1 );
		   else
			 m_pVertices.insert( v2 );
		 }
		m_iter = m_pVertices.begin();
	};
	/*! Destructor */
	~HVertexVertexIterator(){};
	/*! Iterator ++ */
	void operator++()
	{ m_iter++; };
	/*! ++ Iterator */
	void operator++(int) { m_iter ++; };
	/*! Verify if the iterator reaches the end */
	bool end(){ return m_iter == m_pVertices.end(); };
	/*! dereferencing */
	CVertex * operator*() { return *m_iter; };

protected:
	/*! Pointer to the tmesh */
	CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> *  m_pMesh;
	/*! center vertex */
	CVertex *			m_vertex;
	/*! surrounding vertices */
	std::set<CVertex*>	m_pVertices;
	/*! set iterator */
	typename std::set<CVertex*>::iterator m_iter;
  };

/*!
\brief MeshHex Iterator

go through all the hexs in a mesh
*/
template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
class HMeshHexIterator
{
public:
	/*! Consructor */
	HMeshHexIterator( CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> *  pMesh )
	{ 
		m_pMesh = pMesh; 
		//m_id = 0;
		m_iter = pMesh->Hexs().begin();
	};
	/*! Destructor */
	~HMeshHexIterator(){};
	/*! ++iterate */
	void operator++()
	{ m_iter ++; };
	/*! iterator ++ */
	void operator++(int) { m_iter ++; };
	/*! verify if reaches the end */
	 bool end(){ return m_iter == m_pMesh->Hexs().end(); };
	 /*! dereferencing */
	CHex * operator*() { return *m_iter; };

private:
	/*! pointer to a TetMesh */
	CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> *        m_pMesh;
	/*! tet iterator */
	typename std::list<CHex*>::iterator m_iter;
};

/*!
\brief HexHalfFace Iterator

go through all the halfaces of a Hex in a mesh
*/

template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
class HexHalfFaceIterator
{
public:
	/*! Constructor */
	HexHalfFaceIterator( CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> * pMesh, CHex * pT )
	{ 
		m_pMesh = pMesh;
		m_pHex = pT; 
		m_id = 0;
	};
	/*! Destructor */
	~HexHalfFaceIterator(){};
	/*! iterator ++ */
	void operator++()
	{ m_id ++; };
	/*! ++iterator */
	void operator++(int) { m_id ++; };
	/*! verify if the end has been reached */
	bool end(){ return m_id == 6; };
	/*! dereferencing */
	CHalfFace * operator*() { return m_pMesh->HexHalfFace( m_pHex, m_id ); };

private:
	/*! pointer to a TetMesh */
	CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> *  m_pMesh;
	/*! current tet */
	CHex *             m_pHex;
	/*! current half_face */
	int                m_id;
};

/*!
\brief Mesh Edge Iterator
go through all the edges of a Hex in a mesh
*/
template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
class HMeshEdgeIterator
{
public:
	/*! Constructor */
  HMeshEdgeIterator( CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge, CEdge, CHalfFace,CFace, CHex> * pMesh )
  {
    m_pMesh = pMesh;
    m_iter = m_pMesh->edges().begin();
  };
	/*! Destructor */
	~HMeshEdgeIterator(){};
	/*! ++ iterator */
	void operator++()
	{ m_iter ++; };
	/*! iterator++ */
	void operator++(int) { m_iter ++; };
	/*! verify if the end has been reached */
	bool end() { return m_iter == m_pMesh->edges().end();};
	/*! dereferencing */
	CEdge * operator*() { return *m_iter;};
private:
	/*! pointer to the mesh */
	CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> * m_pMesh;
	/*! edge list iterator */
	typename std::list<CEdge*>::iterator m_iter;
};

/*!
\brief Edge->TEdge Iterator
go through all the HEdges of an edge
*/
template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
class EdgeHEdgeIterator
{
public:
	/*! constrructor */
  EdgeHEdgeIterator( CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> * pMesh, CEdge * pE )
  {
    m_pMesh = pMesh;
    m_pEdge = pE;
    m_iter     = m_pMesh->EdgeHEdgeList(pE)->begin();
  };
  /*! Destructor */
  ~EdgeHEdgeIterator(){};
  /*! ++iterator */		
  void operator++(){ m_iter ++;};
  /*! iterator++ */
  void operator++(int unused) { m_iter ++; };
  /*! verify if the end has been reached */	
  bool end(){ return m_iter == m_pMesh->EdgeHEdgeList(m_pEdge)->end();};
  /*! dereferencing */
  CHEdge * operator*() { return *m_iter;};
private:
  /*! pointer to mesh */
  CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> * m_pMesh;
  /*! pointer to the edge */
  CEdge  * m_pEdge;
  /*! TEdge list iterator */
  typename std::list<CHEdge*>::iterator m_iter;
};

/*!
\brief Vertex->Edge Iterator
go through all the edges of a vertex
*/
template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
class HVertexEdgeIterator
{
public:
	/*! constructor */
  HVertexEdgeIterator( CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> * pMesh, CVertex * pV )
  {  m_pMesh = pMesh; m_pV = pV; m_iter = m_pMesh->VertexEdgeList(pV)->begin(); };
  /*! destructor */
  ~HVertexEdgeIterator() {};
  /*! ++ iterator */	
  void operator++(){ m_iter ++; };
  /*! iterator ++ */
  void operator++(int unused) { m_iter ++; };
  /*! verify if the end has been reached */
  bool end() { return m_iter == m_pMesh->VertexEdgeList(m_pV)->end();};
  /*! dereferencing */
  CEdge * operator*() { return *m_iter;};

private:
  /*! pointer to mesh */
  CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> * m_pMesh;
  /*! pointer to the vertex */	
  CVertex * m_pV;
  /*! Edge list iterator */
  typename std::list<CEdge*>::iterator m_iter;
};

/*!
\brief Hex->Edge Iterator
go through all the edges of a Tet
*/
template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
class HexEdgeIterator
{
public:
	/*! constructor */
  HexEdgeIterator( CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> * pMesh,  CHex * pT )
  {
	m_pMesh = pMesh;
    m_pT = pT;  
    for( int i = 0; i < 6; i ++ )
    {
      CHalfFace * pF = m_pMesh->HexHalfFace(m_pT, i);
      CHalfEdge * pH = m_pMesh->HalfFaceHalfEdge( pF );
      for( int j = 0; j < 4; j ++ )
      {
			CHEdge * pHE = m_pMesh->HalfEdgeHEdge( pH );
			CEdge  * pE  = m_pMesh->HEdgeEdge( pHE );
			m_edges.insert( pE );
			pH = m_pMesh->HalfEdgeNext( pH );
		}
	}
	m_iter = m_edges.begin();
   };
    /*! destructor */
	~HexEdgeIterator() {};
	/*! iterator ++ */
	void operator++(){ m_iter ++; };
	/*! ++ iterator */
	void operator++(int unused) { m_iter ++; };
	/*! verify if the end of the list has been reached */
	bool end() { return m_iter == m_edges.end();};
	/*! dereferencing */
	CEdge * operator*() { return *m_iter;};

private:
	/*! Pointer to the mesh */
	CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> * m_pMesh;
	/*! Pointer to the tet */
	CHex * m_pT;
	/*! set of edges adjacent to the central vertex */
	std::set<CEdge*> m_edges;
	/*! edge set iterator */
	typename std::set<CEdge*>::iterator m_iter;
};

/*!
\brief Mesh->Face Iterator
go through all the faces of a Mesh
*/
template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
class HMeshFaceIterator
{
public:
	/*! Constructor */
	HMeshFaceIterator( CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> *  pMesh )
	{ 
		m_pMesh = pMesh; 
		m_iter = pMesh->faces().begin();
	};
	/*! Destructor */
	~HMeshFaceIterator(){};
	/*! ++ iterator */
	void operator++(){ m_iter ++; };
	/*! iterator ++ */
	void operator++(int) { m_iter ++; };
	/*! reach the end of face list */
	bool end(){ return m_iter == m_pMesh->faces().end(); };
	/*! dereferencing */
	CFace * operator*() { return *m_iter; };

private:
	/*! pointer to the mesh */
	CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> *    m_pMesh;
	/*! list of face iterator */
	typename std::list<CFace*>::iterator m_iter;
};

/*!
\brief Edge->Face Iterator
go through all the faces adjacent to an edge
*/
template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
class HEdgeFaceIterator
{
public:
	/*! constructor */
	HEdgeFaceIterator( CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> * pMesh, CEdge * pe )
	{ 
		m_pMesh = pMesh;
		std::list<CHEdge*> * HEdges = m_pMesh->EdgeHEdgeList( pE );
		for( list<CHEdge*>::iterator tit = HEdges->begin();  tit != HEdges->end();  tit ++ )
        {
            CHEdge * pT  = *tit;
            
            CHalfEdge * pL = m_pMesh->HEdgeLeft( pT );
            CHalfEdge * pR = m_pMesh->HEdgeRight( pT );
            
            CHalfFace * pHF = m_pMesh->HalfEdgeHalfFace( pL );
            CFace *      pF = m_pMesh->HalfFaceFace( pHF );
            m_faces.insert( pF );

			pHF = m_pMesh->HalfEdgeHalfFace( pR );
            pF  = m_pMesh->HalfFaceFace( pHF );
            m_faces.insert( pF );
        }

		m_iter = m_faces.begin();
	};
	
	/*! destructor */
	~HEdgeFaceIterator(){};
	/*! ++ iterator */
	void operator++() { m_iter ++; };
	/*! iterator ++ */
	void operator++(int unused) { m_iter ++; };
	/*! reach the end */
	bool end(){ return m_iter == m_faces.end(); };
	/*! dereferencing */
	CFace * operator*() { return *m_iter; };

private:
	/*! pointer to the mesh */
	CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> *    m_pMesh;
	/*! face set */
	std::set<CFace*> m_faces;
	/*! face set iterator */
	typename std::set<CFace*>::iterator m_iter;
};


/*!
\brief HalfFace->vertex Iterator
go through all the vertices of a HalfFace
*/
template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
class HalfFaceVertexIterator
{
public:
	/*! constructor */
  HalfFaceVertexIterator( CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> * pMesh,  CHalfFace * pF )
  {
	m_pMesh = pMesh;
    m_pF = pF;
	m_pH = pMesh->HalfFaceHalfEdge( pF );
   };
    /*! destructor */
	~HalfFaceVertexIterator() {};
	/*! iterator ++ */
	void operator++()
	{ 
		m_pH = m_pMesh->HalfEdgeNext( m_pH ); 
		if( m_pH == m_pMesh->HalfFaceHalfEdge( m_pF ) )
		{
			m_pH = NULL;
		}
	};
	/*! ++ iterator */
	void operator++(int)
	{ 
		m_pH = m_pMesh->HalfEdgeNext( m_pH ); 
		if( m_pH == m_pMesh->HalfFaceHalfEdge( m_pF ) )
		{
			m_pH = NULL;
		}	
	};
	/*! verify if the end of the list has been reached */
	bool end() { return m_pH == NULL; };
	/*! dereferencing */
	CVertex * operator*() { return m_pMesh->HalfEdgeTarget( m_pH );};

private:
	/*! Pointer to the mesh */
	CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> * m_pMesh;
	/*! Pointer to the halfface */
	CHalfFace * m_pF;
	/*! Pointer to the halfedge */
	CHalfEdge * m_pH;
};

/*!
\brief HalfFace->HalfEdge Iterator
go through all the halfedges of a HalfFace
*/
template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
class HalfFaceHalfEdgeIterator
{
public:
	/*! constructor */
   HalfFaceHalfEdgeIterator( CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> * pMesh,  CHalfFace * pF )
   {
	 m_pMesh = pMesh;
     m_pF = pF;
	 m_pH = pMesh->HalfFaceHalfEdge( pF );
   };
    /*! destructor */
	~HalfFaceHalfEdgeIterator() {};
	/*! iterator ++ */
	void operator++()
	{ 
		m_pH = m_pMesh->HalfEdgeNext( m_pH ); 
		if( m_pH == m_pMesh->HalfFaceHalfEdge( m_pF ) )
		{
			m_pH = NULL;
		}
	};
	/*! ++ iterator */
	void operator++(int unused ) 
	{ 
		m_pH = m_pMesh->HalfEdgeNext( m_pH ); 
		if( m_pH == m_pMesh->HalfFaceHalfEdge( m_pF ) )
		{
			m_pH = NULL;
		}	
	};
	/*! verify if the end of the list has been reached */
	bool end() { return m_pH == NULL; };
	/*! dereferencing */
	CHalfEdge * operator*() { return m_pH; };

protected:
	/*! Pointer to the mesh */
	CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> * m_pMesh;
	/*! Pointer to the halfface */
	CHalfFace * m_pF;
	/*! Pointer to the halfedge */
	CHalfEdge * m_pH;
};


/*!
\brief Face->vertex Iterator
go through all the vertices of a Face
*/
template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
class FaceVertexIterator
{
public:
	/*! constructor */
  FaceVertexIterator( CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> * pMesh,  CFace * pF )
  {
	m_pMesh = pMesh;
    m_pF = pF;
	m_pHF = m_pMesh->FaceLeftHalfFace( m_pF );
	m_pH = pMesh->HalfFaceHalfEdge( m_pHF );
   };
    /*! destructor */
	~FaceVertexIterator() {};
	/*! iterator ++ */
	void operator++()
	{ 
		m_pH = m_pMesh->HalfEdgeNext( m_pH ); 
		if( m_pH == m_pMesh->HalfFaceHalfEdge( m_pHF ) )
		{
			m_pH = NULL;
		}
	};
	/*! ++ iterator */
	void operator++(int) 
	{ 
		m_pH = m_pMesh->HalfEdgeNext( m_pH ); 
		if( m_pH == m_pMesh->HalfFaceHalfEdge( m_pHF ) )
		{
			m_pH = NULL;
		}	
	};
	/*! verify if the end of the list has been reached */
	bool end() { return m_pH == NULL; };
	/*! dereferencing */
	CVertex * operator*() { return m_pMesh->HalfEdgeTarget( m_pH );};

private:
	/*! Pointer to the mesh */
	CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> * m_pMesh;
	/*! Pointer to the face */
	CFace * m_pF;
	/*! Pointer to the halfface */
	CHalfFace * m_pHF;
	/*! Pointer to the halfedge */
	CHalfEdge * m_pH;
};


/*!
\brief TVertex->InHalfEdge Iterator
go through all the InHalfEdge of a TVertex
*/
template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
class HVertexInHalfEdgeIterator
{
public:
	/*! constructor */
  HVertexInHalfEdgeIterator( CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> * pMesh,  CHVertex * pTV )
  {
	m_pMesh = pMesh;
    m_pTV = pTV;
	m_pH = pMesh->HVertexHalfEdge( pTV );
   };
    /*! destructor */
	~HVertexInHalfEdgeIterator() {};
	/*! iterator ++ */
	void operator++()
	{ 
		CHalfEdge * pN = m_pMesh->HalfEdgeNext(m_pH);
		m_pH = m_pMesh->HalfEdgeDual(pN);
		if(m_pH == m_pMesh->TVertexHalfEdge(m_pTV))
		{
			m_pH = NULL;
		}
	};
	/*! ++ iterator */
	void operator++(int unused) 
	{ 
		CHalfEdge * pN = m_pMesh->HalfEdgeNext( m_pH );
		m_pH = m_pMesh->HalfEdgeDual( pN ); 
		if( m_pH == m_pMesh->HVertexHalfEdge( m_pTV ) )
		{
			m_pH = NULL;
		}	
	};
	/*! verify if the end of the list has been reached */
	bool end() { return m_pH == NULL; };
	/*! dereferencing */
	CHalfEdge * operator*() { return m_pH; };

private:
	/*! Pointer to the mesh */
	CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> * m_pMesh;
	/*! Pointer to the TVertex */
	CHVertex * m_pTV;
	/*! Pointer to the halfedge */
	CHalfEdge * m_pH;
};

/*!
\brief TVertex->TEdge Iterator
go through all the TEdge of a TVertex
*/
template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
class HVertexHEdgeIterator
{
public:
	/*! constructor */
  HVertexHEdgeIterator( CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> * pMesh,  CHVertex * pTV )
  {
	m_pMesh = pMesh;
    m_pTV = pTV;
	m_pH = pMesh->HVertexHalfEdge( pTV );
   };
    /*! destructor */
	~HVertexHEdgeIterator() {};
	/*! iterator ++ */
	void operator++()
	{ 
		CHalfEdge * pN = m_pMesh->HalfEdgeNext( m_pH );
		m_pH = m_pMesh->HalfEdgeDual( pN ); 
		if( m_pH == m_pMesh->HVertexHalfEdge( m_pTV ) )
		{
			m_pH = NULL;
		}
	};
	/*! ++ iterator */
	void operator++(int unused) 
	{ 
		CHalfEdge * pN = m_pMesh->HalfEdgeNext( m_pH );
		m_pH = m_pMesh->HalfEdgeDual( pN ); 
		if( m_pH == m_pMesh->HVertexHalfEdge( m_pTV ) )
		{
			m_pH = NULL;
		}	
	};
	/*! verify if the end of the list has been reached */
	bool end() { return m_pH == NULL; };
	/*! dereferencing */
	CHEdge * operator*() { return m_pMesh->HalfEdgeTEdge( m_pH ); };

private:
	/*! Pointer to the mesh */
	CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> * m_pMesh;
	/*! Pointer to the TVertex */
	CHVertex * m_pTV;
	/*! Pointer to the halfedge */
	CHalfEdge * m_pH;
};

/*!
\brief Vertex->TVertex Iterator
go through all the tvertices of a vertex
*/
template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
class VertexTVertexIterator
{
public:
	/*! constrructor */
  VertexTVertexIterator( CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> * pMesh, CVertex * pV )
  {
    m_pMesh = pMesh;
    m_pVert = pV;
    m_iter     = m_pMesh->VertexTVertexList(pV)->begin();
  };
  /*! Destructor */
  ~VertexTVertexIterator(){};
  /*! ++iterator */		
  void operator++(){ m_iter ++;};
  /*! iterator++ */
  void operator++(int unused) { m_iter ++; };
  /*! verify if the end has been reached */	
  bool end(){ return m_iter == m_pMesh->VertexTVertexList(m_pVert)->end();};
  /*! dereferencing */
  CHVertex * operator*() { return *m_iter;};
private:
  /*! pointer to mesh */
  CHMesh<CHVertex,CVertex,CHalfEdge,CHEdge,CEdge,CHalfFace,CFace,CHex> * m_pMesh;
  /*! pointer to the vertex */
  CVertex  * m_pVert;
  /*! TVertex list iterator */
  typename std::list<CHVertex*>::iterator m_iter;
};

/*brief Face->Vertex Iterator
* go through all the tvertices of a vertex
*/

}
#endif