/*!
*      \file vertex.h
*      \brief Base vertex Class for all types of Tetrahedron Mesh Classes
*
*		This is the fundamental class for vertex
*	   \author David Gu
*      \date 10/01/2011
*
*/

#ifndef _HMESHLIB_VERTEX_H_
#define _HMESHLIB_VERTEX_H_

#include <list>
#include "../Geometry/Point.h"

namespace HMeshLib
{

class CHVertex;
class CHalfEdge;
class CEdge;
class CHEdge;
class CHalfFace;
class CFace;
class CHex;

class CVertex
{
public:
	CVertex() { m_iID = 0; };
	~CVertex(){ m_pHVertices.clear();}

	CPoint & position() { return m_vPosition; };
	int    & id()       { return m_iID;       };

	std::list<CEdge*>   * edges() { return &m_pEdges; };
	std::list<CHEdge*>  * hedges(){ return &m_pHEdges;};

	std::list<CHalfFace*>  * HalfFaces(){ return &m_pHFaces;};
	std::list<CHVertex*>* hvertices() { return &m_pHVertices;};

	std::string&        string(){ return m_string; };

	virtual void _from_string() { };
	virtual void _to_string() { };

protected:

	CPoint m_vPosition;
	int    m_iID;
	bool   m_bIsBoundary;

	std::list<CHalfFace*>  m_pHFaces;		//temporary HalfFace list, will be empty after loading the whole mesh 
	std::list<CHEdge*>     m_pHEdges;		//temporary TEdge list, will be empty after loading the whole mesh

	std::list<CHVertex*>   m_pHVertices;	//adjacent HVertecies
	std::list<CEdge*>      m_pEdges;	    //adjacent Edges;

	std::string m_string;
};
};

#endif