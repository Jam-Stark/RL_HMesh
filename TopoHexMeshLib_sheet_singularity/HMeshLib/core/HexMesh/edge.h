/*!
*      \file edge.h
*      \brief Base edge Class for all types of Tetrahedron Mesh Classes
*
*		This is the fundamental class for edge
*	   \author David Gu
*      \date 10/01/2011
*
*/

#ifndef _HMESHLIB_EDGE_H_
#define _HMESHLIB_EDGE_H_

#include <list>
#include "../Geometry/Point.h"
#include "../Geometry/Point2.h"

namespace HMeshLib
{

class CVertex;
class CHVertex;
class CHalfEdge;
class CEdge;
class CHalfFace;
class CFace;
class CTet;

/*!
* \brief CEdge, base class for edge
*/
class CEdge	//Edge among tets
{
public:

	CEdge(){ m_vertices[0] = NULL; m_vertices[1] = NULL; };

	~CEdge(){ m_lHEdges.clear();};
	
	std::list<CHEdge*> * edges() { return &m_lHEdges; };

	bool& boundary() { return m_bIsBoundary; };
	bool& sharp() { return m_sharp; };
  
	CVertex *  vertex1() { return m_vertices[0]; };
	CVertex *  vertex2() { return m_vertices[1]; };

	void SetVertex1( CVertex* v ) { m_vertices[0] = v; };
	void SetVertex2( CVertex* v ) { m_vertices[1] = v; };

	virtual void _from_string() { };
	virtual void _to_string() { };

	std::string & string() { return m_string; };

protected:
	std::list<CHEdge*> m_lHEdges;
	CVertex *  m_vertices[2];
	bool     m_bIsBoundary;
	std::string m_string;
	bool m_sharp;
};


}

#endif