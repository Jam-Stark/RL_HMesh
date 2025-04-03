/*!
*      \file Hex.h
*      \brief Base tet Class for all types of Tetrahedron Mesh Classes
*
*		This is the fundamental class for tet
*	   \author David Gu
*      \date 10/01/2011
*
*/

#ifndef _HMESHLIB_TET_H_
#define _HMESHLIB_TET_H_

#include <list>

namespace HMeshLib
{

class CVertex;
class CHVertex;
class CHalfEdge;
class CEdge;
class CHalfFace;
class CFace;
class CHex;

/*!
* \brief CHex, base class for Hex
*/
class CHex
{
public:

	CHex()	{ m_iID  = 0; };
	~CHex(){};

	CHalfFace * half_face(int j) { return m_pHalfFace[j]; };
	CHVertex *  hvertex( int j ) { return m_pHVertex[j]; };
	CVertex * vertex(int j) { return m_pHVertex[j]->vert(); };
	void setHVertex(CHVertex * pHV, int j) { m_pHVertex[j] = pHV; };
	void setHalfFace(CHalfFace * pHF, int j) { m_pHalfFace[j] = pHF; };

	int  & id() { return m_iID; };

	virtual void _from_string() { };
	virtual void _to_string() { };

	std::string & string() { return m_string; };

protected:

	CHalfFace  * m_pHalfFace[6];
	CHVertex   * m_pHVertex[8];
	int          m_iID;

	std::string m_string;
};


}

#endif