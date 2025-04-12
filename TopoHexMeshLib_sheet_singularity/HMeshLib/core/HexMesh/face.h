/*!
*      \file face.h
*      \brief Base Face Class for all types of Tetrahedron Mesh Classes
*
*		This is the fundamental class for Face
*	   \author David Gu
*      \date 10/01/2011
*
*/

#ifndef _HMESHLIB_FACE_H_
#define _HMESHLIB_FACE_H_

#include <list>

namespace HMeshLib
{

class CVertex;
class CHVertex;
class CHalfEdge;
class CEdge;
class CHalfFace;
class CFace;
//class CTet;
class CHex;

/*!
* \brief CFace, base class for Face
*/
class CFace
{
public:
	CFace()
	{
		m_pLeft = NULL;
		m_pRight= NULL;
	};

	~CFace(){};

	CHalfFace * left() { return m_pLeft; };
	CHalfFace * right(){ return m_pRight;};

	void SetLeft ( CHalfFace * pHf ) { m_pLeft  = pHf; };
	void SetRight( CHalfFace * pHf ) { m_pRight = pHf; };
	bool & boundary(){  return m_IsBoundary;};

	virtual void _from_string() { };
	virtual void _to_string() { };

protected:

	CHalfFace * m_pLeft;
	CHalfFace * m_pRight;
	bool		m_IsBoundary;
};

}

#endif