/*!
*      \file halfface.h
*      \brief Base HalfFace Class for all types of Tetrahedron Mesh Classes
*
*		This is the fundamental class for HalfFace
*	   \author David Gu
*      \date 10/01/2011
*
*/

#ifndef _HMESHLIB_HALFFACE_H_
#define _HMESHLIB_HALFFACE_H_

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
* \brief CHalfFace, base class for HalfFace
*/
class CHalfFace
{
public:
	CHalfFace()
	{
		m_pHalfEdge = NULL;
		m_pFace     = NULL;
		m_pHex      = NULL;
		m_pDual    = NULL;
	};

	~CHalfFace(){};

	CHalfEdge * half_edge() { return m_pHalfEdge; };
	CFace     * face()      { return m_pFace;     };
	CHex      * hex()       { return m_pHex;      };
	CHalfFace * dual()      { return m_pDual;     };
	CHalfFace * parallel()  { return m_parallel;  };
	int       & key( int k) { return m_key[k];    };

	void SetHalfEdge( CHalfEdge * pHe ) { m_pHalfEdge = pHe; };
	void SetFace( CFace * pF )          { m_pFace = pF;      };
	void SetHex(  CHex  * pHex )        { m_pHex  = pHex;      };
	void SetDual( CHalfFace   *pF )     { m_pDual = pF;     };
	void SetParallel( CHalfFace  *pF)   { m_parallel = pF;};
	bool operator==( const CHalfFace & f )
	{
		for( int i = 0; i < 4; i ++ )
		  if( m_key[i] != f.m_key[i] ) return false;
		return true;
	};

  
 
protected:

	CHalfEdge * m_pHalfEdge;
	CFace     * m_pFace; 
	CHex      * m_pHex;
	CHalfFace * m_pDual;
	CHalfFace * m_parallel;
	int         m_key[4];
};

}

#endif