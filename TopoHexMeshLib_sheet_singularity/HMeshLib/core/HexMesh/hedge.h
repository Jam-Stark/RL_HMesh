/*!
*      \file tedge.h
*      \brief Base tedge Class for all types of Tetrahedron Mesh Classes
*
*		This is the fundamental class for tedge
*	   \author David Gu
*      \date 10/01/2011
*
*/

#ifndef _HMESHLIB_TEDGE_H_
#define _HMESHLIB_TEDGE_H_

#include <list>
#include "halfedge.h"

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
	* \brief CTEdge, base class for tedge
	*/
	class CHEdge //Edge inside a tetrahedron
	{
	public:

		CHEdge(){ m_pLeft = NULL; m_pRight = NULL; };
		~CHEdge(){};

		CHalfEdge * left() { return m_pLeft; };
		CHalfEdge * right(){ return m_pRight; };
		CEdge     * edge() { return m_pEdge; };
		CHex      * Hex()  { return  m_pHex; };
		CHEdge    * dual() { return  m_pLeft->prev()->prev()->dual()->Hedge();};//hedge in the same hex , against the current one
		std::list<CHEdge*>     parallel() { return m_parallel; };

		void SetLeft(CHalfEdge * pHe) { m_pLeft = pHe; };
		void SetRight(CHalfEdge* pHe) { m_pRight = pHe; };
		void SetEdge(CEdge * pE) { m_pEdge = pE; };
		void SetHex(CHex * pHex)   { m_pHex = pHex; };
		void SetParallel(CHEdge * pHE) { m_parallel.push_back(pHE); };
		
		bool operator==(const CHEdge & e)
		{
			for (int i = 0; i < 2; i++)
			{
				if (m_key[i] != e.m_key[i]) return false;
			}
			return true;
		}

		int & key(int k) { return m_key[k]; };

		virtual void _from_string() { };
		virtual void _to_string() { };

	protected:

		CHalfEdge    * m_pLeft;
		CHalfEdge    * m_pRight;
		CEdge        * m_pEdge;
		CHex         * m_pHex;
		std::list<CHEdge*>  m_parallel;
		int            m_key[2];
	};

}

#endif