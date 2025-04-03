/*!
*      \file tvertex.h
*      \brief Base tvertex Class for all types of hexrahedron Mesh Classes
*
*		This is the fundamental class for tvertex
*	   \author David Gu
*      \date 10/01/2011
*
*/

#ifndef _HMESHLIB_TVERTEX_H_
#define _HMESHLIB_TVERTEX_H_

#include <list>
#include "../Geometry/Point.h"
#include "../Geometry/Point2.h"

namespace HMeshLib
{

	class CVertex;
	class CHalfEdge;
	class CEdge;
	class CHalfFace;
	class CFace;
	class CHex;


	/*!
	* \brief CHVertex, base class for hexrahedron vertex
	*/
	class CHVertex
	{
	public:
		CHVertex(){ m_pVertex = NULL; m_pHex = NULL; m_pHalfedge = NULL; };
		~CHVertex(){};

		int & id() { return m_iID; };

		CVertex         * vert(){ return m_pVertex; };
		CHex            * hex() { return m_pHex; };
		CHalfEdge       * halfedge() { return m_pHalfedge; };

		void set_vert(CVertex * pV) { m_pVertex = pV; };
		void set_hex(CHex        * pT) { m_pHex = pT; };
		void set_halfedge(CHalfEdge * pH) { m_pHalfedge = pH; };

		virtual void _from_string() { };
		virtual void _to_string() { };

	protected:
		//vertex ID
		int            m_iID;
		CVertex      * m_pVertex;
		CHex         * m_pHex;
		//outgoing, halfedge start from this HVertex
		CHalfEdge *    m_pHalfedge;
	};

}

#endif