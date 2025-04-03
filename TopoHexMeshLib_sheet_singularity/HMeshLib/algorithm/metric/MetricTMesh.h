#ifndef _TMESHLIB_METRIC_MESH_H_
#define _TMESHLIB_METRIC_MESH_H_

#include <stdio.h>

#include "..\..\core\TetMesh\BaseTMesh.h"
#include "..\..\core\TetMesh\titerators.h"
#include "..\..\core\Geometry\plane.h"

#ifndef PI
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679
#endif

namespace TMeshLib
{
/*!
* \brief CMetricTVertex, Vertex for Metric TMesh
*/
class CMetricTVertex : public CTVertex
{
public:
	CMetricTVertex() { m_solid_angle = 0; };
	~CMetricTVertex(){};
	double & solid_angle() { return m_solid_angle; };
protected:
	double m_solid_angle;
};

/*!
* \brief CMetricVertex, Vertex for Metric TMesh
*/
class CMetricVertex : public CVertex
{
public:
	CMetricVertex() { m_curvature = 0; };
	~CMetricVertex(){};
	double & k() { return m_curvature; };
protected:
	double m_curvature;
};

/*!
* \brief CMetricHalfEdge, HalfEdge for Metric TMesh
*/
class CMetricHalfEdge : public CHalfEdge
{
public:
	CMetricHalfEdge() { m_corner_angle = 0; };
	~CMetricHalfEdge(){};
	double & corner_angle() { return m_corner_angle; };
protected:
	double m_corner_angle;
};

/*!
* \brief CMetricTEdge, TEdge for Metric TMesh
*/
class CMetricTEdge : public CTEdge
{
public:
	CMetricTEdge() { m_dihedral_angle = 0; };
	~CMetricTEdge(){};
	double & dihedral_angle() { return m_dihedral_angle = 0; };
protected:
	double m_dihedral_angle;
};

/*!
* \brief CMetricEdge, Edge for Metric TMesh
*/
class CMetricEdge : public CEdge
{
public:
	CMetricEdge() { m_curvature = 0; m_length = 0; };
	~CMetricEdge(){};
	double & edge_length() { return m_length; };
	double & k() { return m_curvature; };
protected:
	double m_curvature;
	double m_length;
};

/*!
* \brief CMetricHalfFace, HalfFace for Metric TMesh
*/
class CMetricHalfFace : public CHalfFace
{
public:
	CMetricHalfFace(){};
	~CMetricHalfFace(){};
	CPoint & normal() { return m_normal; };
protected:
	CPoint m_normal;
};
/*!
* \brief CMetricFace, Face for Metric TMesh
*/
class CMetricFace : public CFace
{
};
/*!
* \brief CMetricTet, Tet for Metric TMesh
*/
class CMetricTet : public CTet
{
public:
	CMetricTet() {};
	~CMetricTet(){};
};
/*!
 *	\brief CMetricTMesh class
 *
 *	TMesh class for Metric
 */
template<typename CTVertex, typename CVertex, typename CHalfEdge,typename CTEdge, typename CEdge, typename CHalfFace, typename CFace, typename CTet>
class CMetricTMesh : public CTMesh<CTVertex, CVertex, CHalfEdge,CTEdge, CEdge, CHalfFace, CFace, CTet>
{
public:
	typedef TMeshTetIterator<CTVertex, CVertex, CHalfEdge,CTEdge, CEdge, CHalfFace, CFace, CTet> MeshTetIterator;
	typedef TetHalfFaceIterator<CTVertex, CVertex, CHalfEdge,CTEdge, CEdge, CHalfFace, CFace, CTet> TetHalfFaceIterator;
	typedef TMeshEdgeIterator<CTVertex, CVertex, CHalfEdge,CTEdge, CEdge, CHalfFace, CFace, CTet> MeshEdgeIterator;
	typedef TMeshVertexIterator<CTVertex, CVertex, CHalfEdge,CTEdge, CEdge, CHalfFace, CFace, CTet> MeshVertexIterator;
	typedef TMeshFaceIterator<CTVertex, CVertex, CHalfEdge,CTEdge, CEdge, CHalfFace, CFace, CTet> MeshFaceIterator;
	typedef HalfFaceVertexIterator<CTVertex, CVertex, CHalfEdge,CTEdge, CEdge, CHalfFace, CFace, CTet> HalfFaceVertexIterator;
	typedef HalfFaceHalfEdgeIterator<CTVertex, CVertex, CHalfEdge,CTEdge, CEdge, CHalfFace, CFace, CTet> HalfFaceHalfEdgeIterator;
	typedef TVertexInHalfEdgeIterator<CTVertex, CVertex, CHalfEdge,CTEdge, CEdge, CHalfFace, CFace, CTet> TVertexInHalfEdgeIterator;
	typedef TVertexTEdgeIterator<CTVertex, CVertex, CHalfEdge,CTEdge, CEdge, CHalfFace, CFace, CTet>  TVertexTEdgeIterator;
	typedef EdgeTEdgeIterator<CTVertex, CVertex, CHalfEdge,CTEdge, CEdge, CHalfFace, CFace, CTet>  EdgeTEdgeIterator;
	typedef VertexTVertexIterator<CTVertex, CVertex, CHalfEdge,CTEdge, CEdge, CHalfFace, CFace, CTet> VertexTVertexIterator;
public:
	/*! Compute the edge length from the vertex position */
	void _position_to_edge_length();
	/*! Compute the corner angle from edge length */
	void _edge_length_to_corner_angle();
	/*! Compute the dihedral angle from corner angle */
	void _corner_angle_to_dihedral_angle();
	/*! Compute the solid angle from dihedral angle */
	void _dihedral_angle_to_solid_angle();
	/*! Compute the edge curvature from dihedral angle */
	void _edge_curvature();
	/*! Compute the vertex curvature from solid angle */
	void _vertex_curvature();
};

/* From vertex position to edge length */
template<typename CTVertex, typename CVertex, typename CHalfEdge,typename CTEdge, typename CEdge, typename CHalfFace, typename CFace, typename CTet>
void CMetricTMesh<CTVertex, CVertex, CHalfEdge,CTEdge, CEdge, CHalfFace, CFace, CTet>::_position_to_edge_length()
{
	for( TMeshEdgeIterator<CTVertex, CVertex, CHalfEdge,CTEdge, CEdge, CHalfFace, CFace, CTet> eiter( this ); !eiter.end(); eiter ++ )
	{
		CEdge * pE = *eiter;
		CVertex * pV1 = EdgeVertex1( pE );
		CVertex * pV2 = EdgeVertex2( pE );
		pE->edge_length() = ( pV1->position() - pV2->position() ).norm();
	}
};

/* From edge length to corner angle*/
template<typename CTVertex, typename CVertex, typename CHalfEdge,typename CTEdge, typename CEdge, typename CHalfFace, typename CFace, typename CTet>
void CMetricTMesh<CTVertex, CVertex, CHalfEdge,CTEdge, CEdge, CHalfFace, CFace, CTet>::_edge_length_to_corner_angle()
{
	for( MeshTetIterator titer( this ); !titer.end(); titer ++ )
	{
		CTet * pT = *titer;
		
		for( TetHalfFaceIterator thfiter( this, pT ); !thfiter.end(); thfiter ++ )
		{
			CHalfFace * pF = *thfiter;
			std::vector<CHalfEdge*> hes;
			std::vector<CEdge*> edges;
			for( HalfFaceHalfEdgeIterator hiter(this, pF); !hiter.end(); hiter ++ )
			{
				CHalfEdge * pH = *hiter;
				hes.push_back( pH );
				CTEdge * pTE =  HalfEdgeTEdge( pH );
				CEdge  * pE = TEdgeEdge( pTE );
				edges.push_back( pE );
			}

			for( int i = 0; i < 3; i ++ )
			{
				double a = edges[(i+0)%3]->edge_length();
				double b = edges[(i+1)%3]->edge_length();
				double c = edges[(i+2)%3]->edge_length();
				
				double C = ( a * a + b * b - c * c )/(2 * a * b );
				hes[i]->corner_angle() = C;
			}
		}
	}
};

/* From corner angle to dihedral angle*/
template<typename CTVertex, typename CVertex, typename CHalfEdge,typename CTEdge, typename CEdge, typename CHalfFace, typename CFace, typename CTet>
void CMetricTMesh<CTVertex, CVertex, CHalfEdge,CTEdge, CEdge, CHalfFace, CFace, CTet>::_corner_angle_to_dihedral_angle()
{
	for( MeshTetIterator titer( this ); !titer.end(); titer ++ )
	{
		CTet * pT = *titer;
		for( int k = 0; k < 4; k ++ )
		{
			CTVertex * pTV = TetTVertex( pT, k );
			std::vector<CHalfEdge*> hes;
			for( TVertexInHalfEdgeIterator hiter( this, pTV ); !hiter.end(); hiter ++ )
			{
				CHalfEdge * pH = *hiter;
				hes.push_back( pH );
			}
			for( int i = 0; i < 3; i ++ )
			{
				double a = hes[(i+0)%3]->corner_angle();
				double b = hes[(i+1)%3]->corner_angle();
				double c = hes[(i+2)%3]->corner_angle();
				
				//spherical cosine law
				double C = acos( ( cos( c ) - cos(a) * cos( b )) /( sin(a ) * sin(b) ) );
				CTEdge * pTE = HalfEdgeTEdge( hes[(i+1)%3] );
				pTE->dihedral_angle() = C;
			}
		}
	}
};

/* From dihedral angle to solid angle*/
template<typename CTVertex, typename CVertex, typename CHalfEdge,typename CTEdge, typename CEdge, typename CHalfFace, typename CFace, typename CTet>
void CMetricTMesh<CTVertex, CVertex, CHalfEdge,CTEdge, CEdge, CHalfFace, CFace, CTet>::_dihedral_angle_to_solid_angle()
{
	for( MeshTetIterator titer( this ); !titer.end(); titer ++ )
	{
		CTet * pT = *titer;
		for( int k = 0; k < 4; k ++ )
		{
			CTVertex * pTV = TetTVertex( pT, k );
			double angle = 0;

			for( TVertexTEdgeIterator eiter( this, pTV ); !eiter.end(); eiter ++ )
			{
				CTEdge * pE = *eiter;
				angle += pE->dihedral_angle();
			}
			pTV->solid_angle() = angle - 3.1415926535;
		}
	}
};

/* From dihedral angle to edge curvature*/
template<typename CTVertex, typename CVertex, typename CHalfEdge,typename CTEdge, typename CEdge, typename CHalfFace, typename CFace, typename CTet>
void CMetricTMesh<CTVertex, CVertex, CHalfEdge,CTEdge, CEdge, CHalfFace, CFace, CTet>::_edge_curvature()
{
	for( MeshEdgeIterator eiter( this ); !eiter.end(); eiter ++ )
	{
		CEdge * pE = *eiter;
		pE->k() = 2 * PI;
		for( EdgeTEdgeIterator teiter( this, pE); !teiter.end(); teiter ++ )
		{
			CTEdge * pTE = * teiter;
			pE->k() -= pTE->dihedral_angle();
		}
	}
};

/* From dihedral angle to edge curvature*/
template<typename CTVertex, typename CVertex, typename CHalfEdge,typename CTEdge, typename CEdge, typename CHalfFace, typename CFace, typename CTet>
void CMetricTMesh<CTVertex, CVertex, CHalfEdge,CTEdge, CEdge, CHalfFace, CFace, CTet>::_vertex_curvature()
{
	for( MeshVertexIterator viter( this ); !viter.end(); viter ++ )
	{
		CVertex * pV = *viter;
		pV->k() = 4 * PI;
		for( VertexTVertexIterator vviter( this, pV ); !vviter.end(); vviter ++ )
		{
			CTVertex * pTV = * vviter;
			pV->k() -= pTV->solid_angle();
		}
	}
};

typedef CMetricTMesh<CMetricTVertex, CMetricVertex, CMetricHalfEdge,CMetricTEdge, CMetricEdge, CMetricHalfFace, CMetricFace, CMetricTet> CMTMesh;


};
#endif