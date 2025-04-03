/*!
*      \file BaseHMesh.h
*      \brief Base TetMesh Class for all types of Tetrahedron Mesh Classes
*
*		This is the fundamental class for tetrahedral meshes
*	   \author David Gu, Jerome Jiang
*      \date 10/01/2011
*	   \modified 4/16/2015
*
*/

#ifndef _HMESHLIB_BASE_HEX_MESH_H_
#define _HMESHLIB_BASE_HEX_MESH_H_

#include <assert.h>
#include <list>
#include <math.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <map>

#include "../Geometry/Point.h"
#include "../Geometry/Point2.h"
#include "../Parser/StrUtil.h"

#ifndef MAX_LINE 
#define MAX_LINE 2048
#endif


namespace HMeshLib{

	/*!
	* \brief CBaseTMesh, base class for all types of tet-mesh classes
	*
	*  This is the fundamental class for tet-meshes. All the geometric objects are connected by pointers,
	*  vertex, edge, face, tet are connected by halffaces. The mesh class has file IO functionalities,
	*  supporting .tet file formats. It offers Euler operators, each geometric primative
	*  can access its neighbors freely.
	*
	* \tparam CVertex   vertex   class, derived from TMeshLib::CVertex     class
	* \tparam CHVertex  tetrahedron vertex   class, derived from TMeshLib::CHVertex   class
	* \tparam CHalfEdge halfedge class, derived from TMeshLib::CHalfEdge class
	* \tparam CHEdge	tetrahedron edge class, derived from TMeshLib::CHEdge class
	* \tparam CEdge     edge     class, derived from MeshLib::CEdge     class
	* \tparam CFace     face     class, derived from TMeshLib::CFace     class
	* \tparam CHalfFace half face     class, derived from TMeshLib::CHalfFace     class
	* \tparam CHex      tetrahedron class, derived from TMeshLib::CHex class
	*/

	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	class CHMesh
	{
	public:
		/*!
			CHMesh constructor
			*/
		CHMesh(){};
		/*!
			CHMesh desctructor
			*/
		~CHMesh(){ _clear(); };
		/*!
			Load tet mesh from a ".tet" file
			*/
		void _load(const char *);
		/*!
			Load tet mesh from a ".t" file
			*/
		void _load_t(const char *);
		/*
			Load Hex mesh from a ".Qhex" file
			*/
		void _load_Qhex(const char*);
		/*!
			Write tet mesh to a file
			*/
		void _write(const char *);

		/*!
		Write tet mesh to a .t file
		*/
		void _write_Qhex(const char * filename);

		/*!
			access the list of half faces
			*/
		std::list<CHalfFace*> & half_faces(){ return m_pHalfFaces; };
		/*!
			access the list of edges
			*/
		std::list<CEdge*>       & edges()    { return m_pEdges; };
		/*!
			access list of faces
			*/
		std::list<CFace*>        & faces()    { return m_pFaces; };

		/*!
		access list of vertices
		*/
		std::list<CVertex*> & vertices() { return m_pVertices; };
		/*! Access the array of half edges */
		std::list<CHalfEdge*>& half_edges() { return m_pHalfEdges; };
		/*! Access the array of HEdges */
		std::list<CHEdge*>& hedges() { return m_pHEdges };
		/*! Access the array of Hexs */
		std::list<CHex*>& Hexs() { return m_pHexs; };

		/*! number of Hexs */
		int numHexs()       { return m_nHexs; };

		/*! number of vertices */
		int numVertices() { return m_nVertices; };

		/*! max vertex id*/
		int maxVertexId();

		

		/*! access the vertex with ID */
		CVertex * idVertex(int id) { return m_map_Vertices[id]; };

		/*! access the tet with ID */
		CHex      * idHex(int id)      { return m_map_Hexs[id]; };
		
		/*access the map vertex id*/
		std::map<int, CVertex*>& mapVertices() { return  m_map_Vertices; };
		
		/*access the map hex id*/
		std::map<int, CHex*>& mapHexs() { return  m_map_Hexs; };
			   

		//Access Vertex data members
		/*! Vertex->Edge List */
		std::list<CEdge*> * VertexEdgeList(CVertex * pVertex);
		/*! Vertex->TEdge List */
		std::list<CHEdge*> * VertexHEdgeList(CVertex * pVertex);
		/*! Vertex->HalfFace List */
		std::list<CHalfFace*> * VertexHalfFaceList(CVertex * pVertex);
		/*! Vertex->TVertex List */
		std::list<CHVertex*> * VertexHVertexList(CVertex * pVertex);
		/*! Vertex->Edge */
		CEdge   * VertexEdge(CVertex * v1, CVertex * v2);

		//Access HVertex data memebers
		CVertex * HVertexVertex(CHVertex * pTVertex);
		CHex * HVertexHex(CHVertex * pTVertex);
		CHalfEdge * HVertexHalfEdge(CHVertex * pTVertex);

		//Access HEdge data memebers
		CHalfEdge * HEdgeLeftHalfEdge(CHEdge * pTEdge);
		CHalfEdge * HEdgeRightHalfEdge(CHEdge * pTEdge);
		CEdge *		HEdgeEdge(CHEdge * pTEdge);
		CHEdge *	HEdgeDualHEdge(CHEdge * pTEdge);
		CHex *		HEdgeHex(CHEdge * pTEdge);

		//Access HalfEdge data members
		/*! HalfEdge->source vertex */
		CVertex * HalfEdgeSource(CHalfEdge * pHalfEdge);
		/*! HalfEdge->target vertex */
		CVertex * HalfEdgeTarget(CHalfEdge * pHalfEdge);
		/*! HalfEdge->source tvertex */
		CHVertex * HalfEdgeHSource(CHalfEdge * pHalfEdge);
		/*! HalfEdge->target tvertex */
		CHVertex * HalfEdgeHTarget(CHalfEdge * pHalfEdge);
		/*! HalfEdge->dual halfedge */
		CHalfEdge * HalfEdgeDual(CHalfEdge * pHalfEdge);
		/*! HalfEdge->next HalfEdge */
		CHalfEdge * HalfEdgeNext(CHalfEdge * pHalfEdge);
		/*! HalfEdge->prev HalfEdge */
		CHalfEdge * HalfEdgePrev(CHalfEdge * pHalfEdge);
		/*! HalfEdge->TEdge TEdge */
		CHEdge * HalfEdgeHEdge(CHalfEdge * pHalfEdge);
		/*! HalfEdge->HalfFace */
		CHalfFace * HalfEdgeHalfFace(CHalfEdge * pHalfEdge);
		/* Face->HalfEdge*/
		CHalfEdge* FaceHalfEdge(CFace* pFace);

		//Access Edge data members
		/*! TEdge list of the edge */
		std::list<CHEdge*>* EdgeHEdgeList(CEdge * pEdge);
		/*! Edge->Vertex1 */
		CVertex* EdgeVertex1(CEdge * pEdge);
		/*! Edge->Vertex2 */
		CVertex* EdgeVertex2(CEdge * pEdge);
		/*! length of the edge*/
		double EdgeLength(CEdge * pEdge);

		//Access HalfFace data memebers
		/*! HalfFace->HalfEdge */
		CHalfEdge* HalfFaceHalfEdge(CHalfFace * pHalfFace);
		/*! HalfFace->face */
		CFace* HalfFaceFace(CHalfFace * pHalfFace);
		/*! HalfFace->Tet */
		CHex * HalfFaceHex(CHalfFace * pHalfFace);
		/*! HalfFace->dual half face */
		CHalfFace * HalfFaceDual(CHalfFace * pHalfFace);
		/*! HalfFace->parallel half face*/
		CHalfFace * HalfFaceParallel(CHalfFace * pHalfFace);

		//Face
		/*! access the left half face of a face */
		CHalfFace * FaceLeftHalfFace(CFace * pFace);
		/*! access the right half face of a face */
		CHalfFace * FaceRightHalfFace(CFace * pFace);

		//Tetrahedron

		/*! access the j-th half edge of a tet */
		CHalfFace* HexHalfFace(CHex * pT, int j);
		/*! access the j-th tvertex of a tet */
		CHVertex * HexHVertex(CHex * pT, int j);
		/*! access the j-th vertex of a tet */
		CVertex * HexVertex(CHex * pT, int j);
		



		/*!
		construct Hexahedron
		\tparam v array of vertex ids
		\tparam pT retulting Hexahedron
		*/

		void  _construct_hex(CHex* pT, int tID, int * v);
		/*! construct faces */
		void  _construct_faces();
		/*! construct edges */
		void  _construct_edges();
		/*!  construct half face
		\tparam array of vertex pointers
		*/
		CHalfFace*   _construct_half_face(CHVertex **);

		/*! release all the memory allocations */
		void _clear();

	protected:
		/*!
		list of faces
		*/
		std::list<CFace*>        m_pFaces;
		/*!
		list of half faces
		*/
		std::list<CHalfFace*>	 m_pHalfFaces;
		/*!
		list of half edges
		*/
		std::list<CHalfEdge*>	 m_pHalfEdges;
		/*!
		list of edges
		*/
		std::list<CEdge*>        m_pEdges;
		/*!
		list of tetrahedra
		*/
		std::list<CHEdge*>		 m_pHEdges;

		/*!
		 array of vertices
		 */
		std::list<CVertex*>		 m_pVertices;
		//CVertex *				 m_pVertices;

		/*!
		map of CVertex id and pointer
		*/
		std::map<int, CVertex *> m_map_Vertices;

		/*!
		array of Hexs
		*/
		std::list<CHex*>		 m_pHexs;
		//CHex*                    m_pHexs;

		std::map<int, CHex*>     m_map_Hexs;

		/*! number of vertices */
		int m_nVertices;

		/*! number of Hexs */
		int m_nHexs;

		/*! number of edges */
		int m_nEdges;

		/*! num of sharp edges*/
		int m_nSharpEdges;

		/*! num of corner*/
		int m_nCorners;

		/*! max vertex id */
		int m_maxVertexId;

	};

	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	void CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::_clear()
	{
		for (std::list<CEdge*>::iterator eit = m_pEdges.begin(); eit != m_pEdges.end(); eit++)
		{
			CEdge * pE = *eit;
			delete pE;
		}

		for (std::list<CHEdge*>::iterator eit = m_pHEdges.begin(); eit != m_pHEdges.end(); eit++)
		{
			CHEdge * pE = *eit;
			delete pE;
		}

		for (std::list<CHalfEdge*>::iterator hit = m_pHalfEdges.begin(); hit != m_pHalfEdges.end(); hit++)
		{
			CHalfEdge * pH = *hit;
			delete pH;
		}


		for (std::list<CFace*>::iterator fit = m_pFaces.begin(); fit != m_pFaces.end(); fit++)
		{
			CFace * pF = *fit;
			delete pF;
		}

		for (std::list<CHalfFace*>::iterator fit = m_pHalfFaces.begin(); fit != m_pHalfFaces.end(); fit++)
		{
			CHalfFace * pF = *fit;
			delete pF;
		}

		m_pVertices.clear();
		m_pHexs.clear();
		//delete[] m_phvertices;
	};

	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	void CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::_load(const char * input)
	{
		m_maxVertexId = 0;

		std::fstream is(input, std::fstream::in);

		if (is.fail())
		{
			fprintf(stderr, "Error in opening file %s\n", input);
			return;
		}

		char buffer[MAX_LINE];
		{
			//read in the first line
			is.getline(buffer, MAX_LINE);
			//read in the number of vertices: "number vertices"
			std::string line(buffer);
			line = strutil::trim(line);
			strutil::Tokenizer stokenizer(line, " \r\n");
			stokenizer.nextToken();
			std::string token = stokenizer.getToken();
			m_nVertices = strutil::parseString<int>(token);
			stokenizer.nextToken();
			token = stokenizer.getToken();
			if (token != "vertices")
			{
				fprintf(stderr, "Error in file format %s\n", input);
				return;
			}
		}

		//read in the second line
		is.getline(buffer, MAX_LINE);
		//read in the number of vertices: "number Hexs"
		std::string line(buffer);
		line = strutil::trim(line);
		strutil::Tokenizer stokenizer(line, " \r\n");
		stokenizer.nextToken();
		std::string token = stokenizer.getToken();
		m_nHexs = strutil::parseString<int>(token);
		stokenizer.nextToken();
		token = stokenizer.getToken();
		if (token != "Hexs")
		{
			fprintf(stderr, "Error in file format %s\n", input);
			return;
		}

		//read in the vertices
		for (int i = 0; i < m_nVertices && is.getline(buffer, MAX_LINE); i++)
		{
			std::string line(buffer);
			line = strutil::trim(line);
			strutil::Tokenizer stokenizer(line, " \r\n");

			CPoint p;
			for (int k = 0; k < 3; k++)
			{
				stokenizer.nextToken();
				std::string token = stokenizer.getToken();
				p[k] = strutil::parseString<float>(token);
			}

			CVertex * v = new CVertex();
			//CVertex & v = m_pVertices[i];
			v->id() = i;
			v->position() = p;

			m_pVertices.push_back(v);
			m_map_Vertices.insert(std::pair<int, CVertex*>(i, v));

			if (!stokenizer.nextToken("\t\r\n")) continue;
			std::string token = stokenizer.getToken();

			int sp = (int)token.find("{");
			int ep = (int)token.find("}");

			if (sp >= 0 && ep >= 0)
			{
				v->string() = token.substr(sp + 1, ep - sp - 1);
			}
		}


		//read in Hexs 
		for (int id = 0; id < m_nHexs && is.getline(buffer, MAX_LINE); id++)
		{
			int vId[4];

			std::string line(buffer);
			line = strutil::trim(line);
			strutil::Tokenizer stokenizer(line, " \r\n");

			//skip the first "4" in the line
			stokenizer.nextToken();
			std::string token = stokenizer.getToken();

			for (int k = 0; k < 4; k++)
			{
				stokenizer.nextToken();
				std::string token = stokenizer.getToken();
				vId[k] = strutil::parseString<int>(token);
			}

			CHex * pT = new CHex();
			m_pHexs.push_back(pT);
			m_map_Hexs.insert(std::pair<int, CHex *>(id, pT));

			_construct_hex(pT, id, vId);
		}

		_construct_faces();
		_construct_edges();

		is.close();

		for (std::list<CVertex*>::iterator vIter = m_pVertices.begin(); vIter != m_pVertices.end(); vIter++)
		{
			CVertex * pV = *vIter;
			if (pV->id() > m_maxVertexId)
			{
				m_maxVertexId = pV->id();
			}
		}

	};


	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	void CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::_load_t(const char * input)
	{
		m_maxVertexId = 0;

		std::fstream is(input, std::fstream::in);

		if (is.fail())
		{
			fprintf(stderr, "Error in opening file %s\n", input);
			return;
		}

		char buffer[MAX_LINE];

		m_nVertices = 0;
		m_nHexs = 0;
		m_nEdges = 0;

		while (!is.eof())
		{
			is.getline(buffer, MAX_LINE);
			std::string line(buffer);
			line = strutil::trim(line);
			strutil::Tokenizer stokenizer(line, " \r\n");

			stokenizer.nextToken();
			std::string token = stokenizer.getToken();

			if (token == "Vertex") m_nVertices++;
			if (token == "Tet") m_nHexs++;
			if (token == "Edge") m_nEdges++;
		}

		is.clear();              // forget we hit the end of file
		is.seekg(0, std::ios::beg);   // move to the start of the file

		//read in the vertices
		for (int i = 0; i < m_nVertices && is.getline(buffer, MAX_LINE); i++)
		{
			std::string line(buffer);
			line = strutil::trim(line);
			strutil::Tokenizer stokenizer(line, " \r\n");

			stokenizer.nextToken();
			std::string token = stokenizer.getToken();

			if (token != "Vertex")
			{
				fprintf(stderr, "File Format Error\r\n");
				return;
			}


			stokenizer.nextToken();
			token = stokenizer.getToken();
			int vid = strutil::parseString<int>(token);

			CPoint p;
			for (int k = 0; k < 3; k++)
			{
				stokenizer.nextToken();
				std::string token = stokenizer.getToken();
				p[k] = strutil::parseString<float>(token);
			}

			CVertex * v = new CVertex();
			//CVertex & v = m_pVertices[i];
			v->id() = vid;
			v->position() = p;

			m_pVertices.push_back(v);

			m_map_Vertices.insert(std::pair<int, CVertex *>(vid, v));

			if (!stokenizer.nextToken("\t\r\n")) continue;
			token = stokenizer.getToken();

			int sp = (int)token.find("{");
			int ep = (int)token.find("}");

			if (sp >= 0 && ep >= 0)
			{
				v->string() = token.substr(sp + 1, ep - sp - 1);
			}
		}


		//read in Hexs 
		for (int id = 0; id < m_nHexs && is.getline(buffer, MAX_LINE); id++)
		{
			int vId[4];

			std::string line(buffer);
			line = strutil::trim(line);
			strutil::Tokenizer stokenizer(line, " \r\n");

			stokenizer.nextToken();
			std::string token = stokenizer.getToken();

			if (token != "Tet")
			{
				fprintf(stderr, "File Format Error\r\n");
				return;
			}

			//skip the first "4" in the line
			stokenizer.nextToken();
			token = stokenizer.getToken();
			int tid = strutil::parseString<int>(token);


			for (int k = 0; k < 4; k++)
			{
				stokenizer.nextToken();
				std::string token = stokenizer.getToken();
				vId[k] = strutil::parseString<int>(token);
			}



			CHex * pT = new CHex();
			m_pHexs.push_back(pT);
			m_map_Hexs.insert(std::pair<int, CHex *>(tid, pT));

			_construct_hex(pT, tid, vId);

			// read in string
			if (!stokenizer.nextToken("\t\r\n")) continue;
			token = stokenizer.getToken();

			int sp = (int)token.find("{");
			int ep = (int)token.find("}");

			if (sp >= 0 && ep >= 0)
			{
				pT->string() = token.substr(sp + 1, ep - sp - 1);
			}
		}

		_construct_faces();
		_construct_edges();

		for (int id = 0; id < m_nEdges && is.getline(buffer, MAX_LINE); id++)
		{
			std::string line(buffer);
			line = strutil::trim(line);
			strutil::Tokenizer stokenizer(line, " \r\n");

			stokenizer.nextToken();
			std::string token = stokenizer.getToken();

			if (token != "Edge")
			{
				fprintf(stderr, "File Format Error\r\n");
				return;
			}

			stokenizer.nextToken();
			token = stokenizer.getToken();
			int id1 = strutil::parseString<int>(token);

			stokenizer.nextToken();
			token = stokenizer.getToken();
			int id2 = strutil::parseString<int>(token);

			CVertex * pV1 = idVertex(id1);
			CVertex * pV2 = idVertex(id2);

			CEdge * pE = VertexEdge(pV1, pV2);

			if (!stokenizer.nextToken("\t\r\n"))
			{
				continue;
			}

			token = stokenizer.getToken();

			int sp = (int)token.find("{");
			int ep = (int)token.find("}");

			if (sp >= 0 && ep >= 0)
			{
				pE->string() = token.substr(sp + 1, ep - sp - 1);
			}
		}

		m_nEdges = (int)m_pEdges.size();

		is.close();

		for (std::list<CVertex*>::iterator vIter = m_pVertices.begin(); vIter != m_pVertices.end(); vIter++)
		{
			CVertex * pV = *vIter;
			if (pV->id() > m_maxVertexId)
			{
				m_maxVertexId = pV->id();
			}
		}

	};
	
	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	void CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::_load_Qhex(const char * input)
	{
		m_maxVertexId = 0;

		std::fstream is(input, std::fstream::in);

		if (is.fail())
		{
			fprintf(stderr, "Error in opening file %s\n", input);
			return;
		}

		char buffer[MAX_LINE];

		m_nVertices = 0;
		m_nHexs = 0;
		m_nEdges = 0;
		m_nSharpEdges = 0;
		m_nCorners = 0;

		while (!is.eof())
		{
			is.getline(buffer, MAX_LINE);
			std::string line(buffer);
			line = strutil::trim(line);
			strutil::Tokenizer stokenizer(line, " \r\n");

			stokenizer.nextToken();
			std::string token = stokenizer.getToken();

			if (token == "Vertex") m_nVertices++;
			if (token == "Hex") m_nHexs++;
			if (token == "Edge") m_nSharpEdges++;
			if (token == "corner") m_nCorners++;
			
			
		}

		is.clear();              // forget we hit the end of file
		is.seekg(0, std::ios::beg);   // move to the start of the file

		//read in the vertices
		for (int i = 0; i < m_nVertices && is.getline(buffer, MAX_LINE); i++)
		{
			std::string line(buffer);
			line = strutil::trim(line);
			strutil::Tokenizer stokenizer(line, " \r\n");

			stokenizer.nextToken();
			std::string token = stokenizer.getToken();

			if (token != "Vertex")
			{
				fprintf(stderr, "File Format Error\r\n");
				return;
			}


			stokenizer.nextToken();
			token = stokenizer.getToken();
			int vid = strutil::parseString<int>(token);

			CPoint p;
			for (int k = 0; k < 3; k++)
			{
				stokenizer.nextToken();
				std::string token = stokenizer.getToken();
				p[k] = strutil::parseString<float>(token);
			}

			CVertex * v = new CVertex();
			v->id() = vid;
			v->position() = p;

			m_pVertices.push_back(v);

			m_map_Vertices.insert(std::pair<int, CVertex *>(vid, v));

			if (!stokenizer.nextToken("\t\r\n")) continue;
			token = stokenizer.getToken();

			int sp = (int)token.find("{");
			int ep = (int)token.find("}");

			if (sp >= 0 && ep >= 0)
			{
				v->string() = token.substr(sp + 1, ep - sp - 1);
			}
		}


		//read in Hex 
		for (int id = 0; id < m_nHexs && is.getline(buffer, MAX_LINE); id++)
		{
			int vId[8];

			std::string line(buffer);
			line = strutil::trim(line);
			strutil::Tokenizer stokenizer(line, " \r\n");

			stokenizer.nextToken();
			std::string token = stokenizer.getToken();

			if (token != "Hex")
			{
				fprintf(stderr, "File Format Error\r\n");
				return;
			}

			//skip the first "4" in the line
			stokenizer.nextToken();
			token = stokenizer.getToken();
			int tid = strutil::parseString<int>(token);


			for (int k = 0; k < 8; k++)
			{
				stokenizer.nextToken();
				std::string token = stokenizer.getToken();
				vId[k] = strutil::parseString<int>(token);
			}



			CHex * pHex = new MyH();
			m_pHexs.push_back(pHex);
			m_map_Hexs.insert(std::pair<int, CHex *>(tid, pHex));

			_construct_hex(pHex, tid, vId);

			// read in string
			if (!stokenizer.nextToken("\t\r\n")) continue;
			token = stokenizer.getToken();

			int sp = (int)token.find("{");
			int ep = (int)token.find("}");

			if (sp >= 0 && ep >= 0)
			{
				pHex->string() = token.substr(sp + 1, ep - sp - 1);
			}
		}
		
		_construct_faces();
		_construct_edges();
		m_nEdges = (int)m_pEdges.size();

		//read in sharp edge
		for (int id = 0; id < m_nSharpEdges && is.getline(buffer, MAX_LINE); id++)
		{
			std::string line(buffer);
			line = strutil::trim(line);
			strutil::Tokenizer stokenizer(line, " \r\n");

			stokenizer.nextToken();
			std::string token = stokenizer.getToken();

			if (token != "Edge")
			{
				fprintf(stderr, "File Format Error\r\n");
				return;
			}
			stokenizer.nextToken();
			token = stokenizer.getToken();
			int id0 = strutil::parseString<int>(token);

			stokenizer.nextToken();
			token = stokenizer.getToken();
			int id1 = strutil::parseString<int>(token);

			CVertex* v0 = idVertex(id0);
			CVertex* v1 = idVertex(id1);

			CEdge* edge = vertexEdge(v0, v1);

			if (!stokenizer.nextToken("\t\r\n")) continue;
			token = stokenizer.getToken();

			int sp = (int)token.find("{");
			int ep = (int)token.find("}");

			if (sp >= 0 && ep >= 0)
			{
				edge->string() = token.substr(sp + 1, ep - sp - 1);
			}
		}

		/*read in corner*/
		for (int id = 0; id < m_nCorners && is.getline(buffer, MAX_LINE); id++)
		{
			std::string line(buffer);
			line = strutil::trim(line);
			strutil::Tokenizer stokenizer(line, " \r\n");

			stokenizer.nextToken();
			std::string token = stokenizer.getToken();
			if (token != "Corner")
			{
				fprintf(stderr, "File Format Error\r\n");
				return;
			}
			stokenizer.nextToken();
			token = stokenizer.getToken();
			int id = strutil::parseString<int>(token);
			CVertex* v = idVertex(id0);

			if (!stokenizer.nextToken("\t\r\n")) continue;
			token = stokenizer.getToken();

			int sp = (int)token.find("{");
			int ep = (int)token.find("}");

			if (sp >= 0 && ep >= 0)
			{
				v->string() = token.substr(sp + 1, ep - sp - 1);
			}
		}

		is.close();
		
		for (std::list<CVertex*>::iterator vIter = m_pVertices.begin(); vIter != m_pVertices.end(); vIter++)
		{
			CVertex * pV = *vIter;
			if (pV->id() > m_maxVertexId)
			{
				m_maxVertexId = pV->id();
			}
		}
	}
	
	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	CHalfFace* CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::_construct_half_face(CHVertex ** pTV)
	{
		CHalfFace * pHF = new CHalfFace;
		assert(pHF != NULL);
		m_pHalfFaces.push_back(pHF);

		CVertex * pV[4];

		for (int i = 0; i < 4; i++)
		{
			pV[i] = HVertexVertex(pTV[i]);
		}

		CHalfEdge * pH[4];
		for (int i = 0; i < 4; i++)
		{
			pH[i] = new CHalfEdge;
			assert(pH[i] != NULL);
			m_pHalfEdges.push_back(pH[i]);

			pH[i]->SetHalfFace(pHF);
			pH[i]->SetSource(pTV[i]);
			pH[i]->SetTarget(pTV[(i + 1) % 4]);
			pTV[i]->set_halfedge(pH[i]);
		}

		for (int i = 0; i < 4; i++)
		{
			pH[i]->SetNext(pH[(i + 1) % 4]);
			pH[i]->SetPrev(pH[(i + 3) % 4]);
		}

		pHF->SetHalfEdge(pH[0]);

		for (int i = 0; i < 4; i++)
		{						
			pHF->key(i) = pV[i]->id();				
		}
		
		//bubble

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3 - i; j++)
			{
				if (pHF->key(j) > pHF->key(j + 1))
				{
					int tmp = pHF->key(j);
					pHF->key(j) = pHF->key(j + 1);
					pHF->key(j + 1) = tmp;
				}
			}
		}

		assert(pHF->key(0) <= pHF->key(1));
		assert(pHF->key(1) <= pHF->key(2));
		assert(pHF->key(2) <= pHF->key(3));
		CVertex * keyv = m_map_Vertices[pHF->key(0)];

		keyv->HalfFaces()->push_back(pHF);

		return pHF;
	};

	//construct faces
	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	void CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::_construct_faces()
	{
		CVertex * pV = NULL;
		
		for (std::list<CVertex*>::iterator vIter = m_pVertices.begin(); vIter != m_pVertices.end(); vIter++)
		{
			pV = *vIter;
			
			std::list<CHalfFace*> * pL = VertexHalfFaceList(pV);	
			while (!pL->empty())
			{
				CHalfFace * pF = pL->front();
				pL->remove(pF);
				CFace * f = new CFace;
				assert(f != NULL);
				m_pFaces.push_back(f);
				f->SetLeft(pF);
				pF->SetFace(f);
				
				for (std::list<CHalfFace*>::iterator it = pL->begin(); it != pL->end(); it++)
				{
					CHalfFace * pH = *it;

					if (*pH == *pF)
					{
						pH->SetDual(pF);
						pF->SetDual(pH);
						f->SetRight(pH);
						pH->SetFace(f);
						break;
					}
				}

				if (pF->dual() != NULL)
				{
					
					pL->remove(HalfFaceDual(pF));
				}
				else
				{
					f->boundary() = true;
				}
			}
			

		}
		
	};

	//construct edges
	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	void CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::_construct_edges()
	{
		for (std::list<CVertex*>::iterator vIter = m_pVertices.begin(); vIter != m_pVertices.end(); vIter++)
		{
			CVertex * pV = *vIter;
			std::list<CHEdge*> * pL = VertexHEdgeList(pV);

			while (!pL->empty())
			{
				CHEdge * pTE = pL->front();
				pL->remove(pTE);
				CEdge * e = new CEdge;
				assert(e != NULL);

				int id1 = pTE->key(0);
				CVertex * v1 = m_map_Vertices[id1];
				e->SetVertex1(v1);

				int id2 = pTE->key(1);
				CVertex * v2 = m_map_Vertices[id2];
				e->SetVertex2(v2);

				m_pEdges.push_back(e);

				e->edges()->push_back(pTE);
				pTE->SetEdge(e);

				std::list<CHEdge*> tmp_edges;
				for (std::list<CHEdge*>::iterator it = pL->begin(); it != pL->end(); it++)
				{
					CHEdge * pH = *it;

					if (*pH == *pTE)
					{
						pH->SetEdge(e);
						tmp_edges.push_back(pH);
					}
				}

				for (std::list<CHEdge*>::iterator it = tmp_edges.begin(); it != tmp_edges.end(); it++)
				{
					CHEdge * pH = *it;
					pL->remove(pH);
					e->edges()->push_back(pH);
				}

			}

		}
		for (std::list<CEdge*>::iterator it = m_pEdges.begin(); it != m_pEdges.end(); it++)
		{
			CEdge * pE = *it;
			CVertex * v1 = EdgeVertex1(pE);
			CVertex * v2 = EdgeVertex2(pE);
			v1->edges()->push_back(pE);
			v2->edges()->push_back(pE);
		}
	};


	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	void CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::_construct_hex(CHex* pT, int tId, int * v)
	{
		//set the hex->id

		pT->id() = tId;

		//set HVertices of the Hex

		for (int k = 0; k < 8; k++)
		{
			CHVertex * pTV = new CHVertex();
			pT->setHVertex(pTV, k);
			pTV->id() = k;

			CVertex * pV = m_map_Vertices[v[k]];
			if (pV==NULL)
			{
				std::cout << "can't find v:" << std::endl;
				std::cout << "V->id(): " << v[k] << std::endl;
			}
			pTV->set_vert(pV);
			pV->hvertices()->push_back(pTV);

			pTV->set_hex(pT);
		}

		//set half faces
		int order[6][4] = { { 3,2,1,0 }, { 0,1,5,4 }, { 4,5,6,7 }, { 7,6,2,3 },{1,2,6,5},{3,0,4,7} };

		CHVertex   * pTV[4];
		CHalfFace * pHF[6];

		for (int i = 0; i < 6; i++)
		{
			for (int k = 0; k < 4; k++)
			{
				pTV[k] = HexHVertex(pT, order[i][k]);
			}
			pT->setHalfFace(_construct_half_face(pTV), i);
			pHF[i] = HexHalfFace(pT, i);
		}

		// connect the six half faces

		for (int i = 0; i < 6; i++)
		{
			pHF[i]->SetHex(pT);
		}

		//Seting the dual half edges
		for (int i = 0; i < 4; i++)
		{
			CHalfEdge* pH0 = HalfFaceHalfEdge(pHF[i]);
			CHalfEdge* pH1 = HalfEdgeNext(pH0);
			CHalfEdge* pH2 = HalfEdgeNext(pH1);
			CHalfEdge* pH3 = HalfEdgePrev(pH0);
			
			CHalfEdge * pH2_ = HalfFaceHalfEdge(pHF[(i + 1) % 4]);

			CHalfEdge * pH1_ = HalfFaceHalfEdge(pHF[4]);
			for (int j = 0; j < i; j++)
			{
				pH1_ = HalfEdgePrev(pH1_);
			}

			CHalfEdge * pH3_ = HalfFaceHalfEdge(pHF[5]);
			for (int k = 0; k < i; k++)
			{
				pH3_ = HalfEdgeNext(pH3_);
			}

			pH1->SetDual(pH1_);
			pH1_->SetDual(pH1);
			pH2->SetDual(pH2_);
			pH2_->SetDual(pH2);
			pH3->SetDual(pH3_);
			pH3_->SetDual(pH3);

			/*pH0->SetDual(pH1);
			pH1->SetDual(pH0);*/

			CHEdge * pTE1 = new CHEdge;
			CHEdge * pTE2 = new CHEdge;
			CHEdge * pTE3 = new CHEdge;

			assert(pTE1 != NULL);
			assert(pTE2 != NULL);
			assert(pTE3 != NULL);

			m_pHEdges.push_back(pTE1);
			m_pHEdges.push_back(pTE2);
			m_pHEdges.push_back(pTE3);
			pTE1->SetHex(pT);
			pTE2->SetHex(pT);
			pTE3->SetHex(pT);

			pH1->SetHEdge(pTE1);
			pH1_->SetHEdge(pTE1);
			pH2->SetHEdge(pTE2);
			pH2_->SetHEdge(pTE2);
			pH3->SetHEdge(pTE3);
			pH3_->SetHEdge(pTE3);



			if (pH1->source()->id() < pH1->target()->id())
			{
				pTE1->SetLeft(pH1);
				pTE1->SetRight(pH1_);
			}
			else
			{
				pTE1->SetLeft(pH1_);
				pTE1->SetRight(pH1);
			}

			if (pH2->source()->id() < pH2->target()->id())
			{
				pTE2->SetLeft(pH2);
				pTE2->SetRight(pH2_);
			}
			else
			{
				pTE2->SetLeft(pH2_);
				pTE2->SetRight(pH2);
			}

			if (pH3->source()->id() < pH3->target()->id())
			{
				pTE3->SetLeft(pH3);
				pTE3->SetRight(pH3_);
			}
			else
			{
				pTE3->SetLeft(pH3_);
				pTE3->SetRight(pH3);
			}

			pTE1->key(0) = pTE1->left()->source()->id();
			pTE1->key(1) = pTE1->left()->target()->id();

			pTE2->key(0) = pTE2->left()->source()->id();
			pTE2->key(1) = pTE2->left()->target()->id();

			pTE3->key(0) = pTE3->left()->source()->id();
			pTE3->key(1) = pTE3->left()->target()->id();

			CVertex * v1 = m_map_Vertices[pTE1->key(0)];
			CVertex * v2 = m_map_Vertices[pTE2->key(0)];
			CVertex * v3 = m_map_Vertices[pTE3->key(0)];

			v1->hedges()->push_back(pTE1);
			v2->hedges()->push_back(pTE2);
			v3->hedges()->push_back(pTE3);
		}

		int order_parallel[3][2] = { {0,2},{1,3},{4,5} };
		//setting halfface parallel
		for (int i = 0; i < 3; i++)
		{
			pHF[order_parallel[i][0]]->SetParallel(pHF[order_parallel[i][1]]);
			pHF[order_parallel[i][1]]->SetParallel(pHF[order_parallel[i][0]]);
		}
		//setting HEdge parallel
		for (int i = 0; i < 3; i++)
		{
			//find four parallel hedge each time
			std::vector<CHEdge*> v_hedge;
			CHalfEdge* he1 = HalfFaceHalfEdge(pHF[order_parallel[i][0]]);
			CHalfEdge* he2 = HalfFaceHalfEdge(pHF[order_parallel[i][1]]);
			CHEdge* hexedge1 = HalfEdgeHEdge(he1);
			CHEdge* hexedge2 = HalfEdgeHEdge(he2);
			if (i == 1)
			{
				he1 = HalfEdgeNext(he1);
				he2 = HalfEdgeNext(he2);
			}
			v_hedge.push_back(hexedge1);
			v_hedge.push_back(hexedge2);
			he1 = HalfEdgeNext(HalfEdgeNext(he1));
			he2 = HalfEdgeNext(HalfEdgeNext(he2));
			hexedge1 = HalfEdgeHEdge(he1);
			hexedge2 = HalfEdgeHEdge(he2);
			v_hedge.push_back(hexedge1);
			v_hedge.push_back(hexedge2);
			
			
			//setting parallel hedge
			for (std::vector<CHEdge*>::iterator ite1 = v_hedge.begin();ite1!=v_hedge.end();ite1++)
			{
				CHEdge* hedge1 = *ite1;
				for (std::vector<CHEdge*>::iterator ite2 = v_hedge.begin(); ite2 != v_hedge.end(); ite2++)
				{
					CHEdge* hedge2 = *ite2;
					if (hedge1 != hedge2)
					{
						hedge1->SetParallel(hedge2);
					}
				}
			}
		}
	};

	//write tet mesh to the file

	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	void CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::_write(const char * output)
	{

		std::fstream _os(output, std::fstream::out);
		if (_os.fail())
		{
			fprintf(stderr, "Error is opening file %s\n", output);
			return;
		}
		_os << m_nVertices << " vertices" << std::endl;
		_os << m_nHexs << " Hexs" << std::endl;

		for (std::list<CVertex*>::iterator vIter = m_pVertices.begin(); vIter != m_pVertices.end(); vIter++)
		{
			CVertex * pV = *vIter;
			CPoint  p = pV->position();
			for (int k = 0; k < 3; k++)
			{
				_os << " " << p[k];
			}
			if (pV->string().size() > 0)
			{
				_os << " " << "{" << pV->string() << "}";
			}
			_os << std::endl;
		}

		for (int i = 0; i < m_nHexs; i++)
		{
			CHex * pt = m_map_Hexs[i];
			_os << 4;
			for (int k = 0; k < 4; k++)
			{
				_os << " " << pt->tvertex(k)->vert()->id();
			}
			//if( pt->string().size() > 0 )
			//{
			//	_os << " " << "{"<< pt->string() << "}";
			//}
			_os << std::endl;
		}

		for (std::list<CEdge*>::iterator eIter = m_pEdges.begin(); eIter != m_pEdges.end(); eIter++)
		{
			CEdge * pE = *eIter;
			if (pE->string().size() > 0)
			{
				_os << "Edge " << EdgeVertex1(pE)->id() << " " << EdgeVertex2(pE)->id() << " ";
				_os << "{" << pE->string() << "}" << std::endl;
			}
		}

		_os.close();
	};

	template<typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	void CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::_write_Qhex(const char * output)
	{
		for (std::list<CVertex*>::iterator vIter = m_pVertices.begin(); vIter != m_pVertices.end(); vIter++)
		{
			CVertex * pV = *vIter;
			pV->_to_string();
		}

		for (std::list<CHex *>::iterator tIter = m_pHexs.begin(); tIter != m_pHexs.end(); tIter++)
		{
			CHex * pT = *tIter;
			pT->_to_string();
		}

		for (std::list<CEdge*>::iterator eIter = m_pEdges.begin(); eIter != m_pEdges.end(); eIter++)
		{
			CEdge * pE = *eIter;
			pE->_to_string();
		}
		//write traits end

		std::fstream _os(output, std::fstream::out);
		if (_os.fail())
		{
			fprintf(stderr, "Error while opening file %s\n", output);
			return;
		}

		for (std::list<CVertex*>::iterator vIter = m_pVertices.begin(); vIter != m_pVertices.end(); vIter++)
		{
			CVertex * pV = *vIter;
			CPoint p = pV->position();
			_os << "Vertex " << pV->id();
			for (int k = 0; k < 3; k++)
			{
				_os << " " << p[k];
			}
			if (pV->string().size() > 0)
			{
				_os << " " << "{" << pV->string() << "}";
			}
			_os << std::endl;
		}

		for (std::list<CHex *>::iterator tIter = m_pHexs.begin(); tIter != m_pHexs.end(); tIter++)
		{
			CHex * pT = *tIter;
			_os << "Hex " << pT->id();
			for (int k = 0; k < 8; k++)
			{
				_os << " " << pT->hvertex(k)->vert()->id();
			}
			if (pT->string().size() > 0)
			{
				_os << " " << "{" << pT->string() << "}";
			}
			_os << std::endl;
		}

		/*for (std::list<CEdge*>::iterator eIter = m_pEdges.begin(); eIter != m_pEdges.end(); eIter++)
		{
			CEdge * pE = *eIter;
			if (pE->string().size() > 0)
			{
				_os << "Edge " << EdgeVertex1(pE)->id() << " " << EdgeVertex2(pE)->id() << " ";
				_os << "{" << pE->string() << "}" << std::endl;
			}
		}*/

		_os.close();
	}

	/*------------------------------------------------------------------------------------------------
		Access TVertex data members
		--------------------------------------------------------------------------------------------------*/
	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CVertex * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::HVertexVertex(CHVertex * pTVertex)
	{
		return (CVertex*)pTVertex->vert();
	};
	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CHex * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::HVertexHex(CHVertex * pTVertex)
	{
		return (CHex*)pTVertex->hex();
	};

	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CHalfEdge * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::HVertexHalfEdge(CHVertex * pTVertex)
	{
		return (CHalfEdge*)pTVertex->halfedge();
	};
	/*------------------------------------------------------------------------------------------------
		Access Vertex data members
		--------------------------------------------------------------------------------------------------*/
	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline std::list<CEdge*> * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::VertexEdgeList(CVertex * pVertex)
	{
		return (std::list<CEdge*>*) pVertex->edges();
	};

	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline std::list<CHEdge*> * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::VertexHEdgeList(CVertex * pVertex)
	{
		return (std::list<CHEdge*>*) pVertex->hedges();
	};

	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline std::list<CHalfFace*> * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::VertexHalfFaceList(CVertex * pVertex)
	{
		return (std::list<CHalfFace*>*) pVertex->HalfFaces();
	};

	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline std::list<CHVertex*> * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::VertexHVertexList(CVertex * pVertex)
	{
		return (std::list<CHVertex*>*) pVertex->hvertices();
	};

	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CEdge * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::VertexEdge(CVertex * v1, CVertex * v2)
	{
		std::list<CEdge*> * vEdgeList = VertexEdgeList(v1);

		for (std::list<CEdge*>::iterator titer = (*vEdgeList).begin(); titer != (*vEdgeList).end(); titer++)
		{
			CEdge * pE = *titer;

			CVertex * w1 = EdgeVertex1(pE);
			CVertex * w2 = EdgeVertex2(pE);

			if (w1 == v1 && w2 == v2)
			{
				return pE;
			}
			if (w1 == v2 && w2 == v1)
			{
				return pE;
			}
		}

		return NULL;
	};

	/*------------------------------------------------------------------------------------------------
		Access HEdge data members
		--------------------------------------------------------------------------------------------------*/
	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CHalfEdge * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::HEdgeLeftHalfEdge(CHEdge * pTEdge)
	{
		return (CHalfEdge*)pTEdge->left();
	};
	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CHalfEdge * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::HEdgeRightHalfEdge(CHEdge * pTEdge)
	{
		return (CHalfEdge*)pTEdge->right();
	};
	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CEdge * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::HEdgeEdge(CHEdge * pTEdge)
	{
		return (CEdge*)pTEdge->edge();
	};
	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CHEdge * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::HEdgeDualHEdge(CHEdge * pTEdge)
	{
		return (CHalfEdge*)pTEdge->dual();
	};
	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CHex * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::HEdgeHex(CHEdge * pTEdge)
	{
		return (CHex*)pTEdge->Hex();
	};

	/*------------------------------------------------------------------------------------------------
		Access HalfEdge data members
		--------------------------------------------------------------------------------------------------*/
	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CVertex * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::HalfEdgeSource(CHalfEdge * pHalfEdge)
	{
		return (CVertex*)pHalfEdge->source();
	}
	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CVertex * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::HalfEdgeTarget(CHalfEdge * pHalfEdge)
	{
		return (CVertex*)pHalfEdge->target();
	}

	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CHVertex * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::HalfEdgeHSource(CHalfEdge * pHalfEdge)
	{
		return (CHVertex*)pHalfEdge->HSource();
	};

	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CHVertex * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::HalfEdgeHTarget(CHalfEdge * pHalfEdge)
	{
		return (CHVertex*)pHalfEdge->HTarget();
	};

	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CHalfEdge * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::HalfEdgeDual(CHalfEdge * pHalfEdge)
	{
		return (CHalfEdge*)pHalfEdge->dual();
	};

	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CHalfEdge * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::HalfEdgeNext(CHalfEdge * pHalfEdge)
	{
		return (CHalfEdge*)pHalfEdge->next();
	};

	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CHalfEdge * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::HalfEdgePrev(CHalfEdge * pHalfEdge)
	{
		return (CHalfEdge*)pHalfEdge->prev();
	};

	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CHEdge * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::HalfEdgeHEdge(CHalfEdge * pHalfEdge)
	{
		return (CHEdge*)pHalfEdge->Hedge();
	};

	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CHalfFace * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::HalfEdgeHalfFace(CHalfEdge * pHalfEdge)
	{
		return (CHalfFace*)pHalfEdge->half_face();
	}
	/*------------------------------------------------------------------------------------------------
		Access Edge data members
		--------------------------------------------------------------------------------------------------*/

	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline std::list<CHEdge*> * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::EdgeHEdgeList(CEdge * pEdge)
	{
		return (std::list<CHEdge*>*) pEdge->edges();
	};

	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CVertex * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::EdgeVertex1(CEdge * pEdge)
	{
		return (CVertex*)pEdge->vertex1();
	};

	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline double CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::EdgeLength(CEdge * pEdge)
	{
		CVertex * pV1 = EdgeVertex1(pEdge);
		CVertex * pV2 = EdgeVertex2(pEdge);

		return (pV1->position() - pV2->position()).norm();
	}

	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CVertex * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::EdgeVertex2(CEdge * pEdge)
	{
		return (CVertex*)pEdge->vertex2();
	}
	/*------------------------------------------------------------------------------------------------
		Access HalfFace data members
		--------------------------------------------------------------------------------------------------*/
	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CHalfEdge * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::HalfFaceHalfEdge(CHalfFace * pHalfFace)
	{
		return (CHalfEdge*)pHalfFace->half_edge();
	};

	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CFace * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::HalfFaceFace(CHalfFace * pHalfFace)
	{
		return (CFace*)pHalfFace->face();
	};

	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CHex * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::HalfFaceHex(CHalfFace * pHalfFace)
	{
		return (CHex*)pHalfFace->hex();
	};

	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CHalfFace * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::HalfFaceDual(CHalfFace * pHalfFace)
	{
		return (CHalfFace*)pHalfFace->dual();
	};
	
	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CHalfFace * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::HalfFaceParallel(CHalfFace * pHalfFace)
	{
		return (CHalfFace*)pHalfFace->parallel();
	};
	/*------------------------------------------------------------------------------------------------
		Access Face data members
		--------------------------------------------------------------------------------------------------*/
	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CHalfFace * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::FaceLeftHalfFace(CFace * pFace)
	{
		return (CHalfFace*)pFace->left();
	};

	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CHalfFace * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::FaceRightHalfFace(CFace * pFace)
	{
		return (CHalfFace*)pFace->right();
	};


	/*------------------------------------------------------------------------------------------------
		Access Tetrahedron data members
		--------------------------------------------------------------------------------------------------*/
	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CHalfEdge * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::FaceHalfEdge(CFace * pT)
	{
		return (CHalfEdge*)pT->left()->half_edge();
	};


	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CHalfFace * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::HexHalfFace(CHex * pT, int j)
	{
		return (CHalfFace*)pT->half_face(j);
	};

	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CHVertex * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::HexHVertex(CHex * pT, int j)
	{
		return (CHVertex*)pT->hvertex(j);
	};

	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	inline CVertex * CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::HexVertex(CHex * pT, int j)
	{
		return (CVertex*)pT->hvertex(j)->vert();
	};

	template <typename CHVertex, typename CVertex, typename CHalfEdge, typename CHEdge, typename CEdge, typename CHalfFace, typename CFace, typename CHex>
	int CHMesh<CHVertex, CVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex>::maxVertexId()
	{
		std::list<CVertex*>::iterator vite = m_pVertices.begin();
		for (vite; vite!= m_pVertices.end(); vite++)
		{
			CVertex* v = *vite;
			if (v->id()> m_maxVertexId)
			{
				m_maxVertexId = v->id();
			}
		}
		return m_maxVertexId;
	}
}

#endif _HMESHLIB_BASE_HEX_MESH_H_