#ifndef TOPOLOGY_HEX_H
#define TOPOLOGY_HEX_H
#include<vector>
namespace HMeshLib
{
	class topoH
	{
	public:
		topoH() { m_id = 0; m_boundary = false; m_sheet_id = 0; m_column_id = 0; };
		~topoH() {};
		
		int& id() { return m_id; };
		bool& boundary() { return m_boundary; };
		int& sheetId() { return m_sheet_id; };
		int& columnId() { return m_column_id; };
		double& edge_angle(int eid) 
		{
			int index = edgeIndex(eid);
			assert(index != -1);
			return m_edge_angle[index];	
		};

		std::vector<int> vs;//点的存储按照输入的网格顺序来排序{0,1,2,3,4,5,6,7}
		std::vector<int> es;//边的排序{{3,2},{0,3},{1,2},{0,1},{0,4},{1,5},{4,5},{4,7},{5,6},{6,7},{3,7},{2,6}}
		std::vector<int> fs;//面的排序是{ { 3,2,1,0 }, { 0,1,5,4 }, { 4,5,6,7 }, { 7,6,2,3 },{1,2,6,5},{3,0,4,7} }	
		std::string& string() { return m_string; };
		virtual void _from_string() { };
		virtual void _to_string() { };

		/*table*/
		const int draw_order[6][4] = { { 3, 2, 1, 0 },{ 0,1,5,4 }, { 4,5,6,7 }, { 7,6,2,3 }, { 1,2,6,5 }, { 3,0,4,7 } };

		//  adjcent vertex of vertices in the hexahedron
		const int vadjv[8][3] =
		{
			{1,3,4},
			{0,2,5},
			{1,3,6},
			{0,2,7},
			{0,5,7},
			{1,4,6},
			{2,5,7},
			{3,4,6}
		};
		//adjcent edges of vertices in the hexahedron
		const int vadje[8][3] =
		{
			{1,3,4},
			{2,3,5},
			{0,2,11},
			{0,1,10},
			{4,6,7},
			{5,6,8},
			{9,8,11},
			{9,10,7}
		};
		//adjcent faces of vertices in the hexahedron
		const int vadjf[8][3] =
		{
			{0,1,5},
			{0,1,4},
			{0,3,4},
			{0,3,5},
			{1,2,5},
			{1,2,4},
			{2,3,4},
			{2,3,5}
		};
		//adjcent faces of edge in the hexahedron
		const int eadjf[12][2] =
		{
			{0,3},
			{0,5},
			{0,4},
			{0,1},
			{1,5},
			{1,4},
			{1,2},
			{2,5},
			{2,4},
			{2,3},
			{3,5},
			{3,4}
		};
		//edge adjcent edge in the hexahedron
		const int eadje[12][4] =
		{
			{1,2,10,11},
			{0,3,4,10},
			{0,3,5,11},
			{1,2,4,5},
			{1,3,6,7},
			{2,3,6,8},
			{4,7,5,8},
			{4,6,9,10},
			{5,6,9,11},
			{7,8,10,11},
			{0,1,7,9},
			{0,2,8,9}
		};
		//edge parallel edge in the hexahedron
		const int eparalle[12][3]
		{
			{3,6,9},//0
			{2,7,8},//1
			{1,7,8},//2
			{0,6,9},//3
			{5,10,11},//4
			{4,10,11},//5
			{0,3,9},//6
			{1,2,8},//7
			{1,2,7},//8
			{0,3,6},//9
			{4,5,11},//10
			{4,5,10},//11
		};
		//edge oppsite edge in hexahedron
		const int eoppsitee[12][1]
		{
			{6},
			{8},
			{7},
			{9},
			{11},
			{10},
			{0},
			{2},
			{1},
			{3},
			{5},
			{4}
		};
		//face adjcent face in the hexahedron
		const int fadjf[6][4] =
		{
			{1,3,4,5},
			{0,2,4,5},
			{1,3,4,5},
			{0,2,4,5},
			{0,1,2,3},
			{0,1,2,3}
		};
		// face parallel face in the hexahedron	
		const int fparallelf[6][1] =
		{
			{2},
			{3},
			{0},
			{1},
			{5},
			{4}
		};
		// face adjcent edges in the hexahedron
		const int fadje[6][4] =
		{
			{4,5,11,10},
			{1,2,8,7},
			{4,5,11,10},
			{1,2,8,7},
			{0,3,6,9},
			{0,3,6,9}
		};

		/*funcation*/
		int vertexIndex(int vid);
		int edgeIndex(int eid);
		int faceIndex(int fid);
		int fParallelF(int fid);
		

		/*delete element*/
		bool delete_v(int vid);
		bool delete_e(int eid);
		bool delete_f(int fid);
	protected:
		int m_id;
		bool m_boundary;
		int m_sheet_id;
		int m_column_id;
		double m_edge_angle[12];
		std::string m_string;
	};

	/*
	* get the index of the vertex of this hex, if vertex is not belong to this hex ,return -1
	*/
	int topoH::vertexIndex(int vid)
	{
		int result = -1;
		std::vector<int>::iterator ite = std::find(vs.begin(), vs.end(), vid);
		if (ite != vs.end())

		{
			result = &*ite - &vs[0];
		}
		return result;
	}

	/*
	* get the index of the edge of this hex, if edge is not belong to this hex ,return -1
	*/
	int topoH::edgeIndex(int eid)
	{
		int result = -1;
		std::vector<int>::iterator ite = std::find(es.begin(), es.end(), eid);
		if (ite != es.end())
		{
			result = &*ite - &es[0];
		}
		return result;
	}

	/*
	* get the index of the face of this hex, if face is not belong to this hex ,return -1
	*/
	int topoH::faceIndex(int fid)
	{
		int result = -1;
		std::vector<int>::iterator ite = std::find(fs.begin(), fs.end(), fid);
		if (ite != fs.end())
		{
			result = &*ite - &fs[0];
		}
		return result;
	}

	int topoH::fParallelF(int fid)
	{
		int fIndex = faceIndex(fid);
		if (fIndex == -1)
		{
			return fIndex;
		}
		int paraFIndex = fparallelf[fIndex][0];
		return fs[paraFIndex];

	}

	bool topoH::delete_v(int vid)
	{
		int vIndex = vertexIndex(vid);
		if (vIndex == -1)
		{
			return false;
		}
		else
		{
			vs.erase(vs.begin() + vIndex);
			return true;
		}
	}
	bool topoH::delete_e(int eid)
	{
		int eIndex = edgeIndex(eid);
		if (eIndex == -1)
		{
			return false;
		}
		else
		{
			es.erase(es.begin() + eIndex);
			return true;
		}
	}
	bool topoH::delete_f(int fid)
	{
		int fIndex = faceIndex(fid);
		if (fIndex == -1)
		{
			return false;
		}
		else
		{
			fs.erase(fs.begin() + fIndex);
			return true;
		}
	}
	

}

#endif // !TOPOLOGY_HEX_H