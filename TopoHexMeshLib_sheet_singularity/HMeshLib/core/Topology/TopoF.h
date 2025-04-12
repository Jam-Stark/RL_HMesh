#ifndef TOPOLOGY_FACE_H
#define TOPOLOGY_FACE_H
#include<vector>
namespace HMeshLib
{
	class topoF
	{
	public:
		topoF() { m_id = 0; m_boundary = false; m_singularity_insert = false; m_sheet_id = 0; m_column_id = 0; };
		~topoF() {};
		
		int& id() { return m_id; };
		bool& boundary() { return m_boundary; };
		int& key(int k) { return m_key[k]; };
		CPoint& normal() { return m_normal; };
		int& sheetId() { return m_sheet_id; };
		int& columnId() { return m_column_id; };
		std::string& string() { return m_string; };
		virtual void _from_string() { };
		virtual void _to_string() { };

		std::vector<int> vs;
		std::vector<int> es;
		std::vector<int> neighbor_hs;
			
		int vertexIndex(int vid);
		int edgeIndex(int eid);
		int neighborHIndex(int hid);
		bool is_v_equal(std::vector<int> vid);
		bool delete_v(int vid);
		bool delete_e(int eid);
		bool delete_neighbor_h(int hid);
		void push_back_neighbor_h(int hid);

		/*table*/
		//the vertex adjcent vertex in this face
		const int vadjv[4][2] =
		{
			{3,1},
			{0,2},
			{1,3},
			{2,0}
		};
		const int eadje[4][2] =
		{
			{3,1},
			{0,2},
			{1,3},
			{2,0}
		};
		const int eparalle[4][1] =
		{
			{2},
			{3},
			{0},
			{1}
		};

		/*topology operation mark*/
		bool singularity_insert() { return m_singularity_insert; };

	protected:
		int m_id;
		bool m_boundary;
		int  m_key[4];//????????id??С???????У??????????????ж????
		CPoint m_normal;
		bool m_singularity_insert;
		int m_sheet_id;
		int m_column_id;
		std::string m_string;
	};

	/*	
	* get the index of the vertex of this face, if vertex is not belong to this face ,return -1
	*/
	int topoF::vertexIndex(int vid)
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
	* get the index of the edge of this face, if edge is not belong to this face ,return -1
	*/
	int topoF::edgeIndex(int eid)
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
	* get the index of the hex of neighbor hexs
	*/
	int topoF::neighborHIndex(int hid)
	{
		int result = -1;
		std::vector<int>::iterator ite = std::find(neighbor_hs.begin(), neighbor_hs.end(), hid);
		if (ite != neighbor_hs.end())
		{
			result = &*ite - &neighbor_hs[0];
		}
		return result;
	}

	bool topoF::is_v_equal(std::vector<int> vid)
	{
		for (int i = 0; i < 4; i++)
		{
			int Index = vertexIndex(vid[i]);
			if (Index == -1)
			{
				return false;
			}
		}
		return true;
	}

	/*delete element*/
	bool topoF::delete_v(int vid)
	{
		int fvIndex = vertexIndex(vid);
		if (fvIndex==-1)
		{
			return false;
		}
		else
		{
			vs.erase(vs.begin() + fvIndex);
			return true;
		}	
	}

	bool topoF::delete_e(int eid)
	{
		int feIndex = edgeIndex(eid);
		if (feIndex==-1)
		{
			return false;
		}
		else
		{
			es.erase(es.begin() + feIndex);
			return true;
		}
	}

	bool topoF::delete_neighbor_h(int hid)
	{
		int nhIndex = neighborHIndex(hid);
		if (nhIndex==-1)
		{
			return false;
		}
		else
		{
			neighbor_hs.erase(neighbor_hs.begin() + nhIndex);
			return true;
		}
	}

	void topoF::push_back_neighbor_h(int hid)
	{
		if (neighborHIndex(hid) == -1)
		{
			neighbor_hs.push_back(hid);
		}
	}
}

#endif // !TOPOLOGY_FACE_H