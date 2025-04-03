#ifndef TOPOLOGY_EDGE_H
#define TOPOLOGY_EDGE_H
#include<vector>
#include<string>
namespace HMeshLib
{
	class topoV;
	class topoF;
	class topoH;
	class topoE
	{
	public:
		topoE()
		{
			m_boundary = false;
			m_singularity = false;
			m_sheet_id = 0;
			m_column_id = 0;
			m_id = -1;
		}
		int& id() { return m_id; };
		bool& boundary() { return m_boundary; };
		bool& singularity() { return m_singularity; };
		int& sheetId() { return m_sheet_id; };
		int& columnId(){return m_column_id; };
		bool is_v_equal(std::vector<int> vid);	
		std::string& string() { return m_string; };
		virtual void _from_string() { };
		virtual void _to_string() { };

		std::vector<int> vs;//vertex id is sorted from low to high
		std::vector<int> neighbor_fs;
		std::vector<int> neighbor_hs;
		
		/*funcation*/
		int vertexIndex(int vid);
		int neighborFaceIndex(int fid);
		int neighborHexIndex(int hid);

		/*tool*/
		bool delete_v(int vid);
		bool delete_neighbor_f(int fid);
		bool delete_neighbor_h(int hid);
		
		void push_back_neighbor_f(int fid);
		void push_back_neighbor_h(int hid);

	protected:
		bool m_boundary;
		bool m_singularity;
		int m_sheet_id;
		int m_column_id;
		int m_id;
		std::string m_string;
	};

	int topoE::vertexIndex(int vid)
	{
		int result = -1;
		std::vector<int>::iterator ite = std::find(vs.begin(), vs.end(), vid);
		if (ite != vs.end())
		{
			result = &*ite - &vs[0];
		}
		return result;
	}
	int topoE::neighborFaceIndex(int fid)
	{
		int result = -1;
		std::vector<int>::iterator ite = std::find(neighbor_fs.begin(), neighbor_fs.end(), fid);
		if (ite != neighbor_fs.end())
		{
			result = &*ite - &neighbor_fs[0];
		}
		return result;
	}
	int topoE::neighborHexIndex(int hid)
	{
		int result = -1;
		std::vector<int>::iterator ite = std::find(neighbor_hs.begin(), neighbor_hs.end(), hid);
		if (ite != neighbor_hs.end())
		{
			result = &*ite - &neighbor_hs[0];
		}
		return result;
	}
	bool topoE::is_v_equal(std::vector<int> vid)
	{
		if (vid[0]==vid[1])
		{
			return false;
		}
		for (int i = 0; i < 2; i++)
		{
			int Index = vertexIndex(vid[i]);
			if (Index == -1)
			{
				return false;
			}
			
		}
		return true;
	}

	/*	
	* delete element
	*/
	bool topoE::delete_v(int vid)
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
	bool topoE::delete_neighbor_f(int fid)
	{
		int nfIndex = neighborFaceIndex(fid);
		if (nfIndex == -1)
		{
			return false;
		}
		else
		{
			neighbor_fs.erase(neighbor_fs.begin() + nfIndex);
			return true;
		}
	}
	bool topoE::delete_neighbor_h(int hid)
	{
		int nhIndex = neighborHexIndex(hid);
		if (nhIndex == -1)
		{
			return false;
		}
		else
		{
			neighbor_hs.erase(neighbor_hs.begin() + nhIndex);
			return true;
		}
	}

	void topoE::push_back_neighbor_f(int fid)
	{
		if (neighborFaceIndex(fid) == -1)
		{
			neighbor_fs.push_back(fid);
		}
	}
	void topoE::push_back_neighbor_h(int hid)
	{
		if (neighborHexIndex(hid) == -1)
		{
			neighbor_hs.push_back(hid);
		}
	}
}

#endif // !TOPOLOGY_EDGE_H
