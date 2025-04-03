#ifndef TOPOLOGY_VERTEX_H
#define TOPOLOGY_VERTEX_H
#include<core/Geometry/Point.h>
namespace HMeshLib
{
	class topoE;
	class topoF;
	class topoH;

	class topoV
	{
	public:
		topoV() { m_boundary = false; m_singularity = false; m_sheet_id = 0; m_column_id = 0; };
		~topoV() {};

		int& id() { return m_id; };
		CPoint& position() { return m_position; };
		CPoint& normal() { return m_normal; };
		bool& singularity() { return m_singularity; };
		bool& singularity_node() { return m_singularity_node; };
		bool& boundary() { return m_boundary; };
		int& sheetId() { return m_sheet_id; };
		int& columnId() { return m_column_id; };


		std::vector<int> neighbor_es;
		std::vector<int> neighbor_fs;
		std::vector<int> neighbor_hs;
		std::string& string() { return m_string; };
		CPoint& pre_position() { return m_pre_position; };
		virtual void _from_string() { };
		virtual void _to_string() { };

		
		int neighborEdgeIndex(int eid);
		int neighborFaceIndex(int fid);
		int neighborHexIndex(int hid);

		
		bool delete_neighbor_e(int eid);
		bool delete_neighbor_f(int fid);
		bool delete_neighbor_h(int hid);

		void push_back_neighbor_e(int eid);
		void push_back_neighbor_f(int fid);
		void push_back_neighbor_h(int hid);

	protected:
		int m_id;
		CPoint m_position;
		bool m_boundary;
		bool m_singularity;
		bool m_singularity_node;
		int m_sheet_id;
		int m_column_id;
		std::string m_string;
		CPoint m_normal;
		CPoint m_pre_position;
	};
	
	/*
	* get the index of the edge of this neighbor edges
	*/
	int topoV::neighborEdgeIndex(int eid)
	{
		int result = -1;
		std::vector<int>::iterator ite = std::find(neighbor_es.begin(), neighbor_es.end(), eid);
		if (ite != neighbor_es.end())
		{
			result = &*ite - &neighbor_es[0];
		}
		return result;
	}
	/*
	* get the index of the face of this neighbor faces
	*/
	int topoV::neighborFaceIndex(int fid)
	{
		int result = -1;
		std::vector<int>::iterator ite = std::find(neighbor_fs.begin(), neighbor_fs.end(), fid);
		if (ite != neighbor_fs.end())
		{
			result = &*ite - &neighbor_fs[0];
		}
		return result;
	}
	/*
	* get the index of the hex of this neighbor hexs
	*/
	int topoV::neighborHexIndex(int hid)
	{
		int result = -1;
		std::vector<int>::iterator ite = std::find(neighbor_hs.begin(), neighbor_hs.end(), hid);
		if (ite != neighbor_hs.end())
		{
			result = &*ite - &neighbor_hs[0];
		}
		return result;
	}
	/*	
	* delete element
	*/
	
	bool topoV::delete_neighbor_e(int eid)
	{
		int neIndex = neighborEdgeIndex(eid);
		if (neIndex == -1)
		{
			return false;
		}
		else
		{
			neighbor_es.erase(neighbor_es.begin() + neIndex);
			return true;
		}
	}
	bool topoV::delete_neighbor_f(int fid)
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
	bool topoV::delete_neighbor_h(int hid)
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

	void topoV::push_back_neighbor_e(int eid)
	{
		if (neighborEdgeIndex(eid)==-1)
		{
			neighbor_es.push_back(eid);
		}
	}
	void topoV::push_back_neighbor_f(int fid)
	{
		if (neighborFaceIndex(fid) == -1)
		{
			neighbor_fs.push_back(fid);
		}
	}
	void topoV::push_back_neighbor_h(int hid)
	{
		if (neighborHexIndex(hid) == -1)
		{
			neighbor_hs.push_back(hid);
		}
	}

}

#endif // !TOPOLOGY_VERTEX_H