#pragma once
#include <math.h>
#include "core\HexMesh\hiterators.h"
#include "core\viewer\Arcball.h"                           /*  Arc Ball  Interface         */
#include "core/Geometry/Plane.h"
#include <time.h>
#include "topoMesh.h"
#include <queue>
namespace HMeshLib
{
	template<typename M>
	class Tools
	{
	public:
		typedef typename M::V V;
		typedef typename M::E E;
		typedef typename M::F F;
		typedef typename M::H H;

		Tools(M* input_mesh);
		~Tools();		
		void keyborad_save_files();	
		void load_feature();
		std::vector<E*> get_one_sheet(E* e);
		void get_sheet_hex();
	public:
		std::vector<E*> sharps;
		std::vector<E*> sheet;
		std::set<H*> sheet_hex;
		std::vector<E*> selectEs;
		std::vector<F*> inflate_fs;
	private:
		M* mesh;
		M original_mesh;
	};
	template<typename M>
	Tools<M>::Tools(M* input_mesh)
	{
		mesh = input_mesh;
	}
	template<typename M>
	Tools<M>::~Tools()
	{
	}

	template<typename M>
	void Tools<M>::keyborad_save_files()
	{
		std::string file_path;
		std::cout << "plaese input the file name" << std::endl;
		std::cin >> file_path;
		std::cout << "file name " << file_path << std::endl;
		mesh->write_Qhex(file_path.c_str());
	}

	template<typename M>
	void Tools<M>::load_feature()
	{
		for (M::MVIterator vite(mesh); !vite.end(); vite++)
		{
			V* v = *vite;
			v->_from_string();
			
		}
		for (M::MEIterator eite(mesh); !eite.end(); eite++)
		{
			E* e = *eite;
			e->_from_string();
			if (e->sharp())
			{
				sharps.push_back(e);
			}
			if (e->sheet())
			{
				sheet.push_back(e);
			}
		}
		for (M::MFIterator fite(mesh); !fite.end(); fite++)
		{
			F* f = *fite;
			f->_from_string();
			if (f->sheet_inflate())
			{
				inflate_fs.push_back(f);
			}
		}
		get_sheet_hex();
	}

	template<typename M>
	std::vector<typename M::E*> Tools<M>::get_one_sheet(E* e)
	{
		std::vector<E*> sheet;
		std::queue<E*> eQueue;
		eQueue.push(e);
		e->mark() = true;
		while (!eQueue.empty())
		{
			E* be = eQueue.front();
			sheet.push_back(be);
			eQueue.pop();
			for (int ehIndex = 0; ehIndex < be->neighbor_hs.size(); ehIndex++)
			{
				H* h = mesh->idHexs(be->neighbor_hs[ehIndex]);
				std::vector<E*> parallel_es = mesh->e_parallel_e_in_hex(h, be);
				for (int peIndex = 0; peIndex < parallel_es.size(); peIndex++)
				{
					E* pe = parallel_es[peIndex];
					if (pe->mark() == false)
					{
						eQueue.push(pe);
						pe->mark() = true;
					}
				}
			}
		}
		for (int eIndex = 0; eIndex < sheet.size(); eIndex++)
		{
			E* e = sheet[eIndex];
			e->mark() = false;
		}
		return sheet;
	}

	template<typename M>
	void Tools<M>::get_sheet_hex()
	{
		for (int eIndex = 0; eIndex < sheet.size(); eIndex++)
		{
			E* e = sheet[eIndex];
			for (int ehIndex = 0; ehIndex < e->neighbor_hs.size(); ehIndex++)
			{
				H* h = mesh->idHexs(e->neighbor_hs[ehIndex]);
				sheet_hex.insert(h);
			}
		}

	}
}
