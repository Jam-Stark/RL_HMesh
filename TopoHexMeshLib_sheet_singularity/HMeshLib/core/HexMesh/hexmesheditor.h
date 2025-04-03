#ifndef _HEX_MESH_EDITOR_
#define _HEX_MESH_EDITOR_

#include"basehmesh.h"

namespace HMeshLib
{
	template<typename HM>
	class HexMeshEditor
	{
	public:
		HexMeshEditor(HM* inputMesh) {hmesh = inputMesh};
		~HexMeshEditor() {};

		typedef typename HM::CVertex V;
		typedef typename HM::CHVertex HV;
		typedef typename HM::CHalfEdge HalfE;
		typedef typename HM::CHEdge HE;
		typedef typename HM::CHalfFace HF;
		typedef typename HM::CFace F;
		typedef typename HM::CHex Hex;

		void changeVertex(Hex* hex, V* v,V* changeV);
		void separateFace(F* f);/*seprate the face to two faces*/
		V* creatVertex(CPoint p, int id);
	protected:
		HM* hmesh = NULL;
	};
	
	template<typename HM>
	typename HexMeshEditor<typename HM>::V* HexMeshEditor<typename HM>::creatVertex(CPoint p, int id)
	{
		V* v = new V();
		assert(v != NULL);
		v->position() = p;
		v->id() = id;
		hmesh->m_pVertices().push_back(v);
		hmesh->m_map_Vertices.insert(std::pair<int, CVertex*>(id, v));

	}

	template<typename HM>
	void HexMeshEditor<typename HM>:: changeVertex(Hex* hex, V* v,V* changeV)
	{
		/*
		* 需要修改的内容：
		*	创建一个新的点
		*	找到所有旧的点关联的hex
		*		Hex内修改内容所有关联半边的点
		*		修改vertex中所有的Hvertices和 edges
		*		修改所有的关联的边的值
		*		
		*		修改所有关联的HVertex
		*	删除旧点	
		*/
		//construct a new vertex

		//如果只修改位置，不修改连接关系的话，只修改坐标即可
	}
	template<typename HM>
	void HexMeshEditor<typename HM>::separateFace(F* f)
	{
		/**将一个面分离为两个面*/
		/*
		* 将一个面分离为两个面，并插入一个六面体
		*	复制面上的四个点，将面上点的左右两边关系都分别赋值给两边。
		*	
		*/
		
	}
}

#endif // !_HEX_MESH_EDITOR_
