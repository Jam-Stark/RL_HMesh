#ifndef _MESH_QUALITY_H_
#define _MESH_QUALITY_H_

using namespace std;

namespace HMeshLib
{
	template<typename M>
	class CMeshQuality
	{
	public:
		CMeshQuality(M* pMesh);
		~CMeshQuality();
	public:
		std::vector<CPoint> get_P(typename M::H* h);
		std::vector<CPoint> get_L(std::vector<CPoint> P);
		std::vector<double> get_edgeLengths(std::vector<CPoint> L);
		std::vector<CPoint> get_X(std::vector<CPoint> P);
		std::vector<CPoint> get_D(std::vector<CPoint> P);
		std::vector<std::vector<CPoint>> get_JacobianMatrices(std::vector<CPoint> L, std::vector<CPoint> X);
		std::vector<double> get_JacobianMatricesDet(typename M::H* h);
		std::vector<double> get_NormalizedJacobianMatricesDet(typename M::H* h);

	public:
		double Hex_Jacobian();
		double Hex_ScaledJacobian();
		double Hex_Diagonal();
		double Hex_EdgeRatio();

	protected:
		M* m_pMesh;		
	};

	template<typename M>
	CMeshQuality<M>::CMeshQuality(M* pMesh)
	{
		m_pMesh = pMesh;
			
	}

	template<typename M>
	CMeshQuality<M>::~CMeshQuality()
	{
	
	}	

	
	template<typename M>
	inline std::vector<CPoint> CMeshQuality<M>::get_P(typename M::H* h)
	{
		std::vector<CPoint> P;
		for (typename M::HVIterator hv(m_pMesh, h);!hv.end();hv++)
		{
			typename M::V* v = *hv;
			P.push_back(v->position());
		}
		return P;
	}

	template<typename M>
	inline std::vector<CPoint> CMeshQuality<M>::get_L(std::vector<CPoint> P)
	{
		std::vector<CPoint> L;
		L.push_back(P[1]-P[0]);
		L.push_back(P[2]-P[1]);
		L.push_back(P[3]-P[2]);
		L.push_back(P[3]-P[0]);
		L.push_back(P[4]-P[0]);
		L.push_back(P[5]-P[1]);
		L.push_back(P[6]-P[2]);
		L.push_back(P[7]-P[3]);
		L.push_back(P[5]-P[4]);
		L.push_back(P[6]-P[5]);
		L.push_back(P[7]-P[6]);
		L.push_back(P[7]-P[4]);

		return L;
	}

	template<typename M>
	inline std::vector<CPoint> CMeshQuality<M>::get_D(std::vector<CPoint> P)
	{
		std::vector<CPoint> D;
		D.push_back(P[6] - P[0]);
		D.push_back(P[7] - P[1]);
		D.push_back(P[4] - P[2]);
		D.push_back(P[5] - P[3]);
		return D;
	}

	template<typename M>
	std::vector<double> CMeshQuality<M>::get_edgeLengths(std::vector<CPoint> L)
	{
		std::vector<double> edgeLength;
		for (std::vector<CPoint>::iterator piter = L.begin(); piter != L.end(); piter++)
		{
			CPoint cP = *piter;
			edgeLength.push_back(cP.norm());
		}
		return edgeLength;
	}

	template<typename M>
	inline std::vector<CPoint> CMeshQuality<M>::get_X(std::vector<CPoint> P)
	{
		std::vector<CPoint> X;
		X.push_back(P[1] - P[0] + P[2] - P[3] + P[5] - P[4] + P[6] - P[7]);
		X.push_back(P[3] - P[0] + P[2] - P[1] + P[7] - P[4] + P[6] - P[5]);
		X.push_back(P[4] - P[0] + P[5] - P[1] + P[6] - P[2] + P[7] - P[3]);

		return X;
	}
	
	template<typename M>
	inline std::vector<std::vector<CPoint>> CMeshQuality<M>::get_JacobianMatrices(std::vector<CPoint> L, std::vector<CPoint> X)
	{
		std::vector<std::vector<CPoint>> A;
		std::vector<CPoint> a;
		a.push_back(L[0]);
		a.push_back(L[3]);
		a.push_back(L[4]);
		A.push_back(a);
		a.clear();
		a.push_back(L[1]);
		a.push_back(-L[0]);
		a.push_back(L[5]);
		A.push_back(a);
		a.clear();
		a.push_back(L[2]);
		a.push_back(-L[1]);
		a.push_back(L[6]);
		A.push_back(a);
		a.clear();
		a.push_back(-L[3]);
		a.push_back(-L[2]);
		a.push_back(L[7]);
		A.push_back(a);
		a.clear();
		a.push_back(L[11]);
		a.push_back(L[8]);
		a.push_back(-L[4]);
		A.push_back(a);
		a.clear();
		a.push_back(-L[8]);
		a.push_back(L[9]);
		a.push_back(-L[5]);
		A.push_back(a);
		a.clear();
		a.push_back(-L[9]);
		a.push_back(L[10]);
		a.push_back(-L[6]);
		A.push_back(a);
		a.clear();
		a.push_back(-L[10]);
		a.push_back(-L[11]);
		a.push_back(-L[7]);
		A.push_back(a);
		a.clear();
		a.push_back(X[0]);
		a.push_back(X[1]);
		a.push_back(X[2]);
		A.push_back(a);
		a.clear();
		return A;
	}

	template<typename M>
	inline std::vector<double> CMeshQuality<M>::get_JacobianMatricesDet(typename M::H* h)
	{
		std::vector<CPoint> P = get_P(h);
		std::vector<CPoint> L = get_L(P);
		std::vector<CPoint> X = get_X(P);
		std::vector<std::vector<CPoint>> A = get_JacobianMatrices(L,X);

		std::vector<double> a;
		for (int i = 0;i < A.size();i++)
		{
			a.push_back(A[i][0]*(A[i][1]^A[i][2]));
		}
		return a;
	}

	template<typename M>
	inline std::vector<double> CMeshQuality<M>::get_NormalizedJacobianMatricesDet(typename M::H* h)
	{
		std::vector<CPoint> P = get_P(h);
		std::vector<CPoint> L = get_L(P);
		std::vector<CPoint> X = get_X(P);
		std::vector<std::vector<CPoint>> A = get_JacobianMatrices(L, X);

		std::vector<double> a;
		for (int i = 0; i < A.size(); i++)
		{
			a.push_back((A[i][0]/ A[i][0].norm()) * ((A[i][1] / A[i][1].norm()) ^ (A[i][2] / A[i][2].norm())));
		}
		return a;
	}

	template<typename M>
	inline double CMeshQuality<M>::Hex_Jacobian()
	{
		double minJacobian = 1000000;
		double maxJacobian = 0;
		double sumJacobian = 0;
		int sumHex = 0;
		for (typename M::MHIterator hi(m_pMesh);!hi.end();hi++)
		{
			double singleHexJacobian = 1e10;
			typename M::H* h = *hi;
			std::vector<double> alphas= get_JacobianMatricesDet(h);
			alphas[8] /= 64.0;
			for (int i = 0;i < 9;i++)
			{
				if (singleHexJacobian > alphas[i])
				{
					singleHexJacobian = alphas[i];
				}
			}
			sumHex++;
			sumJacobian += singleHexJacobian;
			if (maxJacobian < singleHexJacobian)maxJacobian = singleHexJacobian;
			if (minJacobian > singleHexJacobian)minJacobian = singleHexJacobian;
			/*  singleHexJacobian为每一个体计算出来的Jacobian  */
			//cout << singleHexJacobian << endl;
		}
		double avgJacobian = sumJacobian / sumHex;

		//计算质量方差
		double VarJacobian = 0;
		for (typename M::MHIterator hi(m_pMesh); !hi.end(); hi++)
		{
			double qmin = 1e10;
			typename M::H* h = *hi;
			std::vector<double> alphas = get_JacobianMatricesDet(h);
			alphas[8] /= 64.0;

			for (int i = 0; i < 9; i++)
			{
				if (qmin > alphas[i])
				{
					qmin = alphas[i];
				}
			}
			VarJacobian += (qmin - avgJacobian) * (qmin - avgJacobian);
		}
		VarJacobian = VarJacobian / sumHex;
		cout << "-----------------------------Jacobian------------------------------" << endl;
		cout  << "Min:  " << minJacobian << "  Max:  " << maxJacobian << "  avg  " << avgJacobian <<"  Var  "<< VarJacobian<< endl;
		cout << endl;
		return 0.0;
	}

	template<typename M>
	inline double CMeshQuality<M>::Hex_ScaledJacobian()
	{
		double minScaledJacobian = 1000000;
		double maxScaledJacobian = 0;
		double sumScaledJacobian = 0;
		int sumHex = 0;
		for (typename M::MHIterator hi(m_pMesh); !hi.end(); hi++)
		{
			double singleHexScaledJacobian = 1e10;
			typename M::H* h = *hi;
			std::vector<double> alphas = get_NormalizedJacobianMatricesDet(h);
			for (int i = 0; i < 9; i++)
			{
				if (singleHexScaledJacobian > alphas[i])
				{
					singleHexScaledJacobian = alphas[i];
				}
			}
			std::vector<CPoint> P = get_P(h);
			std::vector<CPoint> L = get_L(P);
			std::vector<double> Len = get_edgeLengths(L);
			double minlen = 10000000;
			for (int i = 0; i < 12; i++)
			{
				if (Len[i] <  minlen)minlen = Len[i];
			}
			if (minlen * minlen <= DBL_MIN)singleHexScaledJacobian = DBL_MIN;
			sumHex++;
			sumScaledJacobian += singleHexScaledJacobian;
			if (maxScaledJacobian < singleHexScaledJacobian)maxScaledJacobian = singleHexScaledJacobian;
			if (minScaledJacobian > singleHexScaledJacobian)minScaledJacobian = singleHexScaledJacobian;
			/*  singleHexJacobian为每一个体计算出来的ScaledJacobian  */
			//cout << singleHexJacobian << endl;
		}
		double avgScaledJacobian = sumScaledJacobian / sumHex;

		//计算质量方差
		double VarScaledJacobian = 0;
		for (typename M::MHIterator hi(m_pMesh); !hi.end(); hi++)
		{
			double qmin = 1e10;
			typename M::H* h = *hi;
			std::vector<double> alphas = get_NormalizedJacobianMatricesDet(h);

			for (int i = 0; i < 9; i++)
			{
				if (qmin > alphas[i])
				{
					qmin = alphas[i];
				}
			}
			std::vector<CPoint> P = get_P(h);
			std::vector<CPoint> L = get_L(P);
			std::vector<double> Len = get_edgeLengths(L);
			double minlen = 10000000;
			for (int i = 0; i < 12; i++)
			{
				if (Len[i] < minlen)minlen = Len[i];
			}
			if (minlen * minlen <= DBL_MIN)qmin = DBL_MIN;

			VarScaledJacobian += (qmin - avgScaledJacobian) * (qmin - avgScaledJacobian);
		}
		VarScaledJacobian = VarScaledJacobian / sumHex;
		cout << "--------------------------ScaledJacobian---------------------------" << endl;
		cout  << "Min:  " << minScaledJacobian << "  Max:  " << maxScaledJacobian<< "  avg  " << avgScaledJacobian << "  Var  " << VarScaledJacobian << endl;
		cout << endl;
		return 0.0;
	}

	template<typename M>
	inline double CMeshQuality<M>::Hex_Diagonal()
	{
		double minDiagonal = 1000000;
		double maxDiagonal = 0;
		double sumDiagonal = 0;
		int sumHex = 0;
		for (typename M::MHIterator hi(m_pMesh); !hi.end(); hi++)
		{
			typename M::H* h = *hi;
			std::vector<CPoint> P = get_P(h);
			std::vector<CPoint> D = get_D(P);
			double minD = 1e10;
			double maxD = 0;
			for (int i = 0; i < 4; i++)
			{
				if (D[i].norm() > maxD)maxD = D[i].norm();
				if (D[i].norm() < minD)minD = D[i].norm();
			}
			double singleHexDiagonal = minD / maxD;
			if (maxD <= DBL_MIN)singleHexDiagonal = DBL_MIN;
			sumHex++;
			sumDiagonal += singleHexDiagonal;
			if (maxDiagonal < singleHexDiagonal)maxDiagonal = singleHexDiagonal;
			if (minDiagonal > singleHexDiagonal)minDiagonal = singleHexDiagonal;
			/*  singleHexDiagonal为每一个体计算出来的Diagonal  */
			//cout << singleHexDiagonal << endl;
		}
		double avgDiagonal = sumDiagonal / sumHex;

		//计算质量方差
		double VarDiagonal = 0;
		for (typename M::MHIterator hi(m_pMesh); !hi.end(); hi++)
		{
			typename M::H* h = *hi;
			std::vector<CPoint> P = get_P(h);
			std::vector<CPoint> D = get_D(P);
			double minD = 1e10;
			double maxD = 0;
			for (int i = 0; i < 4; i++)
			{
				if (D[i].norm() > maxD)maxD = D[i].norm();
				if (D[i].norm() < minD)minD = D[i].norm();
			}
			double singleHexDiagonal = minD / maxD;
			if (maxD <= DBL_MIN)singleHexDiagonal = DBL_MIN;

			VarDiagonal += (singleHexDiagonal - avgDiagonal) * (singleHexDiagonal - avgDiagonal);
		}
		VarDiagonal = VarDiagonal / sumHex;
		cout << "-----------------------------Diagonal------------------------------" << endl;
		cout << "Min:  " << minDiagonal << "  Max:  " << maxDiagonal << "  avg  " << avgDiagonal << "  Var  " << VarDiagonal << endl;
		cout << endl;
		return 0.0;
	}

	template<typename M>
	inline double CMeshQuality<M>::Hex_EdgeRatio()
	{
		double minEdgeRatio = 1000000;
		double maxEdgeRatio = 0;
		double sumEdgeRatio = 0;
		int sumHex = 0;
		for (typename M::MHIterator hi(m_pMesh); !hi.end(); hi++)
		{
			typename M::H* h = *hi;
			std::vector<CPoint> P = get_P(h);
			std::vector<CPoint> L = get_L(P);
			std::vector<double> Len = get_edgeLengths(L);
			double minlen = 10000000;
			double maxlen = 0;
			for (int i = 0; i < 12; i++)
			{
				if (Len[i] < minlen)minlen = Len[i];
				if (Len[i] > maxlen)maxlen = Len[i];
			}
			double singleHexEdgeRatio =maxlen/minlen ;
			sumHex++;
			sumEdgeRatio += singleHexEdgeRatio;
			if (maxEdgeRatio < singleHexEdgeRatio)maxEdgeRatio = singleHexEdgeRatio;
			if (minEdgeRatio > singleHexEdgeRatio)minEdgeRatio = singleHexEdgeRatio;
			/*  singleHexEdgeRatio为每一个体计算出来的ScaledJacobian  */
			//cout << singleHexEdgeRatio << endl;
		}
		double avgEdgeRatio = sumEdgeRatio / sumHex;

		//计算质量方差
		double VarEdgeRatio = 0;
		for (typename M::MHIterator hi(m_pMesh); !hi.end(); hi++)
		{
			typename M::H* h = *hi;
			std::vector<CPoint> P = get_P(h);
			std::vector<CPoint> L = get_L(P);
			std::vector<double> Len = get_edgeLengths(L);
			double minlen = 10000000;
			double maxlen = 0;
			for (int i = 0; i < 12; i++)
			{
				if (Len[i] < minlen)minlen = Len[i];
				if (Len[i] > maxlen)maxlen = Len[i];
			}
			double singleHexEdgeRatio = maxlen / minlen;

			VarEdgeRatio += (singleHexEdgeRatio - avgEdgeRatio) * (singleHexEdgeRatio - avgEdgeRatio);
		}
		VarEdgeRatio = VarEdgeRatio / sumHex;
		cout << "----------------------------Edge_Ratio-----------------------------" << endl;
		cout << "Min:  " << minEdgeRatio << "  Max:  " << maxEdgeRatio << "  avg  " << avgEdgeRatio << "  Var  " << VarEdgeRatio << endl;
		cout << endl;
		return 0.0;
	}
};
#endif