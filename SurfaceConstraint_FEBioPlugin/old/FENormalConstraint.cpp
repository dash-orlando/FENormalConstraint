

#include "stdafx.h"
#include "FENormalConstraint.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <FECore/FEDataExport.h>
using namespace std;

//-----------------------------------------------------------------------------
//! constructor
FEVolumeSurface::FEVolumeSurface(FEMesh* pm) : FESurface(pm)
{
	m_Lp = 0.0;
	m_p = 0.0;
	m_V0 = 0.0;
	m_Vt = 0.0;

	EXPORT_DATA(PLT_FLOAT, FMT_REGION, &m_p, "volume pressure");
}

// Checks that the new surface's volume is not 0
bool FEVolumeSurface::Init()
{
	m_V0 = Volume();	//Set initial volume
	m_Vt = m_V0;		//Set current volume

	return (m_V0 != 0.0);
}

void FEVolumeSurface::CopyFrom(FEVolumeSurface& s)
{
	m_Node = s.m_Node;		//all nodes in the surface? in the mesh?

	int NE = s.Elements();	// all the elements in the surface ? in the mesh ?
	Create(NE);				// creates the 
}

// serializing is saving, so saving/loading function?
void FEVolumeSurface::Serialize(DumpStream& ar)
{
	FESurface::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_Lp << m_p << m_V0 << m_Vt;
	}
	else
	{
		ar >> m_Lp >> m_p >> m_V0 >> m_Vt;
	}
}

// Calculate initial volume
double FEVolumeSurface::Volume()
{
	FEMesh& mesh = *GetMesh(); // get mesh pointer

	double vol = 0.0;
	int NE = Elements();
	vec3d x[FEElement::MAX_NODES];
	for (int i = 0; i < NE; ++i)
	{
		// get element
		FESurfaceElement& el = Element(i);

		// get nodal coordinates
		int neln = el.Nodes();
		for (int j = 0; j < neln; ++j)
		{
			x[j] = mesh.Node(el.m_node[j]).m_rt;
		}

		double* w = el.GaussWeights();
		int nint = el.GaussPoints(); // integration points
		for (int n = 0; n < nint; ++n)
		{
			//get position vector
			vec3d r = el.eval(x, n);

			// calculate tangent vector
			double* Gr = el.Gr(n);
			double* Gs = el.Gs(n);
			vec3d dxr(0, 0, 0), dxs(0, 0, 0);
			for (int j = 0; j < neln; ++j)
			{
				dxr += x[j] * Gr[j];
				dxs += x[j] * Gs[j];
			}

			vol += w[n] * (r*(dxr^dxs)); // update volume
		}
	}

	return vol / 3.0;
}

BEGIN_PARAMETER_LIST(FEFixedNormalDisplacement, FESurfaceConstraint);
	ADD_PARAMETER(m_blaugon, FE_PARAM_BOOL, "laugon");
	ADD_PARAMETER(m_atol, FE_PARAM_DOUBLE, "augtol");
	ADD_PARAMETER(m_eps, FE_PARAM_DOUBLE, "penalty");
	ADD_PARAMETER(m_maxAug, FE_PARAM_INT, "maxAug");
END_PARAMETER_LIST();

FEFixedNormalDisplacement::FEFixedNormalDisplacement(FEModel* pfem) : FESurfaceConstraint(pfem), m_s(&pfem->GetMesh())
{
	m_eps = 0.0;
	m_atol = 0.0;
	m_maxAug = 10;
	m_blaugon = false;
	m_binit = false;

	m_dofX = pfem->GetDOFIndex("x");
	m_dofY = pfem->GetDOFIndex("y");
	m_dofZ = pfem->GetDOFIndex("z");
}

void FEFixedNormalDisplacement::CopyFrom(FENLConstraint* plc)
{
	FEFixedNormalDisplacement& vc = dynamic_cast<FEFixedNormalDisplacement&>(*plc);

	GetParameterList() = vc.GetParameterList();

	m_s.CopyFrom(vc.m_s);
}

FESurface* FEFixedNormalDisplacement::GetSurface()
{
	return &m_s;
}

void FEFixedNormalDisplacement::Activate()
{
	FENLConstraint::Activate();

	if(m_binit == false)
	{
		m_s.Init();
	}

	m_binit = true;

}

//-----------------------------------------------------------------------------
void FEFixedNormalDisplacement::UnpackLM(FEElement& el, vector<int>& lm)
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	int N = el.Nodes();
	lm.resize(N * 3);
	for (int i = 0; i < N; ++i)
	{
		int n = el.m_node[i];
		FENode& node = mesh.Node(n);
		vector<int>& id = node.m_ID;

		lm[3 * i] = id[m_dofX];
		lm[3 * i + 1] = id[m_dofY];
		lm[3 * i + 2] = id[m_dofZ];
	}
}

//-----------------------------------------------------------------------------
void FEFixedNormalDisplacement::Residual(FEGlobalVector& R, const FETimeInfo& tp)
{
	FEMesh& mesh = *m_s.GetMesh();

	vector<double> fe;
	vector<int> lm;

	//double p = m_s.m_p;															// this is not needed here --06/20/2020--

	int NE = m_s.Elements();
	vec3d x[FEElement::MAX_NODES];
	vec3d un[FEElement::MAX_NODES];													// declaration of nodal dispalcement --06/18/2020--
	for (int i = 0; i < NE; ++i)
	{
		FESurfaceElement& el = m_s.Element(i);

		// get the nodal coordinates
		int neln = el.Nodes();
		for (int j = 0; j < neln; ++j)
		{
			x[j] = mesh.Node(el.m_node[j]).m_rt;
			un[j] = x[j] - mesh.Node(el.m_node[j]).m_r0;							// extracting nodal displacements --06/18/2020--
		}

		int ndof = 3 * neln;
		fe.resize(ndof);
		zero(fe);

		double* w = el.GaussWeights();
		int nint = el.GaussPoints();
		for (int k = 0; k < nint; ++k)												// replaced iteration counter 'n' for 'k' since normal vector 'n' is defined --06/18/2020--
		{
			double* Gr = el.Gr(k);
			double* Gs = el.Gs(k);
			vec3d dxr(0, 0, 0), dxs(0, 0, 0);
			vec3d u(0, 0, 0);														// defining vector 'u' --06/18/2020--
			double* H = el.H(k);													// mode definiton of 'H' prior to evaluation of vector 'u' --06/18/2020--
			for(int j = 0; j < neln; ++j)
			{
				dxr += x[j] * Gr[j];
				dxs += x[j] * Gs[j];
				u += un[j] * H[j];													// evaluation of u at the integration point --06/18/2020--
			}

			vec3d n = (dxr ^ dxs);
			n.unit();																// definition of 'n', the unit vector along 'v' --06/18/2020--
			//evaluate normal vector
			//vec3d v = (dxr ^ dxs)*w[n] * p;										// original to volume-constraint --06/18/2020--
			vec3d v = (dxr ^ dxs)*w[k]*m_eps*(u*n);									// modified evaluation of normal vector 'v', influenced by penalty factor 'm_eps' --06/18/2020--
			
			// evalueate element forces
			//double* H = el.H(k);													// original to volume-constraint --06/18/2020--
			for (int j = 0; j < neln; ++j)
			{
				/*
				fe[3 * j] += H[j] * v.x;
				fe[3 * j+1] += H[j] * v.y;
				fe[3 * j+2] += H[j] * v.z;
				*/

				fe[3 * j]		-= H[j] * v.x;										// element force needs to be subtracted from the global residual --06/20/2020--
				fe[3 * j + 1]	-= H[j] * v.y;
				fe[3 * j + 2]	-= H[j] * v.z;
			}
		}

		UnpackLM(el, lm);

		R.Assemble(el.m_node, lm, fe);
	}
}

//-----------------------------------------------------------------------------
void FEFixedNormalDisplacement::StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp)
{
	FEMesh& mesh = *m_s.GetMesh();

	double p = m_s.m_p;

	// element stiffness matrix
	matrix ke;
	vector<int> lm;
	vector<double> fe;

	int NE = m_s.Elements();
	vec3d x[FEElement::MAX_NODES];
	vec3d un[FEElement::MAX_NODES];													// declaration of nodal dispalcement --06/18/2020--
	for (int l = 0; l < NE; ++l)
	{
		FESurfaceElement& el = m_s.Element(1);

		// get nodal coordinates
		int neln = el.Nodes();
		for (int j = 0; j < neln; ++j)
		{
			x[j] = mesh.Node(el.m_node[j]).m_rt;
			un[j] = x[j] - mesh.Node(el.m_node[j]).m_r0;							// extracting nodal displacements --06/18/2020--
		}

		int ndof = 3 * neln;
		ke.resize(ndof, ndof);
		ke.zero();
		fe.resize(ndof);
		zero(fe);

		//repeat over integration points

		double* w = el.GaussWeights();
		int nint = el.GaussPoints();
		for (int k = 0; k < nint; ++k)												// replaced iteration counter 'n' for 'k' since normal vector 'n' is defined --06/18/2020--
		{
			double* N = el.H(k);
			double* Gr = el.Gr(k);
			double* Gs = el.Gs(k);
			//double* H = el.H(k);													// mode definiton of 'H' prior to evaluation of vector 'u' --06/18/2020--
			vec3d dxr(0, 0, 0), dxs(0, 0, 0);
			vec3d u(0, 0, 0);														// defining vector 'u' --06/18/2020--

			for (int j = 0; j < neln; ++j)
			{
				dxr += x[j] * Gr[j];
				dxs += x[j] * Gs[j];
				u += un[j] * N[j];													// evaluation of u at the integration point --06/18/2020-- is N the right one in this case?
			}

			// evaluate vector normal
			vec3d n = (dxr ^ dxs);
			n.unit();																// definition of 'n', the unit vector along 'v' --06/18/2020--
			//vec3d v = (dxr ^ dxs)*w[n] * m_eps;
			vec3d v = (dxr ^ dxs)*w[k] * m_eps*(u*n);								// modified evaluation of normal vector 'v', influenced by penalty factor 'm_eps' --06/18/2020--

			for (int i = 0; i < neln; ++i)
			{
				fe[3 * i] += N[i] * v.x;
				fe[3 * i + 1] += N[i] * v.y;
				fe[3 * i + 1] += N[i] * v.z;
			}
			for(int i=0; i<neln; ++i)
				for (int j = 0; j < neln; ++j)
				{
					ke[3 * i][3 * j] += fe[3 * i] * fe[3 * j];
					ke[3 * i][3 * j + 1] += fe[3 * i] * fe[3 * j + 1];
					ke[3 * i][3 * j + 2] += fe[3 * i] * fe[3 * j + 2];

					ke[3 * i + 1][3 * j] += fe[3 * i + 1] * fe[3 * j];
					ke[3 * i + 1][3 * j + 1] += fe[3 * i + 1] * fe[3 * j + 1];
					ke[3 * i + 1][3 * j + 2] += fe[3 * i + 1] * fe[3 * j + 2];

					ke[3 * i + 2][3 * j] += fe[3 * i + 2] * fe[3 * j];
					ke[3 * i + 2][3 * j + 1] += fe[3 * i + 2] * fe[3 * j + 1];
					ke[3 * i + 2][3 * j + 2] += fe[3 * i + 2] * fe[3 * j + 2];
				}

			// vec3d kab;															// for this code, we will replace the vec3d definition for a nested mat3d --06/20/2020--
			vec3d nu = dxr ^ dxs;													// cross-product of 'dxr' and 'dxs' [[note that this may be the same as 'n']] --06/20/2020--
			double da = nu.unit();													// magnitude of vector croos product 'nu' --06/20/2020--
			mat3ds Nu = dyad(nu);													// dyadic product of unit vector 'nu' with itself --06/20/2020--
			mat3ds ImN = mat3dd(1) - Nu;											// identity minus dyadic product 'N' (I-(nXn)) --06/20/2020--
			mat3d U = (mat3dd(1)*(u*nu) + (nu & u))*ImN;							// other.. need to compare to notes --06/20/2020-- [[why '&' instead of '^']]

			for(int i=0; i<neln; ++i)
				for (int j = 0; j < neln; ++j)
				{
					mat3d kab;														// 2nd-order tensor definition of Kab --06/20/2020--
					mat3da A(dxr*Gs[j] - dxs*Gr[j]);								// asymmetric tensor 'A' --06/20/2020--

					/*
					kab = (dxr*(N[j] * Gs[i] - N[i] * Gs[j])
						- dxs*(N[j] * Gr[i] - N[i] * Gr[j]))*w[k] * 0.5*p;

					ke[3 * i][3 * j] += 0;
					ke[3 * i][3 * j + 1] += -kab.z;
					ke[3 * i][3 * j + 2] += kab.y;

					ke[3 * i + 1][3 * j] += kab.z;
					ke[3 * i + 1][3 * j + 1] += 0;
					ke[3 * i + 1][3 * j + 2] += -kab.x;

					ke[3 * i + 2][3 * j] += -kab.y;
					ke[3 * i + 2][3 * j + 1] += kab.x;
					ke[3 * i + 2][3 * j + 2] += 0;
					*/

					kab = (U*A*N[j] + N[j]*N[i]*Nu*da)*w[k]*0.5*p*m_eps;			// new form of tensor 'Kab', including all new definitions

					ke[3 * i][3 * j] += kab(0, 0);
					ke[3 * i][3 * j + 1] += kab(0, 1);
					ke[3 * i][3 * j + 2] += kab(0, 2);

					ke[3 * i + 1][3 * j] += kab(1, 0);
					ke[3 * i + 1][3 * j + 1] += kab(1, 1);
					ke[3 * i + 1][3 * j + 2] += kab(1, 2);

					ke[3 * i + 2][3 * j] += kab(2, 0);
					ke[3 * i + 2][3 * j + 1] += kab(2, 1);
					ke[3 * i + 2][3 * j + 2] += kab(2, 2);

				}
		}

		UnpackLM(el, lm);

		psolver->AssembleStiffness(el.m_node, lm, ke);

	}
}

//-----------------------------------------------------------------------------
bool FEFixedNormalDisplacement::Augment(int naug, const FETimeInfo& tp)
{

	if (!m_blaugon || m_atol <= 0.0) return true;

	if (naug > m_maxAug) return true;

	double Dp = m_eps*(m_s.m_Vt - m_s.m_V0);
	double Lp = m_s.m_p;
	double err = fabs(Dp / Lp);

	if (err < m_atol) return true;

	m_s.m_Lp = Lp;
	m_s.m_p = Lp + Dp;

	// here we need to print the augmented tolerance value...

	return false;
}

void FEFixedNormalDisplacement::Serialize(DumpStream& ar)
{
	FENLConstraint::Serialize(ar);
	m_s.Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_binit;
	}
	else
	{
		ar >> m_binit;
	}
}

void FEFixedNormalDisplacement::Reset()
{

}

//-----------------------------------------------------------------------------
void FEFixedNormalDisplacement::BuildMatrixProfile(FEGlobalMatrix& M)
{
	// nothing to do here
}

void FEFixedNormalDisplacement::Update(int niter, const FETimeInfo& tp)
{
	m_s.m_Vt = m_s.Volume();

	m_s.m_p = m_s.m_Lp + m_eps*(m_s.m_Vt - m_s.m_V0);
}