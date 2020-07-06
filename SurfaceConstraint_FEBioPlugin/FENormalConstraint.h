#pragma once

#include <FECore/FESurfaceConstraint.h>
#include <FECore/FESurface.h>
#include <iostream>							// Added variables for file writing [06/25/2020]
#include <fstream>							// Added variables for file writing [06/25/2020]


class FEVolumeSurface : public FESurface
{
public:
	//! constructor
	FEVolumeSurface(FEMesh* pm);

	//! Initialization
	bool Init();

	//! copy data
	void CopyFrom(FEVolumeSurface& s);

	//! serialization
	void Serialize(DumpStream& ar);

public:
	double Volume();

public:
	double	m_Lp;	//!< Lagrange multipler pressure
	double	m_p;	//!< applied pressure (= Lp + eps*DV)
	double	m_V0;	//!< Initial volume
	double	m_Vt;	//!< current volume
};

class FEFixedNormalDisplacement : public FESurfaceConstraint
{
public:
	// Constructor / destructor
	FEFixedNormalDisplacement(FEModel* pfem);
	~FEFixedNormalDisplacement();
	
	void Activate() override;
	void Residual(FEGlobalVector& R, const FETimeInfo& tp) override;
	void StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp) override;
	bool Augment(int naug, const FETimeInfo& tp) override { return true; }
	void Serialize(DumpStream& ar) override;
	void CopyFrom(FENLConstraint* plc) override;

	void Reset() override;
	void Update(int niter, const FETimeInfo& tp) override;

	void UnpackLM(FEElement& el, vector<int>& lm);
	
	void BuildMatrixProfile(FEGlobalMatrix& M) override;

	FESurface* GetSurface() override;

	// output log function
	virtual void write_variable(char* s);


public:
	FEVolumeSurface m_s;

public:
	double	m_eps;		//!< penalty parameter
	double	m_atol;		//!< augmented Lagrangian tolerance
	bool	m_blaugon;	//!< augmentation flag
	int		m_maxAug;	//!< maximum augmentation calculations
	int		m_minAug;	//!< minimum augmentation calculations

	// output variable log variables [06/25/2020]
	ofstream	var_log;	///!< The pointer to the log file
	char		line[100];
	double		oldTime;


private:
	bool	m_binit;	//!< flag indicating whether the constraint is initialized

	int m_dofX;
	int m_dofY;
	int m_dofZ;

	DECLARE_PARAMETER_LIST();
};