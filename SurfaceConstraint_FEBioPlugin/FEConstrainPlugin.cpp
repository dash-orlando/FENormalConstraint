// FEConstrainPlugin.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include <FECore\FECoreKernel.h>
#include "FENormalConstraint.h"


FECoreKernel* pFEBio;

//! Factory Classes for this plugin
FEPluginFactory_T<FEFixedNormalDisplacement, FENLCONSTRAINT_ID> normal_constraint_factory ("surface-normal-constraint");

//-----------------------------------------------------------------------------
FECORE_EXPORT unsigned int GetSDKVersion()
{
	return FE_SDK_VERSION;
}

//-----------------------------------------------------------------------------
FECORE_EXPORT void PluginInitialize(FECoreKernel& febio)
{
	FECoreKernel::SetInstance(&febio);
	pFEBio = &febio;
}

//-----------------------------------------------------------------------------
FECORE_EXPORT FECoreFactory* PluginGetFactory(int i)
{
	switch (i)
	{
	case 0: return &normal_constraint_factory;
	default: return 0;
	}
}

//-----------------------------------------------------------------------------
FECORE_EXPORT void PluginCleanup()
{

}

//-----------------------------------------------------------------------------
FECORE_EXPORT void GetPluginVersion(int& major, int& minor, int& patch)
{
	major = 0;
	minor = 1;
	patch = 0;
}