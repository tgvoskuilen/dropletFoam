#include "readTimeControls.H"

volSourceLimit =
    runTime.controlDict().lookupOrDefault<scalar>("volSourceLimit", 1e-2);

bool correctPhi = 
	pimple.dict().lookupOrDefault<Switch>("correctPhi", true);

bool checkMeshCourantNo = 
	pimple.dict().lookupOrDefault<Switch>("checkMeshCourantNo", false);
