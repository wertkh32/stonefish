#pragma once
#include "Mesh.h"
#include "QuadTetMesh.h"
#include "ConstrainedRows.h"
#include "GPUDataStructs.cuh"
#include "GPUIntegratorFuncs.cuh"


class GPUIntegrator
{
	MESH* mesh;
	int numnodes, numelements;
	int aligned_numnodes;

	GPUElement* gpuElements;
	GPUNode*   gpuNodes;

	float*   xt;//dynamic
	float*   vt;//dynamic
	float*	 extforces;//dynamic
	float*	 mass;
	char*	 allowed;

public:
	GPUIntegrator(MESH* _mesh, ConstrainedRows* r=0);
	
	void assembleGPUElements();
	void assembleGPUNodes();

	void assembleXt();
	void assembleVt();
	void assembleExtForce();

	void assembleLumpedMass();

	void initVars();
	void copyVarstoGPU();
	
	void timeStep();
	void updatePositions();

	~GPUIntegrator(void);
};