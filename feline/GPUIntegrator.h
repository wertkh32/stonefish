#pragma once
#include "Mesh.h"
#include "ConstrainedRows.h"
#include "GPUDataStructs.cuh"

class GPUIntegrator
{
	Mesh* mesh;
	int numnodes, numelements;

	GPUElement* gpuElements;
	GPUNode*   gpuNodes;

	float*   xt;//dynamic
	float*   vt;//dynamic
	float*	 extforces;//dynamic

public:
	GPUIntegrator(Mesh* _mesh, ConstrainedRows* r=0);
	
	void assembleGPUElements();
	void assembleGPUNodes();

	void assembleXt();
	void assembleVt();
	void assembleExtForce();

	void initVars();
	void copyVarstoGPU();
	
	void timeStep();
	~GPUIntegrator(void);
};