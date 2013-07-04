#pragma once
#include "Mesh.h"
#include "ConstrainedRows.h"

class GPUIntegrator
{
	Mesh* mesh;
	float *extforces, *x0;
	int n;
public:
	GPUIntegrator(Mesh* _mesh, ConstrainedRows* r=0);
	void initVars();
	
	void timeStep();
	~GPUIntegrator(void);
};