#pragma once

#include "ConjugateGradientSolver.h"
#include "SymSparseMatrix.h"
#include "PolarDecompose.h"

class Integrator
{
	QuadTetMesh* mesh;
	int n;
	float dt;
	bool **matmap;
	float **globalStiffness, ** globalMass, **globalDamping, **RK, **RKRT, **A;	
	float *extforces, *intforces, *x0, *xt, *fu, *b, *v, *mass, *kxt;
	bool *allowed;
	
	ConstrainedRows* rowSet;
	ConjugateGradientSolver solver;

public:
	Integrator(QuadTetMesh* _mesh, ConstrainedRows* r=0);
	
	//precomputation
	void assembleLumpedMassVec();
	void assembleX0();

	//timestep funcs
	void assembleExtForces();
	void assembleDisplacement();
	void assembleUndeformForces();
	void assembleRotations();
	void assembleA();

	//matfree timestep funcs
	void computeElementMatrices();
	void mulRK(float* in, float* out);
	void mulRKRT(float* in, float* out);
	void assembleKxt();

	void timeStep();
	void updateNodes();
	void debug();
	void sysMul(float* in, float* out);

	~Integrator(void);
};

