#pragma once
#include "Mesh.h"
#include "ConjugateGradientSolver.h"
#include "SymSparseMatrix.h"
#include "PolarDecompose.h"


class Integrator
{
	Mesh* mesh;
	int n;
	float dt;
	bool **matmap;
	float **globalStiffness, ** globalMass, **globalDamping, **RK, **RKRT, **A;	
	float *extforces, *intforces, *x0, *xt, *fu, *b, *v, *mass;
	
	ConstrainedRows* rowSet;

public:
	Integrator(Mesh* _mesh, ConstrainedRows* r=0);
	void assembleLumpedMassVec();
	void assembleX0();
	void assembleExtForces();
	void assembleDampingMat();
	void assembleDisplacement();
	void assembleUndeformForces();
	void assembleRotations();
	void assembleA();
	void timeStep();
	void updateNodes();
	void debug();

	~Integrator(void);
};

