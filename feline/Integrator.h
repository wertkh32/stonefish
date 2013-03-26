#pragma once
#include "Mesh.h"
#include "ConjugateGradientSolver.h"

class Integrator
{
	Mesh* mesh;
	int n;
	float dt;
	float **globalStiffness, ** globalMass, **globalDamping, **A;	
	float *extforces, *intforces, *x0, *xt, *fu, *b, *v;
	ConjugateGradientSolver* solver;

public:
	Integrator(Mesh* _mesh);
	void assembleExtForces();
	void assembleDampingMat();
	void assembleDisplacement();
	void assembleUndeformForces();
	void assembleA();
	void timeStep();
	void updateNodes();

	~Integrator(void);
};

