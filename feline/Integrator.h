#pragma once
#include "Mesh.h"
#include "ConjugateGradientSolver.h"

class Integrator
{
	Mesh* mesh;
	int n;
	double dt;
	double **globalStiffness, ** globalMass, **globalDamping, **A;	
	double *extforces, *intforces, *x0, *xt, *fu, *b, *v;
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

