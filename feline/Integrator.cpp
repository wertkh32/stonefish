#include "Integrator.h"


Integrator::Integrator(Mesh* _mesh)
{
	mesh = _mesh;
	n = mesh->getNoOfNodes();
	//K is constant -> linear FEM. nonlinear FEM later.
	globalStiffness = mesh->assembleGlobalStiffness();
	globalMass = mesh->assembleGlobalMass();

	dt = 1./FPS;
	//printf("%lf ",dt);

	extforces = (double*)malloc(sizeof(double) * n * 3);
	x0 = (double*)malloc(sizeof(double) * n * 3);
	xt = (double*)malloc(sizeof(double) * n * 3);
	v = (double*)malloc(sizeof(double) * n * 3);
	fu = (double*)malloc(sizeof(double) * n * 3);
	b = (double*)malloc(sizeof(double) * n * 3);

	globalDamping = (double**)malloc(sizeof(double*) * n * 3);
	for(int i=0;i<n * 3;i++)
		globalDamping[i] = (double*)malloc(sizeof(double) * n * 3);

	A = (double**)malloc(sizeof(double*) * n * 3);
	for(int i=0;i<n * 3;i++)
		A[i] = (double*)malloc(sizeof(double) * n * 3);

	Rot = (double**)malloc(sizeof(double*) * n * 3);
	for(int i=0;i<n * 3;i++)
		Rot[i] = (double*)malloc(sizeof(double) * n * 3);

	solver = new ConjugateGradientSolver(n * 3, A);

	assembleDampingMat();
	assembleUndeformForces();
	assembleA();
}

void
Integrator::assembleUndeformForces()
{
	for(int i=0;i<n;i++)
	{
		x0[i * 3] = mesh->nodes[i]->pos.x;
		x0[i * 3 + 1] = mesh->nodes[i]->pos.y;
		x0[i * 3 + 2] = mesh->nodes[i]->pos.z;
	}

	for(int i=0;i<n *3;i++)
	{
		fu[i] = 0;
		for(int j=0;j<n*3;j++)
		{
			fu[i] += globalStiffness[i][j] *x0[j];
			printf("[%.0lf] ", globalStiffness[i][j]);
		}
	}
}

void
Integrator::assembleDampingMat()
{
	//damping mat. constant values for now
	double alpha = 0.5, beta = 0.5;

	for(int i=0;i<n*3;i++)
		for(int j=0;j<n*3;j++)
		{
			globalDamping[i][j] = globalStiffness[i][j] * alpha;// + globalMass[i][j] * beta;
			//printf("la %lf ", globalDamping[i][j]);
		}

}

void
Integrator::assembleRotations()
{
	for(int i=0;i<n;i++);
}

void
Integrator::assembleDisplacement()
{
	for(int i=0;i<n;i++)
	{
		xt[i * 3] = mesh->nodes[i]->pos_t.x;
		xt[i * 3 + 1] = mesh->nodes[i]->pos_t.y;
		xt[i * 3 + 2] = mesh->nodes[i]->pos_t.z;

		v[i * 3] = mesh->nodes[i]->vec_t.x;
		v[i * 3 + 1] = mesh->nodes[i]->vec_t.y;
		v[i * 3 + 2] = mesh->nodes[i]->vec_t.z;
	}
}

void
Integrator::assembleExtForces()
{
	//just gravity for now
		for(int i=0;i<n;i++)
	{
		extforces[i * 3] = mesh->nodes[i]->force.x;
		extforces[i * 3 + 1] = mesh->nodes[i]->force.y;
		extforces[i * 3 + 2] = mesh->nodes[i]->force.z;
	}

/*
	for(int i=0;i<n * 3; i++)
	//	if(i%3==1)
	//	{
	//		extforces[i] += GRAVITY * mesh->nodes[i]->mass;
	//	}
	//	else
		{
			extforces[i] = 0;
		}
		*/
}

void
Integrator::assembleA()
{
	for(int i=0;i<n*3;i++)
		for(int j=0;j<n*3;j++)
		{
			A[i][j] = globalMass[i][j] + globalDamping[i][j] * dt + globalStiffness[i][j] * dt * dt;
			//printf("%lf ",A[i][j]);
		}
}

void
Integrator::updateNodes()
{
	for(int i=0;i<n;i++)
	{
		vector3<double> temp(v[i * 3],v[i * 3 + 1],v[i * 3 + 2]);
		mesh->nodes[i]->pos_t += (temp * dt);
		mesh->nodes[i]->vec_t = temp;

		//printf("%lf, %lf, %lf\n", temp.x,temp.y,temp.z);
	}
}

void
Integrator::timeStep()
{
	//Av(t + 1) = b; (3.13)
	//where
	//A = (M + tC + t2K); 
	//b = Mv(t) - dt(Kx(t) + fu - fext)
	assembleExtForces();
	assembleDisplacement();

	for(int i=0;i<n*3;i++)
	{
		b[i] = 0;
		for(int j=0;j<n*3;j++)
		{
			b[i] += (globalMass[i][j] * v[j] - dt * (globalStiffness[i][j] * (xt[j] - x0[j])));
		}

		//b[i] += dt * (fu[i] - extforces[i]);
		b[i] = dt* extforces[i];
	}
	
	solver->solve(v,b);
	updateNodes();
}

Integrator::~Integrator(void)
{
}
