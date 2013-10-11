#include "Integrator.h"


Integrator::Integrator(MESH* _mesh, ConstrainedRows* r)
{
	rowSet = r;
	mesh = _mesh;
	n = mesh->getNoOfNodes();

	//K is constant -> linear FEM. nonlinear FEM later.
	//globalStiffness = mesh->assembleGlobalStiffness();
	//NON-LINEAR TIME!!!
	//globalMass = mesh->assembleGlobalMass();

	dt = 1./FPS;
	//printf("%lf ",dt);

	extforces = (float*)malloc(sizeof(float) * n * 3);
	x0 = (float*)malloc(sizeof(float) * n * 3);
	xt = (float*)malloc(sizeof(float) * n * 3);
	v = (float*)malloc(sizeof(float) * n * 3);
	fu = (float*)malloc(sizeof(float) * n * 3);
	b = (float*)malloc(sizeof(float) * n * 3);
	mass = (float*)malloc(sizeof(float) * n * 3);
	kxt = (float*)malloc(sizeof(float) * n * 3);

	allowed = (bool*)malloc(sizeof(bool) * n * 3);
	
	for(int i=0;i<n * 3;i++)
	{
		allowed[i] = true;
	}

	for(int i=0;i<rowSet->list.size();i++)
	{
		allowed[rowSet->list[i]] = false;
	}


	/*
	globalDamping = (float**)malloc(sizeof(float*) * n * 3);
	for(int i=0;i<n * 3;i++)
		globalDamping[i] = (float*)malloc(sizeof(float) * n * 3);
	*/


	/*
	A = (float**)malloc(sizeof(float*) * n * 3);
	for(int i=0;i<n * 3;i++)
	{
		A[i] = (float*)malloc(sizeof(float) * n * 3);
		for(int j=0;j<n*3;j++)
			A[i][j] = 0.0;
	}
	*/
	//RK = (float**)malloc(sizeof(float*) * n * 3);
	//for(int i=0;i<n * 3;i++)
	//	RK[i] = (float*)malloc(sizeof(float) * n * 3);

	//RKRT = (float**)malloc(sizeof(float*) * n * 3);
	//for(int i=0;i<n * 3;i++)
	//RKRT[i] = (float*)malloc(sizeof(float) * n * 3);

	solver.initSolver(n*3,A);
	assembleX0();
	//mass lumping
	assembleLumpedMassVec();



}

void
Integrator::assembleX0()
{
	for(int i=0;i<n;i++)
	{
		x0[i * 3] = mesh->nodes[i]->pos.x;
		x0[i * 3 + 1] = mesh->nodes[i]->pos.y;
		x0[i * 3 + 2] = mesh->nodes[i]->pos.z;
	}
}

void
Integrator::assembleLumpedMassVec()
{
	for(int i=0;i<n*3;i++)
		mass[i] = 0;

	for(int i=0;i<mesh->elements.size();i++)
	{
		//float elenodemass = (mesh->elements[i]->getDensity() * mesh->elements[i]->getVolume()) /4;
		for(int j=0;j<NUM_NODES_PER_ELE;j++)
		{
			mass[mesh->nodeIndices[i][j] * 3] += mesh->elements[i]->nodemass[j];
			mass[mesh->nodeIndices[i][j] * 3 + 1] += mesh->elements[i]->nodemass[j];
			mass[mesh->nodeIndices[i][j] * 3 + 2] += mesh->elements[i]->nodemass[j];
		}
	}
}

/*
void
Integrator::assembleDampingMat()
{
	//ignore
	//damping mat. constant values for now
	float alpha = 0.1, beta = 0.3;

	for(int i=0;i<n*3;i++)
		for(int j=0;j<n*3;j++)
		{
			globalDamping[i][j] = RKRT[i][j] * alpha + globalMass[i][j] * beta;
			//printf("la %lf ", globalDamping[i][j]);
			//printf("la %lf ", RKRT[i][j]);
		}

}
*/

void
Integrator::computeElementMatrices()
{
	for(int i=0;i<mesh->elements.size();i++)
	{
		//mesh->elements[i]->computeRKRTandRK();
		//mesh->elements[i]->computeMatFreeVars();
		mesh->elements[i]->computeRotation();
	}
}

void
Integrator::mulRK(float* in, float* out)
{

	for(int i=0;i<mesh->elements.size();i++)
		{
			ELEMENT& ele = *(mesh->elements[i]);
			GenMatrix<float,NUM_NODES_PER_ELE * 3,NUM_NODES_PER_ELE * 3>& A = ele.getStiffnessMat();
			Matrix3d& R = ele.getRotation();

			float pos[NUM_NODES_PER_ELE * 3] = {0};
			float pos2[NUM_NODES_PER_ELE * 3] = {0};

			for(int a=0;a<NUM_NODES_PER_ELE;a++)
			{
				int x = mesh->nodeIndices[i][a];

				pos[a * 3] = in[x * 3];
				pos[a * 3 + 1] = in[x * 3 + 1];
				pos[a * 3 + 2] = in[x * 3 + 2];
			}

			for(int a=0;a<NUM_NODES_PER_ELE * 3;a++)
				for(int b=0;b<NUM_NODES_PER_ELE * 3;b++)
					pos2[a] += A(a,b) * pos[b];


			for(int a=0;a<NUM_NODES_PER_ELE;a++)
			{
				pos[a*3] = 0;
				pos[a*3+1] = 0;
				pos[a*3+2] = 0;

				for(int b=0;b<3;b++)
					for(int c=0;c<3;c++)
						pos[a*3+b] += R(b,c) * pos2[a*3+c];
			}

			for(int a=0;a<NUM_NODES_PER_ELE;a++)
			{
				int x = mesh->nodeIndices[i][a];

				out[x * 3] += pos[a * 3];
				out[x * 3 + 1] += pos[a * 3 + 1];
				out[x * 3 + 2] += pos[a * 3 + 2];
			}
		}
}

void
Integrator::mulRKRT(float* in, float* out)
{
		for(int i=0;i<mesh->elements.size();i++)
		{
			ELEMENT& ele = *(mesh->elements[i]);
			GenMatrix<float,NUM_NODES_PER_ELE * 3,NUM_NODES_PER_ELE * 3>& A = (ele.getStiffnessMat());
			Matrix3d& R = ele.getRotation();

			float pos[NUM_NODES_PER_ELE * 3] = {0};
			float pos2[NUM_NODES_PER_ELE * 3] = {0};
			float pos3[NUM_NODES_PER_ELE * 3] = {0};

			for(int a=0;a<NUM_NODES_PER_ELE;a++)
			{
				int x = mesh->nodeIndices[i][a];

				pos[a * 3] = in[x * 3];
				pos[a * 3 + 1] = in[x * 3 + 1];
				pos[a * 3 + 2] = in[x * 3 + 2];
			}

			for(int a=0;a<NUM_NODES_PER_ELE;a++)
			{
				for(int b=0;b<3;b++)
					for(int c=0;c<3;c++)
						pos2[a*3+b] += R(c,b) * pos[a*3+c];
			}

			for(int a=0;a<NUM_NODES_PER_ELE * 3;a++)
				for(int b=0;b<NUM_NODES_PER_ELE * 3;b++)
					pos3[a] += A(a,b) * pos2[b];

			for(int a=0;a<NUM_NODES_PER_ELE;a++)
			{
				pos[a*3] = 0;
				pos[a*3+1] = 0;
				pos[a*3+2] = 0;

				for(int b=0;b<3;b++)
					for(int c=0;c<3;c++)
						pos[a*3+b] += R(b,c) * pos3[a*3+c];
			}

			for(int a=0;a<NUM_NODES_PER_ELE;a++)
			{
				int x = mesh->nodeIndices[i][a];

				out[x * 3] += pos[a * 3];
				out[x * 3 + 1] += pos[a * 3 + 1];
				out[x * 3 + 2] += pos[a * 3 + 2];
			}

		}
}


void
Integrator::assembleUndeformForces()
{
	
	for(int i=0;i<n *3;i++)
	{
		fu[i] = 0;
		//for(int j=0;j<n*3;j++)
		//{
		//	fu[i] += RK[i][j] * x0[j];
		//}
	}

	//mesh->mulRK(x0,fu);
	mulRK(x0,fu);
}

void
Integrator::assembleKxt()
{
	for(int i=0;i<n *3;i++)
	{
		kxt[i] = 0;
		//for(int j=0;j<n*3;j++)
		//{
		//	fu[i] += RK[i][j] * x0[j];
		//}
	}

	//mesh->mulRKRT(xt,kxt);
	mulRKRT(xt,kxt);
}

void
Integrator::assembleRotations()
{
	/*
	for(int i=0;i<n*3;i++)
		for(int j=0;j<n*3;j++)
		{
			RK[i][j] = 0.0;
			RKRT[i][j] = 0.0;
		}

	for(int i=0;i<mesh->elements.size();i++)
	{
		GenMatrix<float,12,12> *rk=0, *rkrt=0;
		mesh->elements[i]->getRKRTandRK(rk,rkrt);

		GenMatrix<float,12,12> &_rk = *rk, &_rkrt = *rkrt;

		for(int a=0;a<4;a++)
			for(int b=0;b<4;b++)
				for(int j=0;j<3;j++)
					for(int k=0;k<3;k++)
					{
						RK[mesh->nodeIndices[i][a] * 3 + j][mesh->nodeIndices[i][b] * 3 + k] += _rk(a * 3 + j, b * 3 + k);
						RKRT[mesh->nodeIndices[i][a] * 3 + j][mesh->nodeIndices[i][b] * 3 + k] += _rkrt(a * 3 + j, b * 3 + k);
					}
			
	}
	*/
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

	for(int i=0;i<n;i++)
	{
		extforces[i * 3] = mesh->nodes[i]->force.x;
		extforces[i * 3 + 1] = mesh->nodes[i]->force.y;
		extforces[i * 3 + 2] = mesh->nodes[i]->force.z;
	}

}

void
Integrator::assembleA()
{
	static const float alpha = 0.1, beta = 0.3;
	static const float coeffK = dt * beta + dt * dt, coeffM = 1 + dt * alpha;
	
	/*
	for(int i=0;i<n*3;i++)
		for(int j=0;j<n*3;j++)
		{
			{
				A[i][j] =  RKRT[i][j] * coeffK; //globalMass[i][j] * coeffM
				if(i==j)
					A[i][j] += mass[i] * coeffM;
			}
		}
	*/

	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
		{
			if(matmap[i][j])
			{
				for(int a=0;a<3;a++)
					for(int b=0;b<3;b++)
					{	
						int x = i*3+a, y = j*3+b;
						A[x][y] = RKRT[x][y] * coeffK;
						if(x==y)
							A[x][y] += mass[x] * coeffM;
					}
			}

		}

}

void
Integrator::updateNodes()
{
	for(int i=0;i<n;i++)
	{
		vector3<float> temp(v[i * 3],v[i * 3 + 1],v[i * 3 + 2]);
		mesh->nodes[i]->pos_t += (temp * dt);
		mesh->nodes[i]->vec_t = temp;

		//printf("%lf, %lf, %lf\n", temp.x,temp.y,temp.z);
	}
}

void
Integrator::sysMul(float* in, float* out)
{
	for(int i=0;i<n*3;i++)
		out[i] = 0;

	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
		{
			if(matmap[i][j])
			{
				for(int a=0;a<3;a++)
					for(int b=0;b<3;b++)
					{
						int x = i * 3 + a;
						int y = j * 3 + b;
						if(allowed[x] && allowed[y])
							out[x] += A[x][y] * in[y];
					}
			}
		}

}

void
Integrator::timeStep()
{
	assembleDisplacement();
	//assembleRotations();
	//assembleDampingMat();
	//assembleRK();
	computeElementMatrices();
	assembleUndeformForces();//fu
	assembleKxt();
	//assembleA();
	//Av(t + 1) = b; (3.13)
	//where
	//A = (M + tC + t2K); 
	//b = Mv(t) - dt(Kx(t) + fu - fext)
	assembleExtForces();


	//ConjugateGradientSolver solver(n*3,A);
	//printf("\nB\n");
	for(int i=0;i<n*3;i++)
	{
		//b[i] = 0;
		//for(int j=0;j<n*3;j++)
		//{
		//	b[i] += globalMass[i][j] * v[j];
		//}
		//dont sum separately. causes 0.0/-0.0 error in solver.
		b[i] = mass[i] * v[i] + dt * (fu[i] + extforces[i] - kxt[i]);

		//printf("%f\n", b[i]);
	}
	//system("pause");

	if(rowSet)
	{
		solver.solveWithConstraints(v,b,allowed,mesh);
	}
	else
	{
		solver.solve(v,b);
	}
	updateNodes();
}

void
Integrator::debug()
{
	printf("RKRT\n");
	for(int i=0;i<n*3;i++)
	{
		for(int j=0;j<n*3;j++)
		{
			printf("%.2lf ", RKRT[i][j]);
		}
		printf("\n");
	}
	printf("\n");

	printf("Mass\n");
	for(int i=0;i<n*3;i++)
	{
		for(int j=0;j<n*3;j++)
		{
			printf("%.2lf ", globalMass[i][j]);
		}
		printf("\n");
	}
	printf("\n");

	/*
	printf("Stiffness\n");
	for(int i=0;i<n*3;i++)
	{
		for(int j=0;j<n*3;j++)
		{
			printf("%.2lf ", globalStiffness[i][j]);
		}
		printf("\n");
	}
	*/

	printf("A\n");
	for(int i=0;i<n*3;i++)
	{
		for(int j=0;j<n*3;j++)
		{
			printf("%.2lf ", A[i][j]);
		}
		printf("\n");
	}

		printf("x\n");
	for(int i=0;i<n*3;i++)
		printf("%.2lf ",xt[i]);
			printf("\n");
				printf("V\n");
	for(int i=0;i<n*3;i++)
		printf("%.2lf ",v[i]);
			printf("\n");

	printf("b\n");
	for(int i=0;i<n*3;i++)
		printf("%.2lf ",b[i]);
			printf("\n");
}

Integrator::~Integrator(void)
{
}
