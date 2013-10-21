#include "GPUIntegrator.h"


GPUIntegrator::GPUIntegrator(MESH* _mesh, ConstrainedRows* r)
{
	mesh = _mesh;
	numnodes = mesh->getNoOfNodes();
	numelements = mesh->getNoOfElements();

	int numblocksperele = (numelements / BLOCK_SIZE) + 1;
	int numblockpernode = (numnodes / NODE_BLOCK_SIZE) + 1;

	gpuElements = (GPUElement*)malloc(sizeof(GPUElement) * numblocksperele);
	gpuNodes = (GPUNode*)malloc(sizeof(GPUNode) * numblockpernode);

	xt = (float*)malloc(sizeof(float) * numnodes * 3);
	vt = (float*)malloc(sizeof(float) * numnodes * 3);
	extforces = (float*)malloc(sizeof(float) * numnodes * 3);
	mass = (float*)malloc(sizeof(float) * numnodes * 3);
	allowed = (char*)malloc(sizeof(bool) * numnodes);

	for(int i=0;i<numnodes;i++)
		allowed[i] = 0;

	for(int i=0;i<r->list.size();i++)
	{
		printf("BOOOOOO %d | %d\n",r->list[i] / 3, r->list[i] % 3);
		allowed[r->list[i] / 3] |= (1<<(r->list[i] % 3));
	}

	initVars();
	copyVarstoGPU();
}

void 
GPUIntegrator::assembleGPUElements()
{
	for(int i=0;i<numelements;i++)
	{
		int tid = i % BLOCK_SIZE;
		int bid = i / BLOCK_SIZE;

		float x0[(NUM_NODES_PER_ELE * 3)];
		float f0[(NUM_NODES_PER_ELE * 3)] = {0};

		GenMatrix<float,(NUM_NODES_PER_ELE * 3),(NUM_NODES_PER_ELE * 3)>& stiff = (mesh->elements[i]->getStiffnessMat());

		for(int a=0;a<3;a++)
			for(int b=0;b<3;b++)
				gpuElements[bid].B[a][b][tid] = mesh->elements[i]->B[a][b];

		gpuElements[bid].c1[tid] = mesh->elements[i]->con1;
		gpuElements[bid].c2[tid] = mesh->elements[i]->con2;

		for(int a=0;a<NUM_NODES_PER_ELE;a++)
		{
			x0[a*3] = mesh->elements[i]->nodes[a]->pos.x;
			x0[a*3 + 1] = mesh->elements[i]->nodes[a]->pos.y;
			x0[a*3 + 2] = mesh->elements[i]->nodes[a]->pos.z;
		}

		for(int a=0;a<(NUM_NODES_PER_ELE * 3);a++)
			for(int b=0;b<(NUM_NODES_PER_ELE * 3);b++)
			{
				f0[a] += stiff(a,b) * x0[b];
			}

		#if defined(_QUAD_TET_) && defined(_BERSTEIN_POLY_)

		for(int a=0;a<(NUM_NODES_PER_ELE * 3);a++)
			for(int b=0;b<(NUM_NODES_PER_ELE * 3);b++)
				gpuElements[bid].system[a][b][tid] = stiff(a,b);
		
		#endif

		for(int a=0;a<(NUM_NODES_PER_ELE * 3);a++)
		{
			gpuElements[bid].f0[a][tid] = f0[a];
		}

		for(int a=0;a<NUM_NODES_PER_ELE;a++)
		{
			gpuElements[bid].nodeindex[a][tid] = mesh->nodeIndices[i][a];
		}

		//Matrix3d& inv = mesh->elements[i]->getUndeformShapeMatInv();

		/*
		for(int a=0;a<3;a++)
			for(int b=0;b<3;b++)
			{
				gpuElements[bid].undefShapeMatInv[a][b][tid] = inv(a,b);
			}
		*/

		//gpuElements[bid].nodalmass[tid] = (mesh->elements[i]->getDensity() * mesh->elements[i]->getVolume())/4.0;

	}
}

void
GPUIntegrator::assembleGPUNodes()
{
	for(int i=0;i<numnodes;i++)
	{
		int tid = i % NODE_BLOCK_SIZE;
		int bid = i / NODE_BLOCK_SIZE;
		
		int n = 0;
		
		for(int b=0;b<NODE_THREADS;b++)
			gpuNodes[bid].n[b][tid] = 0;

		for(int a=0;a<numelements;a++)
		{
			for(int b=0;b<NUM_NODES_PER_ELE;b++)
			{
				if(mesh->nodeIndices[a][b] == i)
				{
					int t = n%NODE_THREADS;
					gpuNodes[bid].elementindex[gpuNodes[bid].n[t][tid]][0][t][tid] = a;
					gpuNodes[bid].elementindex[gpuNodes[bid].n[t][tid]][1][t][tid] = b;
					n++;
					(gpuNodes[bid].n[t][tid])++;
					break;
				}
			}
		}
	}
}

void
GPUIntegrator::assembleLumpedMass()
{
	for(int i=0;i<numnodes*3;i++)
		mass[i] = 0;

	for(int i=0;i<mesh->elements.size();i++)
	{
		//float elenodemass = (mesh->elements[i]->getDensity() * mesh->elements[i]->getVolume()) /10.0;
		for(int j=0;j<NUM_NODES_PER_ELE;j++)
		{
			mass[mesh->nodeIndices[i][j]] += mesh->elements[i]->nodemass[j];
			mass[mesh->nodeIndices[i][j] + numnodes] += mesh->elements[i]->nodemass[j];
			mass[mesh->nodeIndices[i][j] + numnodes * 2] += mesh->elements[i]->nodemass[j];
		}
	}
}

void GPUIntegrator::assembleXt()
{
	for(int i=0;i<numnodes;i++)
	{
		xt[i] = mesh->nodes[i]->pos_t.x;
		xt[i + numnodes] = mesh->nodes[i]->pos_t.y;
		xt[i + numnodes * 2] = mesh->nodes[i]->pos_t.z;	
	}
}

void GPUIntegrator::assembleVt()
{
	for(int i=0;i<numnodes;i++)
	{	
		vt[i] = mesh->nodes[i]->vec_t.x;
		vt[i + numnodes] = mesh->nodes[i]->vec_t.y;
		vt[i + numnodes * 2] = mesh->nodes[i]->vec_t.z;
	}
}

void GPUIntegrator::assembleExtForce()
{
	for(int i=0;i<numnodes;i++)
	{
		extforces[i] = mesh->nodes[i]->force.x;
		extforces[i + numnodes] = mesh->nodes[i]->force.y;
		extforces[i + numnodes * 2] = mesh->nodes[i]->force.z;
	}
}

void
GPUIntegrator::initVars()
{
	gpuInitVars(numelements, numnodes);
	assembleGPUElements();
	assembleGPUNodes();
	assembleXt();
	assembleVt();
	assembleExtForce();
	assembleLumpedMass();
}

void
GPUIntegrator::copyVarstoGPU()
{
	gpuUploadVars(gpuElements, gpuNodes, xt, vt, extforces, mass, allowed, numnodes, numelements);
}

void
GPUIntegrator::timeStep()
{
	assembleExtForce();
	gpuUploadExtForces(extforces, numnodes);
	gpuTimeStep(numelements, numnodes);
	updatePositions();
}

void
GPUIntegrator::updatePositions()
{
	gpuDownloadVars(xt, numnodes);

	for(int i=0;i<numnodes;i++)
	{
		mesh->nodes[i]->pos_t.x = xt[i];
		mesh->nodes[i]->pos_t.y = xt[i + numnodes];
		mesh->nodes[i]->pos_t.z = xt[i + numnodes * 2];
	}
}


GPUIntegrator::~GPUIntegrator(void)
{
	gpuDestroyVars();
}
