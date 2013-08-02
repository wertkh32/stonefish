#include "GPUIntegrator.h"


GPUIntegrator::GPUIntegrator(Mesh* _mesh, ConstrainedRows* r)
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

		GenMatrix<float,12,12>& stiff = *(mesh->elements[i]->getStiffnessMat());
		for(int a=0;a<12;a++)
			for(int b=0;b<12;b++)
				gpuElements[bid].unwarpK[a][b][tid] = stiff(a,b);

		for(int a=0;a<4;a++)
		{
			gpuElements[bid].x0[a * 3][tid] = mesh->elements[i]->nodes[a]->pos.x;
			gpuElements[bid].x0[a * 3 + 1][tid] = mesh->elements[i]->nodes[a]->pos.y;
			gpuElements[bid].x0[a * 3 + 2][tid] = mesh->elements[i]->nodes[a]->pos.z;

			gpuElements[bid].nodeindex[a][tid] = mesh->nodeIndices[i][a];
		}

		Matrix3d& inv = mesh->elements[i]->getUndeformShapeMatInv();

		for(int a=0;a<3;a++)
			for(int b=0;b<3;b++)
			{
				gpuElements[bid].undefShapeMatInv[a][b][tid] = inv(a,b);
			}
		
		gpuElements[bid].nodalmass[tid] = (mesh->elements[i]->getDensity() * mesh->elements[i]->getVolume())/4.0;

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
		
		for(int b=0;b<4;b++)
			gpuNodes[bid].n[tid][b] = 0;

		for(int a=0;a<numelements;a++)
		{
			for(int b=0;b<4;b++)
			{
				if(mesh->nodeIndices[a][b] == i)
				{
					int t = n%4;
					gpuNodes[bid].elementindex[gpuNodes[bid].n[tid][t]][0][tid][t] = a;
					gpuNodes[bid].elementindex[gpuNodes[bid].n[tid][t]][1][tid][t] = b;
					n++;
					(gpuNodes[bid].n[tid][t])++;
					break;
				}
			}
		}
	}
}
/*
void
GPUIntegrator::assembleLumpedMass()
{
	for(int i=0;i<n*3;i++)
		mass[i] = 0;

	for(int i=0;i<mesh->elements.size();i++)
	{
		float elenodemass = (mesh->elements[i]->getDensity() * mesh->elements[i]->getVolume()) /4;
		for(int j=0;j<4;j++)
		{
			mass[mesh->nodeIndices[i][j] * 3] += elenodemass;
			mass[mesh->nodeIndices[i][j] * 3 + 1] += elenodemass;
			mass[mesh->nodeIndices[i][j] * 3 + 2] += elenodemass;
		}
	}
}
*/
void GPUIntegrator::assembleXt()
{
	for(int i=0;i<numnodes;i++)
	{
		xt[i * 3] = mesh->nodes[i]->pos_t.x;
		xt[i * 3 + 1] = mesh->nodes[i]->pos_t.y;
		xt[i * 3 + 2] = mesh->nodes[i]->pos_t.z;	
	}
}

void GPUIntegrator::assembleVt()
{
	for(int i=0;i<numnodes;i++)
	{	
		vt[i * 3] = mesh->nodes[i]->vec_t.x;
		vt[i * 3 + 1] = mesh->nodes[i]->vec_t.y;
		vt[i * 3 + 2] = mesh->nodes[i]->vec_t.z;
	}
}

void GPUIntegrator::assembleExtForce()
{
	for(int i=0;i<numnodes;i++)
	{
		extforces[i * 3] = mesh->nodes[i]->force.x;
		extforces[i * 3 + 1] = mesh->nodes[i]->force.y;
		extforces[i * 3 + 2] = mesh->nodes[i]->force.z;
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
}

void
GPUIntegrator::copyVarstoGPU()
{
	gpuUploadVars(gpuElements, gpuNodes, xt, vt, extforces, numnodes, numelements);
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
		mesh->nodes[i]->pos_t.x = xt[i * 3];
		mesh->nodes[i]->pos_t.y = xt[i * 3 + 1];
		mesh->nodes[i]->pos_t.z = xt[i * 3 + 2];
	}
}


GPUIntegrator::~GPUIntegrator(void)
{
	gpuDestroyVars();
}
