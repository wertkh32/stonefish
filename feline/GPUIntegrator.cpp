#include "GPUIntegrator.h"


GPUIntegrator::GPUIntegrator(Mesh* _mesh, ConstrainedRows* r)
{
	mesh = _mesh;
	numnodes = mesh->getNoOfNodes();
	numelements = mesh->getNoOfElements();

	gpuElements = (GPUElement*)malloc(sizeof(GPUElement) * numelements * 3);
	gpuNodes = (GPUNode*)malloc(sizeof(GPUNode) * numnodes);

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
		GenMatrix<float,12,12>& stiff = *(mesh->elements[i]->getStiffnessMat());
		for(int a=0;a<12;a++)
			for(int b=0;b<12;b++)
				gpuElements[i].unwarpK[a][b] = stiff(a,b);

		for(int a=0;a<4;a++)
		{
			gpuElements[i].x0[a * 3] = mesh->elements[i]->nodes[a]->pos.x;
			gpuElements[i].x0[a * 3 + 1] = mesh->elements[i]->nodes[a]->pos.y;
			gpuElements[i].x0[a * 3 + 2] = mesh->elements[i]->nodes[a]->pos.z;

			gpuElements[i].nodeindex[a] = mesh->nodeIndices[i][a];
		}
		
		gpuElements[i].nodalmass = (mesh->elements[i]->getDensity() * mesh->elements[i]->getVolume())/4.0;

	}
}

void
GPUIntegrator::assembleGPUNodes()
{
	for(int i=0;i<numnodes;i++)
	{
		gpuNodes[i].n = 0;
		for(int a=0;a<numelements;a++)
		{
			for(int b=0;b<4;b++)
			{
				if(mesh->nodeIndices[a][b] == i)
				{
					gpuNodes[i].elementindex[gpuNodes[i].n][0] = a;
					gpuNodes[i].elementindex[gpuNodes[i].n][1] = b;
					gpuNodes[i].n++;
					break;
				}
			}
		}
	}
}
	
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
	//assembleExtForce();
	//gpuUploadExtForces(extforces, numnodes);
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
