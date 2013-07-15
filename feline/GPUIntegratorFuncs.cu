#pragma once
#include <cuda.h>
#include <stdio.h>
#include "defines.h"
#include "GPUDataStructs.cuh"
#include "GPUPolarDecompose.cu"

#define BLOCK_SIZE 256
#define ALPHA 0.1
#define BETA 0.1

#define MAX_ITER 20
#define EPSIL 0.01

__constant__ float COEFFK, COEFFM, dt;

GPUElement* gpuptrElements;
GPUNode*   gpuptrNodes;
mulData*	gpuptrMulData;
float*   gpuptr_xt;//dynamic
float*   gpuptr_vt;//dynamic
float*	 gpuptr_extforces;//dynamic
float*	 gpuptr_b;//dynamic

//for CG
float* gpuptrR;
float* gpuptrD;
float* gpuptrQ;
CGVars* gpuptrVars;

__host__
void
gpuInitVars(int numele, int numnodes)
{
	cudaMalloc(&gpuptrElements, numele * sizeof(GPUElement));
	cudaMalloc(&gpuptrNodes, numnodes * sizeof(GPUNode));
	cudaMalloc(&gpuptrMulData, numele * sizeof(mulData));
	cudaMalloc(&gpuptr_xt, numnodes * 3 * sizeof(float));
	cudaMalloc(&gpuptr_vt, numnodes * 3 * sizeof(float));
	cudaMalloc(&gpuptr_extforces, numnodes * 3 * sizeof(float));
	cudaMalloc(&gpuptr_b, numnodes * 3 * sizeof(float));

	cudaMalloc(&gpuptrR, numnodes * 3 * sizeof(float));
	cudaMalloc(&gpuptrD, numnodes * 3 * sizeof(float));
	cudaMalloc(&gpuptrQ, numnodes * 3 * sizeof(float));
	cudaMalloc(&gpuptrVars, sizeof(CGVars));

	float coeffK = dt * BETA + dt * dt, coeffM = 1 + dt * ALPHA;
	float dt = 1.0/FPS;

	cudaMemcpyToSymbol("COEFFK", &coeffK, sizeof(float));
	cudaMemcpyToSymbol("COEFFM", &coeffM, sizeof(float));
	cudaMemcpyToSymbol("dt", &dt, sizeof(float));

	//testing
	float F[3][3] = { {1,2,3}, {3,2,1}, {1,3,2} };
	float R[3][3];
	gpuComputePolarDecomposition(F,R);

	for(int i=0;i<3;i++, printf("\n"))
		for(int j=0; j<3; j++)
			printf("%f ", R[i][j]);
	system("pause");
}

__host__
void
gpuUploadVars(GPUElement* gpuElements, GPUNode* gpuNodes,float* xt, 
			  float* vt, float* extforces, int numnodes, int numelements)
{
	cudaMemcpy(gpuptrElements, gpuElements, numelements * sizeof(GPUElement), cudaMemcpyHostToDevice);
	cudaMemcpy(gpuptrNodes, gpuNodes, numnodes * sizeof(GPUNode), cudaMemcpyHostToDevice);
	cudaMemcpy(gpuptr_xt, xt, numnodes * 3 * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(gpuptr_vt, vt, numnodes * 3 * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(gpuptr_extforces, extforces, numnodes * 3 * sizeof(float), cudaMemcpyHostToDevice);
}

__host__
void
gpuDownloadVars(float* xt, int numnodes)
{
	cudaMemcpy(xt, gpuptr_xt, numnodes * 3 * sizeof(float), cudaMemcpyDeviceToHost);
}

__host__
void
gpuUploadExtForces(float* extforces, int numnodes)
{
	cudaMemcpy(gpuptr_extforces, extforces, numnodes * 3 * sizeof(float), cudaMemcpyHostToDevice);
}


__host__
void
gpuDestroyVars()
{
	cudaFree(gpuptrElements);
	cudaFree(gpuptrNodes);
	cudaFree(gpuptrMulData);
	cudaFree(gpuptr_xt);
	cudaFree(gpuptr_vt);
	cudaFree(gpuptr_extforces);
	cudaFree(gpuptr_b);
	cudaFree(gpuptrR);
	cudaFree(gpuptrD);
	cudaFree(gpuptrQ);
	cudaFree(gpuptrVars);
}


__device__
void makeRK(float mat[12][12], float R[3][3])
{
	float RK[12][12];
	for(int i=0;i<4;i++)
	{
		for(int j=0;j<4;j++)
		{
			for(int a=0;a<3;a++)
			{
				for(int b=0;b<3;b++)
				{
					RK[a + i * 3][b + j * 3]=0;
						
					for(int c=0;c<3;c++)
						RK[a + i * 3][b + j * 3] += R[a][c] * mat[c + i * 3][b + j * 3];
				}
			}
		}
	}

	for(int i=0;i<12;i++)
		for(int j=0;j<12;j++)
			mat[i][j] = RK[i][j];
}


__device__
void makeRKRT(float mat[12][12], float R[3][3])
{
	float RKRT[12][12];

	for(int i=0;i<4;i++)
	{
		for(int j=i;j<4;j++)
		{
			for(int a=0;a<3;a++)
				for(int b=0;b<3;b++)
				{
					RKRT[a + i * 3][b + j * 3] = 0;
					for(int c=0;c<3;c++)
						RKRT[a + i * 3][b + j * 3] += mat[a + i * 3][c + j * 3] * R[b][c]; //R(c, b) but its RT so, R(b , c)
				}
		}
	}
	//lower triangle
	for(int i=1;i<4;i++)
		for(int j=0;j<i;j++)
		{
			for(int a=0;a<3;a++)
				for(int b=0;b<3;b++)
						RKRT[a + i * 3][b + j * 3] = RKRT[b + j * 3][a + i * 3];
		}

	for(int i=0;i<12;i++)
		for(int j=0;j<12;j++)
			mat[i][j] = RKRT[i][j];
}


__device__
void mulSystem(GPUElement* elements, mulData* solverData, float* x)
{
	int tid = threadIdx.x + blockIdx.x * BLOCK_SIZE;
	GPUElement* t_ele = &(elements[tid]);
	mulData* t_solvedata = &(solverData[tid]);

	float nodes[12];

	for(int i=0;i<4;i++)
	{
		nodes[i * 3] = x[t_ele->nodeindex[i] * 3];
		nodes[i * 3 + 1] = x[t_ele->nodeindex[i] * 3 + 1];
		nodes[i * 3 + 2] = x[t_ele->nodeindex[i] * 3 + 2];
	}

	for(int i=0;i<12;i++)
	{
		t_solvedata->product[i] = 0;
		for(int j=0;j<12;j++)
			t_solvedata->product[i] += t_solvedata->system[i][j] * nodes[j];
	}
}

__device__
void mulSystemGather(GPUNode* nodes, mulData* solverData, float* x)
{
	int tid = threadIdx.x + blockIdx.x * BLOCK_SIZE;
	GPUNode* node = &(nodes[tid]);

	int n = node->n;
	
	x[tid * 3] = 0;
	x[tid * 3 + 1] = 0;
	x[tid * 3 + 2] = 0;

	for(int i=0;i<n;i++)
	{
		int tetindex = node->elementindex[i][0];
		int nodeindex = node->elementindex[i][1];

		x[tid * 3] += solverData[tetindex].product[nodeindex * 3];
		x[tid * 3 + 1] += solverData[tetindex].product[nodeindex * 3 + 1];
		x[tid * 3 + 2] += solverData[tetindex].product[nodeindex * 3 + 2];
	}
}

__device__
void dot(float*a, float*b, float* out, int n) 
{
	__shared__ float temp[BLOCK_SIZE];
	int index = threadIdx.x;
	int element = index;

	float tmp = 0;

	while(element < n)
	{
		tmp += a[element] * b[element];
		element += BLOCK_SIZE;
	}

	temp[index] = tmp;

	__syncthreads();


	int i = BLOCK_SIZE >> 1;
	while(i>0)
	{
		if(index < i)
			temp[index] += temp[index + i];
		__syncthreads();
		i>>=1;
	}

	if(index == 0)
		*out = temp[0];
}

//step 1
//precompute
__global__
void precompute(GPUElement* elements, mulData* solverData, float* xt, float* vt, float* extforces, int numelements)
{
	int tid = threadIdx.x + blockIdx.x * BLOCK_SIZE;

	if(tid < numelements)
	{
		GPUElement* t_ele = &(elements[tid]);
		mulData* t_solvedata = &(solverData[tid]);

		float nodalmass = t_ele->nodalmass;

		float nodes[12], vel[12], F[3][3], R[3][3] = {{1,0,0},{0,1,0},{0,0,1}};

		for(int i=0;i<4;i++)
		{
			nodes[i * 3] = xt[t_ele->nodeindex[i] * 3];
			nodes[i * 3 + 1] = xt[t_ele->nodeindex[i] * 3 + 1];
			nodes[i * 3 + 2] = xt[t_ele->nodeindex[i] * 3 + 2];
		}

		for(int i=0;i<4;i++)
		{
			vel[i * 3] = vt[t_ele->nodeindex[i] * 3];
			vel[i * 3 + 1] = vt[t_ele->nodeindex[i] * 3 + 1];
			vel[i * 3 + 2] = vt[t_ele->nodeindex[i] * 3 + 2];
		}

		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
			{
				F[i][j] = nodes[j*3 + i] - nodes[9 + i];
			}

		//gpuComputePolarDecomposition(F,R);


		for(int i=0;i<12;i++)
			for(int j=0;j<12;j++)
				t_solvedata->system[i][j] = t_ele->unwarpK[i][j];
	
		for(int i=0;i<4;i++)
		{
			t_solvedata->b[i * 3] = extforces[t_ele->nodeindex[i] * 3];
			t_solvedata->b[i * 3 + 1] = extforces[t_ele->nodeindex[i] * 3 + 1];
			t_solvedata->b[i * 3 + 2] = extforces[t_ele->nodeindex[i] * 3 + 2];
		}

		makeRK(t_solvedata->system, R);

		for(int i=0;i<12;i++)
			for(int j=0;j<12;j++)
				t_solvedata->b[i] += t_solvedata->system[i][j] * t_ele->x0[j];
	
		makeRKRT(t_solvedata->system, R);

		for(int i=0;i<12;i++)
			for(int j=0;j<12;j++)
				t_solvedata->b[i] -= t_solvedata->system[i][j] * nodes[j];

		for(int i=0;i<12;i++)
			t_solvedata->b[i] = t_solvedata->b[i] * dt + nodalmass * vel[i];

		//final system matrix
		for(int i=0;i<12;i++)
			for(int j=0;j<12;j++)
			{
				t_solvedata->system[i][j] *= COEFFK;
				if(i==j)
					t_solvedata->system[i][i] += COEFFM * nodalmass;
			}
	}
}

//step 2
//precompute
__global__
void gatherB(GPUNode* nodes, mulData* solverData, float* b, int numnodes)
{
	int tid = threadIdx.x + blockIdx.x * BLOCK_SIZE;

	if(tid < numnodes)
	{
		GPUNode* node = &(nodes[tid]);

		int n = node->n;

		b[tid * 3] = 0;
		b[tid * 3 + 1] = 0;
		b[tid * 3 + 2] = 0;

		for(int i=0;i<n;i++)
		{
			int tetindex = node->elementindex[i][0];
			int nodeindex = node->elementindex[i][1];

			b[tid * 3] += solverData[tetindex].b[nodeindex * 3];
			b[tid * 3 + 1] += solverData[tetindex].b[nodeindex * 3 + 1];
			b[tid * 3 + 2] += solverData[tetindex].b[nodeindex * 3 + 2];
		}
	}
}

//step 1
//init CG
// x = velocity
__global__
void
initAx(GPUElement* elements, mulData* solverData, float* x, int numelements)
{
	int tid = threadIdx.x + blockIdx.x * BLOCK_SIZE;

	if(tid < numelements)
	{
		mulSystem(elements, solverData, x);
	}
}

//step2
//init CG
__global__
void
initRandD(GPUNode* nodes, mulData* solverData, float* r, float* d, float* b, int numnodes)
{
	int tid = threadIdx.x + blockIdx.x * BLOCK_SIZE;

	if(tid < numnodes)
	{
		mulSystemGather(nodes, solverData, r);
	
		//r = b-Ax
		r[tid * 3] = b[tid * 3] - r[tid * 3];
		r[tid * 3 + 1] = b[tid * 3 + 1] - r[tid * 3 + 1];
		r[tid * 3 + 2] = b[tid * 3 + 2] - r[tid * 3 + 2];

		//d=r
		d[tid * 3] = r[tid * 3];
		d[tid * 3 + 1] = r[tid * 3 + 1];
		d[tid * 3 + 2] = r[tid * 3 + 2];
	}
}

//step3
//init CG
//1 block, BLOCK_SIZE threads
__global__
void
initDeltaVars(CGVars* vars, float* r, int numnodes)
{
	__shared__ float rr;
	dot(r, r, &rr, numnodes * 3);
	
	if(threadIdx.x == 0)
	{
		vars->deltaNew = rr;
		vars->delta0 = vars->deltaNew;
	}
}

//step 4
//CG loop
//q = Ad
__global__
void
makeQprod(GPUElement* elements, mulData* solverData, float* d, int numelements)
{
	int tid = threadIdx.x + blockIdx.x * BLOCK_SIZE;

	if(tid < numelements)
	{
		mulSystem(elements, solverData, d);
	}
}

//step 5
//CG loop
//q = Ad
__global__
void
gatherQprod(GPUNode* nodes, mulData* solverData, float* q, int numnodes)
{
	int tid = threadIdx.x + blockIdx.x * BLOCK_SIZE;

	if(tid < numnodes)
	{
		mulSystemGather(nodes, solverData, q);
	}
}

//step 6
//CG Loop
//make vars
//1 block, BLOCK_SIZE threads
__global__
void
makeVars(CGVars* vars, float* d, float* q, float* r, int numnodes)
{
	float dq, rq, qq;
	dot(d,q,&dq,numnodes * 3);
	dot(r,q,&rq,numnodes * 3);
	dot(q,q,&qq,numnodes * 3);

	__syncthreads();

	if(threadIdx.x == 0)
	{
		vars->alpha = vars->deltaNew / dq;
		vars->deltaOld = vars->deltaNew;

		//r.r = r'.r' - 2*alpha*(r'.q) + alpha * alpha * (q.q)
		vars->deltaNew = vars->deltaNew - (2 * vars->alpha) * rq + (vars->alpha * vars->alpha) * qq;
		vars->beta = vars->deltaNew / vars->deltaOld;
	}
}

//step 7
//CG Loop
//make x, r, d
//x = velocity
__global__
void
makeXRandD(GPUNode* nodes, CGVars* vars, float *x, float* r, float* d, float* q, int numnodes)
{
	int tid = threadIdx.x + blockIdx.x * BLOCK_SIZE;
	if(tid < numnodes)
	{
		x[tid * 3] = x[tid * 3] + vars->alpha * d[tid * 3];
		x[tid * 3 + 1] = x[tid * 3 + 1] + vars->alpha * d[tid * 3 + 1];
		x[tid * 3 + 2] = x[tid * 3 + 2] + vars->alpha * d[tid * 3 + 2];

		r[tid * 3] = r[tid * 3] - vars->alpha * q[tid * 3];
		r[tid * 3 + 1] = r[tid * 3 + 1] - vars->alpha * q[tid * 3 + 1];
		r[tid * 3 + 2] = r[tid * 3 + 2] - vars->alpha * q[tid * 3 + 2];

		d[tid * 3] = r[tid * 3] + vars->beta * d[tid * 3];
		d[tid * 3 + 1] = r[tid * 3 + 1] + vars->beta * d[tid * 3 + 1];
		d[tid * 3 + 2] = r[tid * 3 + 2] + vars->beta * d[tid * 3 + 2];
	}
} 

//step 8
//make x(t+1)
__global__
void
integrate(float *x, float* v, int numnodes)
{
	int tid = threadIdx.x + blockIdx.x * BLOCK_SIZE;
	if(tid < numnodes)
	{
		x[tid * 3] = x[tid * 3] + dt * v[tid * 3];
		x[tid * 3 + 1] = x[tid * 3 + 1] + dt * v[tid * 3 + 1];
		x[tid * 3 + 2] = x[tid * 3 + 2] + dt * v[tid * 3 + 2];
	}
}

__host__
void
gpuTimeStep(int numelements, int numnodes)
{
	const int num_blocks_ele = (numelements/BLOCK_SIZE) + 1;
	const int num_blocks_node = (numnodes/BLOCK_SIZE) + 1;

	printf("Started\n");
	
	//printf("Precompute\n");
	precompute<<<num_blocks_ele, BLOCK_SIZE>>>(gpuptrElements, gpuptrMulData, gpuptr_xt, gpuptr_vt, gpuptr_extforces, numelements);
	
	//printf("GatherB\n");
	gatherB<<<num_blocks_node, BLOCK_SIZE>>>(gpuptrNodes, gpuptrMulData, gpuptr_b, numnodes);

	//printf("InitAx\n");
	initAx<<<num_blocks_ele, BLOCK_SIZE>>>(gpuptrElements, gpuptrMulData, gpuptr_vt, numelements);

	//printf("InitRandD\n");
	initRandD<<<num_blocks_node, BLOCK_SIZE>>>(gpuptrNodes, gpuptrMulData, gpuptrR, gpuptrD, gpuptr_b, numnodes);

	//printf("InitDeltaVars\n");
	initDeltaVars<<<1, BLOCK_SIZE>>>(gpuptrVars, gpuptrR, numnodes);

	int i=0;

	CGVars vars;
	cudaMemcpy(&vars, gpuptrVars, sizeof(CGVars), cudaMemcpyDeviceToHost);

	printf("%f %f %f %f %f\n", vars.delta0, vars.deltaNew, vars.deltaOld, vars.alpha, vars.beta);
	//system("pause");

	printf("Loop Started");

	while(i < MAX_ITER && vars.deltaNew > (EPSIL * EPSIL) * vars.delta0)
	{
		//printf("MakeQProd\n");
		makeQprod<<<num_blocks_ele, BLOCK_SIZE>>>(gpuptrElements, gpuptrMulData, gpuptrD, numelements);

		//printf("gatherQprod\n");
		gatherQprod<<<num_blocks_node, BLOCK_SIZE>>>(gpuptrNodes, gpuptrMulData, gpuptrQ, numnodes);

		//printf("MakeVars\n");
		makeVars<<<1, BLOCK_SIZE>>>(gpuptrVars, gpuptrD, gpuptrQ, gpuptrR, numnodes);

		//printf("MakeXRandD\n");
		makeXRandD<<<num_blocks_node, BLOCK_SIZE>>>(gpuptrNodes, gpuptrVars, gpuptr_vt, gpuptrR, gpuptrD, gpuptrQ, numnodes);

		//printf("Memcpy\n");
		cudaMemcpy(&vars, gpuptrVars, sizeof(CGVars), cudaMemcpyDeviceToHost);
		//printf("%f %f %f %f %f\n", vars.delta0, vars.deltaNew, vars.deltaOld, vars.alpha, vars.beta);
		i++;
		//printf("End iter\n");
	}

	printf("Loop Ended: %d ", i);

	//printf("Integrate\n");
	integrate<<<num_blocks_node, BLOCK_SIZE>>>(gpuptr_xt, gpuptr_vt, numnodes);
}

