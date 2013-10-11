#pragma once
#include <cuda.h>
#include <stdio.h>
#include "defines.h"

#ifdef _QUAD_TET_

#include "GPUDataStructs.cuh"
#include "GPUPolarDecompose.cu"

//#define BLOCK_SIZE 512

#define ALPHA 0.4
#define BETA 0.5

#define MAX_ITER 20
#define EPSIL 0.07

__constant__ float COEFFK, COEFFM, dt;

GPUElement* gpuptrElements;
GPUNode*   gpuptrNodes;
mulData*	gpuptrMulData;
float*   gpuptr_xt;//dynamic
float*   gpuptr_vt;//dynamic
float*	 gpuptr_extforces;//dynamic
float*	 gpuptr_mass;//static
float*	 gpuptr_b;//dynamic
char*	 gpuptr_allowed;

//for CG
float* gpuptrR;
float* gpuptrD;
float* gpuptrQ;
CGVars* gpuptrVars;

static void HandleError( cudaError_t err,
                         const char *file,
                         int line ) {
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ),
                file, line );
		system("pause");
        exit( EXIT_FAILURE );
    }
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))


__host__
void
gpuInitVars(int numele, int numnodes)
{
	int numblocksperele = (numele / BLOCK_SIZE) + 1;
	int numblockpernode = (numnodes / NODE_BLOCK_SIZE) + 1;

	//cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);

	HANDLE_ERROR( cudaMalloc(&gpuptrElements, numblocksperele * sizeof(GPUElement)) );
	HANDLE_ERROR( cudaMalloc(&gpuptrMulData, numblocksperele * sizeof(mulData)) );
	HANDLE_ERROR( cudaMalloc(&gpuptrNodes, numblockpernode * sizeof(GPUNode)) );
	HANDLE_ERROR( cudaMalloc(&gpuptr_xt, numnodes * 3 * sizeof(float)) );
	HANDLE_ERROR( cudaMalloc(&gpuptr_vt, numnodes * 3 * sizeof(float)) );
	HANDLE_ERROR( cudaMalloc(&gpuptr_extforces, numnodes * 3 * sizeof(float)) );
	HANDLE_ERROR( cudaMalloc(&gpuptr_mass, numnodes * 3 * sizeof(float)) );
	HANDLE_ERROR( cudaMalloc(&gpuptr_b, numnodes * 3 * sizeof(float)) );
	HANDLE_ERROR( cudaMalloc(&gpuptr_allowed, numnodes * sizeof(char)) );


	HANDLE_ERROR( cudaMalloc(&gpuptrR, numnodes * 3 * sizeof(float)) );
	HANDLE_ERROR( cudaMalloc(&gpuptrD, numnodes * 3 * sizeof(float)) );
	HANDLE_ERROR( cudaMalloc(&gpuptrQ, numnodes * 3 * sizeof(float)) );
	HANDLE_ERROR( cudaMalloc(&gpuptrVars, sizeof(CGVars)) );

	float ddt = 1.0/FPS;
	float coeffK = ddt * BETA + ddt * ddt, coeffM = 1 + ddt * ALPHA;
	

	HANDLE_ERROR( cudaMemcpyToSymbol("COEFFK", &coeffK, sizeof(float)) );
	HANDLE_ERROR( cudaMemcpyToSymbol("COEFFM", &coeffM, sizeof(float)) );
	HANDLE_ERROR( cudaMemcpyToSymbol("dt", &ddt, sizeof(float)) );

		cudaDeviceSynchronize();
		cudaError_t error = cudaGetLastError();
		if(error != cudaSuccess)
		{
			// print the CUDA error message and exit
			printf("CUDA error: %s\n", cudaGetErrorString(error));
			//exit(-1);
			system("pause");
		}
}

__host__
void
gpuUploadVars(GPUElement* gpuElements, GPUNode* gpuNodes,float* xt, 
			  float* vt, float* extforces, float* mass, char* allowed, int numnodes, int numelements)
{
	int numblocksperele = (numelements / BLOCK_SIZE) + 1;
	int numblockpernode = (numnodes / NODE_BLOCK_SIZE) + 1;

	HANDLE_ERROR( cudaMemcpy(gpuptrElements, gpuElements, numblocksperele * sizeof(GPUElement), cudaMemcpyHostToDevice) );
	HANDLE_ERROR( cudaMemcpy(gpuptrNodes, gpuNodes, numblockpernode * sizeof(GPUNode), cudaMemcpyHostToDevice) );
	HANDLE_ERROR( cudaMemcpy(gpuptr_xt, xt, numnodes * 3 * sizeof(float), cudaMemcpyHostToDevice) );
	HANDLE_ERROR( cudaMemcpy(gpuptr_vt, vt, numnodes * 3 * sizeof(float), cudaMemcpyHostToDevice) );
	HANDLE_ERROR( cudaMemcpy(gpuptr_extforces, extforces, numnodes * 3 * sizeof(float), cudaMemcpyHostToDevice) );
	HANDLE_ERROR( cudaMemcpy(gpuptr_mass, mass, numnodes * 3 * sizeof(float), cudaMemcpyHostToDevice) );
	HANDLE_ERROR( cudaMemcpy(gpuptr_allowed, allowed, numnodes * sizeof(char), cudaMemcpyHostToDevice) );

		cudaDeviceSynchronize();
		cudaError_t error = cudaGetLastError();
		if(error != cudaSuccess)
		{
			// print the CUDA error message and exit
			printf("CUDA error: %s\n", cudaGetErrorString(error));
			//exit(-1);
			system("pause");
		}
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

	cudaDeviceSynchronize();
	cudaError_t error = cudaGetLastError();
	if(error != cudaSuccess)
	{
		// print the CUDA error message and exit
		printf("CUDA error: %s\n", cudaGetErrorString(error));
		//exit(-1);
		system("pause");
	}
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
	cudaFree(gpuptr_mass);
	cudaFree(gpuptr_b);
	cudaFree(gpuptrR);
	cudaFree(gpuptrD);
	cudaFree(gpuptrQ);
	cudaFree(gpuptrVars);
}

__device__
void makeFU(float f0[30][BLOCK_SIZE], float R[3][3], float out[30])
{
	int ltid = threadIdx.x;
	float x[30];

	#pragma unroll 30
	for(int i=0;i<30;i++)
		x[i] = f0[i][ltid];

	#pragma unroll 10
	for(int i=0;i<10;i++)
		#pragma unroll 3
		for(int j=0;j<3;j++)
		{
			out[i*3 + j] = 0;
			#pragma unroll 3
			for(int k=0;k<3;k++)
			out[i*3+j] += R[j][k] * x[i*3 + k];
		}		
}


__device__
void mulK(float temp[3], float x[30], float b[4][3], float con1, float con2)
{
	int etid = threadIdx.x / BLOCK_SIZE;
	float con3 = (con1 - con2)/2.0;
	/*
	switch(etid)
	{

	}*/

}

//xt will be malnipulated and used as a temp array. do not further malnipulate xt
__device__
void makeRKRT(float system[30][30][BLOCK_SIZE], float R[3][3], float xt[30], float b[30])
{
	float temp[30];

	#pragma unroll 10
	for(int i=0;i<10;i++)
		#pragma unroll 3
		for(int j=0;j<3;j++)
		{
			temp[i*3 + j] = 0;
			#pragma unroll 3
			for(int k=0;k<3;k++)
			temp[i*3+j] += R[k][j] * xt[i*3 + k]; //RT first
		}

	//mulK(temp,xt, B,c1,c2);
	#pragma unroll 30
	for(int i=0;i<30;i++)
	{
		xt[i] = 0;
		#pragma unroll 30
		for(int j=0;j<30;j++)
		xt[i] += system[i][j][threadIdx.x] * temp[j];
	}		

	#pragma unroll 10
	for(int i=0;i<10;i++)
		#pragma unroll 3
		for(int j=0;j<3;j++)
		{
			#pragma unroll 3
			for(int k=0;k<3;k++)
			b[i*3+j] -= R[j][k] * xt[i*3 + k];
		}

}

__device__
void mulSystem(GPUElement* elements, mulData* solverData, float* x, int numelements)
{
	int bid = blockIdx.x;
	int ltid = threadIdx.x % BLOCK_SIZE;
	int etid = threadIdx.x / BLOCK_SIZE;
	int tid = ltid + blockIdx.x * BLOCK_SIZE;
	
	GPUElement* t_ele = &(elements[bid]);
	mulData* t_solvedata = &(solverData[bid]);

	__shared__ float nodes[30][BLOCK_SIZE]; 
	__shared__ float R[3][3][BLOCK_SIZE];
	__shared__ float B[4][3][BLOCK_SIZE];

	float temp[3];

	if(tid < numelements)
	{
		if(etid == 0)
		{
			#pragma unroll 3
			for(int i=0;i<3;i++)
				#pragma unroll 3
				for(int j=0;j<3;j++)
					R[i][j][ltid] = t_solvedata->R[i][j][ltid];

			#pragma unroll 3
			for(int i=0;i<3;i++)
				#pragma unroll 3
				for(int j=0;j<3;j++)
					B[i][j][ltid] = t_ele->B[i][j][ltid];

			B[3][0][ltid] = -B[0][0][ltid] -B[1][0][ltid] -B[2][0][ltid];
			B[3][1][ltid] = -B[0][1][ltid] -B[1][1][ltid] -B[2][1][ltid];
			B[3][2][ltid] = -B[0][2][ltid] -B[1][2][ltid] -B[2][2][ltid];
		}
	}

	__syncthreads();

	if(tid < numelements)
	{
		//rotate by x by RT first
		int index = t_ele->nodeindex[etid][ltid];

		float temp2[3];
		temp2[0] = x[index * 3];
		temp2[1] = x[index * 3 + 1];
		temp2[2] = x[index * 3 + 2];

		#pragma unroll 3
		for(int j=0;j<3;j++)
		{
			temp[j] = 0;
			#pragma unroll 3
			for(int k=0;k<3;k++)
			temp[j] += R[k][j][ltid] * temp2[k];
		}

		nodes[etid * 3][ltid] = temp[0];
		nodes[etid * 3 + 1][ltid] = temp[1];
		nodes[etid * 3 + 2][ltid] = temp[2];
	}

	__syncthreads();

	if(tid < numelements)
	{
		temp[0] = 0;
		temp[1] = 0;	
		temp[2] = 0;

		#pragma unroll 30
		for(int j=0;j<30;j++)
		{
			#pragma unroll 3
			for(int i=0;i<3;i++)
				temp[i] += t_ele->system[etid * 3 + i][j][ltid] * nodes[j][ltid];
		}		


		#pragma unroll 3
		for(int j=0;j<3;j++)
		{
			float temp3 = 0;
			#pragma unroll 3
			for(int k=0;k<3;k++)
				temp3 += R[j][k][ltid] * temp[k];
			temp3 *= COEFFK;
			//temp3 += nodalmass * nodes[i*3 + j];
			t_solvedata->product[etid*3 + j][ltid] = temp3;
		}
	}
}

__device__
void dot(float*a, float*b, float* out, int n) 
{
	__shared__ float temp[DOT_BLOCK_SIZE];
	int index = threadIdx.x;
	int element = index;

	float tmp = 0;

	while(element < n)
	{
		tmp += a[element] * b[element];
		element += DOT_BLOCK_SIZE;
	}

	temp[index] = tmp;

	__syncthreads();


	int i = DOT_BLOCK_SIZE >> 1;
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
void precompute(GPUElement* elements, mulData* solverData, float* xt, int numelements)
{
	int tid = threadIdx.x + blockIdx.x * BLOCK_SIZE;
	int bid = blockIdx.x;
	int ltid = threadIdx.x;

	if(tid < numelements)
	{
		GPUElement* t_ele = &(elements[bid]);
		mulData* t_solvedata = &(solverData[bid]);

		float nodes[30], b[30], R[3][3]={0}, D[3][3];
		int index[10];
			

		#pragma unroll 10
		for(int i=0;i<10;i++)
		{
			index[i] = t_ele->nodeindex[i][ltid];
			nodes[i * 3] = xt[index[i] * 3];
			nodes[i * 3 + 1] = xt[index[i] * 3 + 1];
			nodes[i * 3 + 2] = xt[index[i] * 3 + 2];
		}

		#pragma unroll 3
		for(int i=0;i<3;i++)
			#pragma unroll 3
			for(int j=0;j<3;j++)
				D[i][j] = t_ele->B[i][j][ltid];

		#pragma unroll 3
		for(int i=0;i<3;i++)
			#pragma unroll 3
			for(int j=0;j<3;j++)
					R[i][j] = (nodes[i] - nodes[9 + i]) * D[0][j] + 
							  (nodes[3 + i] - nodes[9 + i]) * D[1][j] + 
							  (nodes[6 + i] - nodes[9 + i]) * D[2][j];

		gpuComputePolarDecomposition(R);
	
		#pragma unroll 3
		for(int i=0;i<3;i++)
			#pragma unroll 3
			for(int j=0;j<3;j++)
				t_solvedata->R[i][j][ltid] =  R[i][j];

		makeFU(t_ele->f0,R,b);
	
		makeRKRT(t_ele->system, R, nodes, b);

		#pragma unroll 30
		for(int i=0;i<30;i++)
			t_solvedata->b[i][ltid] = b[i] * dt;// + nodalmass * vt[index[(i/3)] * 3 + (i%3)];

	}
}

//step 2
//precompute
__global__
void gatherB(GPUNode* nodes, mulData* solverData, float* b, float* mass, float* vt, float* extforces, char* allowed,int numnodes)
{
	int groupid = threadIdx.x % NODE_BLOCK_SIZE;// / NODE_THREADS;
	int grouptid = threadIdx.x / NODE_BLOCK_SIZE; //% NODE_THREADS;
	int nodeno = blockIdx.x * NODE_BLOCK_SIZE + groupid;

	__shared__ float cache[NODE_THREADS][NODE_BLOCK_SIZE][3];
	GPUNode* node = &(nodes[blockIdx.x]);
	int n = node->n[grouptid][groupid];
	
	if(nodeno < numnodes)
	{

		cache[grouptid][groupid][0] = 0;
		cache[grouptid][groupid][1] = 0;
		cache[grouptid][groupid][2] = 0;


		for(int i=0;i<n;i++)
		{
			int tetindex = node->elementindex[i][0][grouptid][groupid] / BLOCK_SIZE;
			int tetindex2 = node->elementindex[i][0][grouptid][groupid] % BLOCK_SIZE;
			int nodeindex = node->elementindex[i][1][grouptid][groupid];

			cache[grouptid][groupid][0] += solverData[tetindex].b[nodeindex * 3][tetindex2];
			cache[grouptid][groupid][1] += solverData[tetindex].b[nodeindex * 3 + 1][tetindex2];
			cache[grouptid][groupid][2] += solverData[tetindex].b[nodeindex * 3 + 2][tetindex2];
		}
	}

	__syncthreads();

	if(nodeno < numnodes)
	{
		if(grouptid == 0)
		{
			b[nodeno * 3]     = cache[0][groupid][0] + cache[1][groupid][0] /*+ cache[2][groupid][0] + cache[3][groupid][0]*/ + mass[nodeno * 3] * vt[nodeno * 3] + extforces[nodeno * 3] * dt;
			b[nodeno * 3 + 1] = cache[0][groupid][1] + cache[1][groupid][1] /*+ cache[2][groupid][1] + cache[3][groupid][1]*/ + mass[nodeno * 3 + 1] * vt[nodeno * 3 + 1] + extforces[nodeno * 3 + 1] * dt;
			b[nodeno * 3 + 2] = cache[0][groupid][2] + cache[1][groupid][2] /*+ cache[2][groupid][2] + cache[3][groupid][2]*/ + mass[nodeno * 3 + 2] * vt[nodeno * 3 + 2] + extforces[nodeno * 3 + 2] * dt;

			char bitsy = allowed[nodeno];
			if(bitsy & 1)
				vt[nodeno * 3] = 0;
			if(bitsy & 2)
				vt[nodeno * 3 + 1] = 0;
			if(bitsy & 4)
				vt[nodeno * 3 + 2] = 0;

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
		mulSystem(elements, solverData, x, numelements);
}

//step2
//init CG
__global__
void
initRandD(GPUNode* nodes, mulData* solverData, float* r, float* d, float* b, float* mass, float* vt, char* allowed, int numnodes)
{
	int groupid = threadIdx.x % NODE_BLOCK_SIZE;// / NODE_THREADS;
	int grouptid = threadIdx.x / NODE_BLOCK_SIZE; //% NODE_THREADS;
	int nodeno = blockIdx.x * NODE_BLOCK_SIZE + groupid;

	__shared__ float cache[NODE_THREADS][NODE_BLOCK_SIZE][3];
	GPUNode* node = &(nodes[blockIdx.x]);
	int n = node->n[grouptid][groupid];
	
	if(nodeno < numnodes)
	{

		cache[grouptid][groupid][0] = 0;
		cache[grouptid][groupid][1] = 0;
		cache[grouptid][groupid][2] = 0;

		for(int i=0;i<n;i++)
		{
			int tetindex = node->elementindex[i][0][grouptid][groupid] / BLOCK_SIZE;
			int tetindex2 = node->elementindex[i][0][grouptid][groupid] % BLOCK_SIZE;
			int nodeindex = node->elementindex[i][1][grouptid][groupid];

			cache[grouptid][groupid][0] += solverData[tetindex].product[nodeindex * 3][tetindex2];
			cache[grouptid][groupid][1] += solverData[tetindex].product[nodeindex * 3 + 1][tetindex2];
			cache[grouptid][groupid][2] += solverData[tetindex].product[nodeindex * 3 + 2][tetindex2];
		}
	}

	__syncthreads();

	if(nodeno < numnodes)
	{
		if(grouptid == 0)
		{	
			char bitsy = allowed[nodeno];

			//r = b-Ax
			float r0 =  (bitsy & 1) ? 0 : (b[nodeno * 3] - (cache[0][groupid][0] + cache[1][groupid][0]/* + cache[2][groupid][0] + cache[3][groupid][0]*/ + mass[nodeno * 3] * vt[nodeno * 3] * COEFFM));
			float r1 =  (bitsy & 2) ? 0 : (b[nodeno * 3 + 1] - (cache[0][groupid][1] + cache[1][groupid][1]/* + cache[2][groupid][1] + cache[3][groupid][1]*/ + mass[nodeno * 3 + 1] * vt[nodeno * 3 + 1] * COEFFM));
			float r2 =  (bitsy & 4) ? 0 : (b[nodeno * 3 + 2] - (cache[0][groupid][2] + cache[1][groupid][2]/* + cache[2][groupid][2] + cache[3][groupid][2]*/ + mass[nodeno * 3 + 2] * vt[nodeno * 3 + 2] * COEFFM));

			r[nodeno * 3] = r0;
			r[nodeno * 3 + 1] = r1;
			r[nodeno * 3 + 2] = r2;

			//d=r
			d[nodeno * 3] = r0;
			d[nodeno * 3 + 1] = r1;
			d[nodeno * 3 + 2] = r2;
		}
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
		mulSystem(elements, solverData, d, numelements);
}

//step 5
//CG loop
//q = Ad
__global__
void
gatherQprod(GPUNode* nodes, mulData* solverData, float* q, float* mass, float* d, char* allowed, int numnodes)
{
	int groupid = threadIdx.x % NODE_BLOCK_SIZE;// / NODE_THREADS;
	int grouptid = threadIdx.x / NODE_BLOCK_SIZE; //% NODE_THREADS;
	int nodeno = blockIdx.x * NODE_BLOCK_SIZE + groupid;

	__shared__ float cache[NODE_THREADS][NODE_BLOCK_SIZE][3];
	GPUNode* node = &(nodes[blockIdx.x]);
	int n = node->n[grouptid][groupid];
	
	if(nodeno < numnodes)
	{

		cache[grouptid][groupid][0] = 0;
		cache[grouptid][groupid][1] = 0;
		cache[grouptid][groupid][2] = 0;

		for(int i=0;i<n;i++)
		{
			int tetindex = node->elementindex[i][0][grouptid][groupid] / BLOCK_SIZE;
			int tetindex2 = node->elementindex[i][0][grouptid][groupid] % BLOCK_SIZE;
			int nodeindex = node->elementindex[i][1][grouptid][groupid];

			cache[grouptid][groupid][0] += solverData[tetindex].product[nodeindex * 3][tetindex2];
			cache[grouptid][groupid][1] += solverData[tetindex].product[nodeindex * 3 + 1][tetindex2];
			cache[grouptid][groupid][2] += solverData[tetindex].product[nodeindex * 3 + 2][tetindex2];
		}
	}

	__syncthreads();

	if(nodeno < numnodes)
	{
		if(grouptid == 0)
		{
			char bitsy = allowed[nodeno];
			q[nodeno * 3]     = (bitsy & 1) ? 0 : (cache[0][groupid][0] + cache[1][groupid][0]/* + cache[2][groupid][0] + cache[3][groupid][0]*/ + mass[nodeno * 3] * d[nodeno * 3] * COEFFM);
			q[nodeno * 3 + 1] = (bitsy & 2) ? 0 : (cache[0][groupid][1] + cache[1][groupid][1]/* + cache[2][groupid][1] + cache[3][groupid][1]*/ + mass[nodeno * 3 + 1] * d[nodeno * 3 + 1] * COEFFM);
			q[nodeno * 3 + 2] = (bitsy & 4) ? 0 : (cache[0][groupid][2] + cache[1][groupid][2]/* + cache[2][groupid][2] + cache[3][groupid][2]*/ + mass[nodeno * 3 + 2] * d[nodeno * 3 + 2] * COEFFM);
		}
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
makeXRandD(CGVars* vars, float *x, float* r, float* d, float* q, int numnodes)
{
	int tid = threadIdx.x + blockIdx.x * VECTOR_BLOCK_SIZE;
	if(tid < numnodes)
	{
		float alpha = vars->alpha;
		float beta = vars->beta;

		x[tid * 3] = x[tid * 3] + alpha * d[tid * 3];
		x[tid * 3 + 1] = x[tid * 3 + 1] + alpha * d[tid * 3 + 1];
		x[tid * 3 + 2] = x[tid * 3 + 2] + alpha * d[tid * 3 + 2];

		r[tid * 3] = r[tid * 3] - alpha * q[tid * 3];
		r[tid * 3 + 1] = r[tid * 3 + 1] - alpha * q[tid * 3 + 1];
		r[tid * 3 + 2] = r[tid * 3 + 2] - alpha * q[tid * 3 + 2];

		d[tid * 3] = r[tid * 3] + beta * d[tid * 3];
		d[tid * 3 + 1] = r[tid * 3 + 1] + beta * d[tid * 3 + 1];
		d[tid * 3 + 2] = r[tid * 3 + 2] + beta * d[tid * 3 + 2];
	}
} 

//step 8
//make x(t+1)
__global__
void
integrate(float *x, float* v, int numnodes)
{
	int tid = threadIdx.x + blockIdx.x * VECTOR_BLOCK_SIZE;
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
	const int num_blocks_node = (numnodes/NODE_BLOCK_SIZE) + 1;
	const int num_blocks_vec = (numnodes/VECTOR_BLOCK_SIZE) + 1;

	cudaError_t error;

	printf("Started\n");
	
	precompute<<<num_blocks_ele, BLOCK_SIZE>>>(gpuptrElements, gpuptrMulData, gpuptr_xt, numelements);
	
	cudaDeviceSynchronize();
	error = cudaGetLastError();
	if(error != cudaSuccess)
	{
		printf("1");
		// print the CUDA error message and exit
		printf("CUDA error: %s\n", cudaGetErrorString(error));
		//exit(-1);
	}


	gatherB<<<num_blocks_node, GATHER_THREAD_NO>>>(gpuptrNodes, gpuptrMulData, gpuptr_b, gpuptr_mass, gpuptr_vt, gpuptr_extforces, gpuptr_allowed, numnodes);

	cudaDeviceSynchronize();
	error = cudaGetLastError();
	if(error != cudaSuccess)
	{
		printf("2");
		// print the CUDA error message and exit
		printf("CUDA error: %s\n", cudaGetErrorString(error));
		//exit(-1);
	}

	initAx<<<num_blocks_ele, THREADS_PER_BLOCK>>>(gpuptrElements, gpuptrMulData, gpuptr_vt, numelements);

	cudaDeviceSynchronize();
	error = cudaGetLastError();
	if(error != cudaSuccess)
	{
		printf("3");
		// print the CUDA error message and exit
		printf("CUDA error: %s\n", cudaGetErrorString(error));
		//exit(-1);
	}

	initRandD<<<num_blocks_node, GATHER_THREAD_NO>>>(gpuptrNodes, gpuptrMulData, gpuptrR, gpuptrD, gpuptr_b, gpuptr_mass, gpuptr_vt,  gpuptr_allowed,numnodes);

	cudaDeviceSynchronize();
	error = cudaGetLastError();
	if(error != cudaSuccess)
	{
		printf("4");
		// print the CUDA error message and exit
		printf("CUDA error: %s\n", cudaGetErrorString(error));
		//exit(-1);
	}

	initDeltaVars<<<1, DOT_BLOCK_SIZE>>>(gpuptrVars, gpuptrR, numnodes);

	cudaDeviceSynchronize();
	error = cudaGetLastError();
	if(error != cudaSuccess)
	{
		printf("5");
		// print the CUDA error message and exit
		printf("CUDA error: %s\n", cudaGetErrorString(error));
		//exit(-1);
	}

	int i=0;

	CGVars vars;
	cudaMemcpy(&vars, gpuptrVars, sizeof(CGVars), cudaMemcpyDeviceToHost);

	printf("Loop Started");

	while(i < MAX_ITER && vars.deltaNew > (EPSIL * EPSIL) * vars.delta0)
	{
		makeQprod<<<num_blocks_ele, THREADS_PER_BLOCK>>>(gpuptrElements, gpuptrMulData, gpuptrD, numelements);

		cudaDeviceSynchronize();
		error = cudaGetLastError();
		if(error != cudaSuccess)
		{
			printf("6");
			// print the CUDA error message and exit
			printf("CUDA error: %s\n", cudaGetErrorString(error));
			//exit(-1);
		}

		gatherQprod<<<num_blocks_node, GATHER_THREAD_NO>>>(gpuptrNodes, gpuptrMulData, gpuptrQ, gpuptr_mass, gpuptrD, gpuptr_allowed,numnodes);

		cudaDeviceSynchronize();
		error = cudaGetLastError();
		if(error != cudaSuccess)
		{
			printf("7");
			// print the CUDA error message and exit
			printf("CUDA error: %s\n", cudaGetErrorString(error));
			//exit(-1);
		}

		makeVars<<<1, DOT_BLOCK_SIZE>>>(gpuptrVars, gpuptrD, gpuptrQ, gpuptrR, numnodes);

		cudaDeviceSynchronize();
		error = cudaGetLastError();
		if(error != cudaSuccess)
		{
			printf("8");
			// print the CUDA error message and exit
			printf("CUDA error: %s\n", cudaGetErrorString(error));
			//exit(-1);
		}

		makeXRandD<<<num_blocks_vec, VECTOR_BLOCK_SIZE>>>(gpuptrVars, gpuptr_vt, gpuptrR, gpuptrD, gpuptrQ, numnodes);

		cudaDeviceSynchronize();
		error = cudaGetLastError();
		if(error != cudaSuccess)
		{
			printf("9");
			// print the CUDA error message and exit
			printf("CUDA error: %s\n", cudaGetErrorString(error));
			//exit(-1);
		}

		cudaMemcpy(&vars, gpuptrVars, sizeof(CGVars), cudaMemcpyDeviceToHost);
		i++;

	}

	printf("Loop Ended: %d\n", i);

	integrate<<<num_blocks_vec, VECTOR_BLOCK_SIZE>>>(gpuptr_xt, gpuptr_vt, numnodes);
}


#endif