#pragma once
#include <cuda.h>
#include <stdio.h>
#include "defines.h"
#include "GPUDataStructs.cuh"
#include "GPUPolarDecompose.cu"

//#define BLOCK_SIZE 512

#define ALPHA 0.3
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

	cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);

	HANDLE_ERROR( cudaMalloc(&gpuptrElements, numblocksperele * sizeof(GPUElement)) );
	HANDLE_ERROR( cudaMalloc(&gpuptrMulData, numblocksperele * sizeof(mulData)) );
	HANDLE_ERROR( cudaMalloc(&gpuptrNodes, numblockpernode * sizeof(GPUNode)) );
	HANDLE_ERROR( cudaMalloc(&gpuptr_xt, numnodes * 3 * sizeof(float)) );
	HANDLE_ERROR( cudaMalloc(&gpuptr_vt, numnodes * 3 * sizeof(float)) );
	HANDLE_ERROR( cudaMalloc(&gpuptr_extforces, numnodes * 3 * sizeof(float)) );
	HANDLE_ERROR( cudaMalloc(&gpuptr_b, numnodes * 3 * sizeof(float)) );

	HANDLE_ERROR( cudaMalloc(&gpuptrR, numnodes * 3 * sizeof(float)) );
	HANDLE_ERROR( cudaMalloc(&gpuptrD, numnodes * 3 * sizeof(float)) );
	HANDLE_ERROR( cudaMalloc(&gpuptrQ, numnodes * 3 * sizeof(float)) );
	HANDLE_ERROR( cudaMalloc(&gpuptrVars, sizeof(CGVars)) );

	float coeffK = dt * BETA + dt * dt, coeffM = 1 + dt * ALPHA;
	float dt = 1.0/FPS;

	HANDLE_ERROR( cudaMemcpyToSymbol("COEFFK", &coeffK, sizeof(float)) );
	HANDLE_ERROR( cudaMemcpyToSymbol("COEFFM", &coeffM, sizeof(float)) );
	HANDLE_ERROR( cudaMemcpyToSymbol("dt", &dt, sizeof(float)) );

		cudaDeviceSynchronize();
		cudaError_t error = cudaGetLastError();
		if(error != cudaSuccess)
		{
			// print the CUDA error message and exit
			printf("CUDA error: %s\n", cudaGetErrorString(error));
			//exit(-1);
			system("pause");
		}

	//testing
	//float F[3][3] = { {1,2,3}, {3,2,1}, {1,3,2} };
	//float R[3][3];
	//gpuComputePolarDecomposition(F,R);

	//for(int i=0;i<3;i++, printf("\n"))
	//	for(int j=0; j<3; j++)
	//		printf("%f ", R[i][j]);
	//system("pause");
}

__host__
void
gpuUploadVars(GPUElement* gpuElements, GPUNode* gpuNodes,float* xt, 
			  float* vt, float* extforces, int numnodes, int numelements)
{
	int numblocksperele = (numelements / BLOCK_SIZE) + 1;
	int numblockpernode = (numnodes / NODE_BLOCK_SIZE) + 1;

	HANDLE_ERROR( cudaMemcpy(gpuptrElements, gpuElements, numblocksperele * sizeof(GPUElement), cudaMemcpyHostToDevice) );
	HANDLE_ERROR( cudaMemcpy(gpuptrNodes, gpuNodes, numblockpernode * sizeof(GPUNode), cudaMemcpyHostToDevice) );
	HANDLE_ERROR( cudaMemcpy(gpuptr_xt, xt, numnodes * 3 * sizeof(float), cudaMemcpyHostToDevice) );
	HANDLE_ERROR( cudaMemcpy(gpuptr_vt, vt, numnodes * 3 * sizeof(float), cudaMemcpyHostToDevice) );
	HANDLE_ERROR( cudaMemcpy(gpuptr_extforces, extforces, numnodes * 3 * sizeof(float), cudaMemcpyHostToDevice) );

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
	cudaFree(gpuptr_b);
	cudaFree(gpuptrR);
	cudaFree(gpuptrD);
	cudaFree(gpuptrQ);
	cudaFree(gpuptrVars);
}

__device__
void makeFU(float x0[12][BLOCK_SIZE], float mat[12][12], float R[3][3], float out[12])
{
	int ltid = threadIdx.x;
	float kx[12];
	float x[12];

	#pragma unroll 12
	for(int i=0;i<12;i++)
		x[i] = x0[i][ltid];
	
	for(int i=0;i<12;i++)
	{
		kx[i] = 0;
		for(int j=0;j<12;j++)
			kx[i] += mat[i][j] * x[j];
	}

	for(int i=0;i<4;i++)
		for(int j=0;j<3;j++)
		{
			out[i*3 + j] = 0;
			for(int k=0;k<3;k++)
			out[i*3+j] += R[j][k] * kx[i*3 + k];
		}		
				
}

__device__
void makeRKRT(float mat[12][12], float R[3][3])
{
	float temp[12][12];
	//int ltid = threadIdx.x;

	for(int i=0;i<4;i++)
	{
		for(int j=i;j<4;j++)
		{
			for(int a=0;a<3;a++)
				for(int b=0;b<3;b++)
				{
					temp[a + i * 3][b + j * 3]= 0.0f;
						
					for(int c=0;c<3;c++)
						temp[a + i * 3][b + j * 3] += R[a][c] * mat[c + i * 3][b + j * 3];
				}
		}
	}

	for(int i=0;i<4;i++)
	{
		for(int j=i;j<4;j++)
		{
			for(int a=0;a<3;a++)
				for(int b=0;b<3;b++)
				{
					mat[a + i * 3][b + j * 3] = 0.0f;
					for(int c=0;c<3;c++)
						mat[a + i * 3][b + j * 3] += temp[a + i * 3][c + j * 3] * R[b][c]; //R(c, b) but its RT so, R(b , c)
				}
		}
	}



	//lower triangle
	for(int i=1;i<4;i++)
		for(int j=0;j<i;j++)
		{
			for(int a=0;a<3;a++)
				for(int b=0;b<3;b++)
						mat[a + i * 3][b + j * 3] = mat[b + j * 3][a + i * 3];
		}

}

__device__
void mulSystem(GPUElement* elements, mulData* solverData, float* x)
{
	int bid = blockIdx.x;
	int ltid = threadIdx.x;

	GPUElement* t_ele = &(elements[bid]);
	mulData* t_solvedata = &(solverData[bid]);

	float nodes[12];

	#pragma unroll 4
	for(int i=0;i<4;i++)
	{
		int index = t_ele->nodeindex[i][ltid];
		nodes[i * 3] = x[index * 3];
		nodes[i * 3 + 1] = x[index * 3 + 1];
		nodes[i * 3 + 2] = x[index * 3 + 2];
	}

	for(int i=0;i<12;i++)
	{
		float temp = 0;
		#pragma unroll 12
		for(int j=0;j<12;j++)
			 temp += t_solvedata->system[i][j][ltid] * nodes[j];
		t_solvedata->product[i][ltid] = temp;
	}
}

__device__
void mulSystemGather(GPUNode* nodes, mulData* solverData, float* x, int numnodes)
{
	int groupid = threadIdx.x / NODE_THREADS;
	int grouptid = threadIdx.x % NODE_THREADS;
	int nodeno = blockIdx.x * NODE_BLOCK_SIZE + groupid;

	__shared__ float cache[NODE_BLOCK_SIZE][NODE_THREADS][3];
	GPUNode* node = &(nodes[blockIdx.x]);
	int n = node->n[groupid][grouptid];
	
	if(nodeno < numnodes)
	{

		cache[groupid][grouptid][0] = 0;
		cache[groupid][grouptid][1] = 0;
		cache[groupid][grouptid][2] = 0;

		for(int i=0;i<n;i++)
		{
			int tetindex = node->elementindex[i][0][groupid][grouptid] / BLOCK_SIZE;
			int tetindex2 = node->elementindex[i][0][groupid][grouptid] % BLOCK_SIZE;
			int nodeindex = node->elementindex[i][1][groupid][grouptid];

			cache[groupid][grouptid][0] += solverData[tetindex].product[nodeindex * 3][tetindex2];
			cache[groupid][grouptid][1] += solverData[tetindex].product[nodeindex * 3 + 1][tetindex2];
			cache[groupid][grouptid][2] += solverData[tetindex].product[nodeindex * 3 + 2][tetindex2];
		}
	}

	__syncthreads();

	if(nodeno < numnodes)
	{
		if(grouptid == 0)
		{
			x[nodeno * 3]     = cache[groupid][0][0] + cache[groupid][1][0] + cache[groupid][2][0] + cache[groupid][3][0];
			x[nodeno * 3 + 1] = cache[groupid][0][1] + cache[groupid][1][1] + cache[groupid][2][1] + cache[groupid][3][1];
			x[nodeno * 3 + 2] = cache[groupid][0][2] + cache[groupid][1][2] + cache[groupid][2][2] + cache[groupid][3][2];
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
void precompute(GPUElement* elements, mulData* solverData, float* xt, float* vt, float* extforces, int numelements)
{
	int tid = threadIdx.x + blockIdx.x * BLOCK_SIZE;
	int bid = blockIdx.x;
	int ltid = threadIdx.x;

	if(tid < numelements)
	{
		GPUElement* t_ele = &(elements[bid]);
		mulData* t_solvedata = &(solverData[bid]);

		float nodalmass = t_ele->nodalmass[ltid];

		float nodes[12], b[12], R[3][3]={0};// = {{1,0,0},{0,1,0},{0,0,1}};
		float K[12][12];
		int index[4];

		#pragma unroll 4
		for(int i=0;i<4;i++)
			index[i] = t_ele->nodeindex[i][ltid];

		#pragma unroll 4
		for(int i=0;i<4;i++)
		{
			nodes[i * 3] = xt[index[i] * 3];
			nodes[i * 3 + 1] = xt[index[i] * 3 + 1];
			nodes[i * 3 + 2] = xt[index[i] * 3 + 2];
		}

		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
				for(int k=0;k<3;k++)
					R[i][j] += (nodes[k*3 + i] - nodes[9 + i]) * t_ele->undefShapeMatInv[k][j][ltid];

		gpuComputePolarDecomposition(R,R);


		for(int i=0;i<12;i++)
			for(int j=0;j<12;j++)
			{
				K[i][j] = t_ele->unwarpK[i][j][ltid];
			}
	
		makeFU(t_ele->x0,K,R,b);
	
		makeRKRT(K, R);


		for(int i=0;i<12;i++)
			for(int j=0;j<12;j++)
				b[i] -= K[i][j] * nodes[j];

		#pragma unroll 4
		for(int i=0;i<4;i++)
		{
			b[i * 3] += extforces[index[i] * 3];
			b[i * 3 + 1] += extforces[index[i] * 3 + 1];
			b[i * 3 + 2] += extforces[index[i] * 3 + 2];
		}

		//final system matrix
		for(int i=0;i<12;i++)
			for(int j=0;j<12;j++)
			{
				K[i][j] *= COEFFK;
				if(i==j)
					K[i][i] += COEFFM * nodalmass;
			}

		#pragma unroll 12
		for(int i=0;i<12;i++)
			t_solvedata->b[i][ltid] = b[i] * dt + nodalmass * vt[t_ele->nodeindex[(i/3)][ltid] * 3 + (i%3)];

		for(int i=0;i<12;i++)
			for(int j=0;j<12;j++)
				t_solvedata->system[i][j][ltid] =  K[i][j];
	}
}

//step 2
//precompute
__global__
void gatherB(GPUNode* nodes, mulData* solverData, float* b, int numnodes)
{
	int groupid = threadIdx.x / NODE_THREADS;
	int grouptid = threadIdx.x % NODE_THREADS;
	int nodeno = blockIdx.x * NODE_BLOCK_SIZE + groupid;

	__shared__ float cache[NODE_BLOCK_SIZE][NODE_THREADS][3];
	GPUNode* node = &(nodes[blockIdx.x]);
	int n = node->n[groupid][grouptid];
	
	if(nodeno < numnodes)
	{

		cache[groupid][grouptid][0] = 0;
		cache[groupid][grouptid][1] = 0;
		cache[groupid][grouptid][2] = 0;

		for(int i=0;i<n;i++)
		{
			int tetindex = node->elementindex[i][0][groupid][grouptid] / BLOCK_SIZE;
			int tetindex2 = node->elementindex[i][0][groupid][grouptid] % BLOCK_SIZE;
			int nodeindex = node->elementindex[i][1][groupid][grouptid];

			cache[groupid][grouptid][0] += solverData[tetindex].b[nodeindex * 3][tetindex2];
			cache[groupid][grouptid][1] += solverData[tetindex].b[nodeindex * 3 + 1][tetindex2];
			cache[groupid][grouptid][2] += solverData[tetindex].b[nodeindex * 3 + 2][tetindex2];
		}
	}

	__syncthreads();

	if(nodeno < numnodes)
	{
		if(grouptid == 0)
		{
			b[nodeno * 3]     = cache[groupid][0][0] + cache[groupid][1][0] + cache[groupid][2][0] + cache[groupid][3][0];
			b[nodeno * 3 + 1] = cache[groupid][0][1] + cache[groupid][1][1] + cache[groupid][2][1] + cache[groupid][3][1];
			b[nodeno * 3 + 2] = cache[groupid][0][2] + cache[groupid][1][2] + cache[groupid][2][2] + cache[groupid][3][2];
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
	int groupid = threadIdx.x / NODE_THREADS;
	int grouptid = threadIdx.x % NODE_THREADS;
	int nodeno = blockIdx.x * NODE_BLOCK_SIZE + groupid;

	__shared__ float cache[NODE_BLOCK_SIZE][NODE_THREADS][3];
	GPUNode* node = &(nodes[blockIdx.x]);
	int n = node->n[groupid][grouptid];
	
	if(nodeno < numnodes)
	{

		cache[groupid][grouptid][0] = 0;
		cache[groupid][grouptid][1] = 0;
		cache[groupid][grouptid][2] = 0;

		for(int i=0;i<n;i++)
		{
			int tetindex = node->elementindex[i][0][groupid][grouptid] / BLOCK_SIZE;
			int tetindex2 = node->elementindex[i][0][groupid][grouptid] % BLOCK_SIZE;
			int nodeindex = node->elementindex[i][1][groupid][grouptid];

			cache[groupid][grouptid][0] += solverData[tetindex].product[nodeindex * 3][tetindex2];
			cache[groupid][grouptid][1] += solverData[tetindex].product[nodeindex * 3 + 1][tetindex2];
			cache[groupid][grouptid][2] += solverData[tetindex].product[nodeindex * 3 + 2][tetindex2];
		}
	}
	
	__syncthreads();

	if(nodeno < numnodes)
	{
		if(grouptid == 0)
		{
			r[nodeno * 3]     = cache[groupid][0][0] + cache[groupid][1][0] + cache[groupid][2][0] + cache[groupid][3][0];
			r[nodeno * 3 + 1] = cache[groupid][0][1] + cache[groupid][1][1] + cache[groupid][2][1] + cache[groupid][3][1];
			r[nodeno * 3 + 2] = cache[groupid][0][2] + cache[groupid][1][2] + cache[groupid][2][2] + cache[groupid][3][2];
		
			//r = b-Ax
			r[nodeno * 3] = b[nodeno * 3] - r[nodeno * 3];
			r[nodeno * 3 + 1] = b[nodeno * 3 + 1] - r[nodeno * 3 + 1];
			r[nodeno * 3 + 2] = b[nodeno * 3 + 2] - r[nodeno * 3 + 2];

			//d=r
			d[nodeno * 3] = r[nodeno * 3];
			d[nodeno * 3 + 1] = r[nodeno * 3 + 1];
			d[nodeno * 3 + 2] = r[nodeno * 3 + 2];
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
		mulSystemGather(nodes, solverData, q, numnodes);

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
	
	precompute<<<num_blocks_ele, BLOCK_SIZE>>>(gpuptrElements, gpuptrMulData, gpuptr_xt, gpuptr_vt, gpuptr_extforces, numelements);
	
	cudaDeviceSynchronize();
	error = cudaGetLastError();
	if(error != cudaSuccess)
	{
		printf("1");
		// print the CUDA error message and exit
		printf("CUDA error: %s\n", cudaGetErrorString(error));
		//exit(-1);
	}


	gatherB<<<num_blocks_node, GATHER_THREAD_NO>>>(gpuptrNodes, gpuptrMulData, gpuptr_b, numnodes);

	cudaDeviceSynchronize();
	error = cudaGetLastError();
	if(error != cudaSuccess)
	{
		printf("2");
		// print the CUDA error message and exit
		printf("CUDA error: %s\n", cudaGetErrorString(error));
		//exit(-1);
	}

	initAx<<<num_blocks_ele, BLOCK_SIZE>>>(gpuptrElements, gpuptrMulData, gpuptr_vt, numelements);

	cudaDeviceSynchronize();
	error = cudaGetLastError();
	if(error != cudaSuccess)
	{
		printf("3");
		// print the CUDA error message and exit
		printf("CUDA error: %s\n", cudaGetErrorString(error));
		//exit(-1);
	}

	initRandD<<<num_blocks_node, GATHER_THREAD_NO>>>(gpuptrNodes, gpuptrMulData, gpuptrR, gpuptrD, gpuptr_b, numnodes);

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
		makeQprod<<<num_blocks_ele, BLOCK_SIZE>>>(gpuptrElements, gpuptrMulData, gpuptrD, numelements);

		cudaDeviceSynchronize();
		error = cudaGetLastError();
		if(error != cudaSuccess)
		{
			printf("6");
			// print the CUDA error message and exit
			printf("CUDA error: %s\n", cudaGetErrorString(error));
			//exit(-1);
		}

		gatherQprod<<<num_blocks_node, GATHER_THREAD_NO>>>(gpuptrNodes, gpuptrMulData, gpuptrQ, numnodes);

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



















//EXTRAS

//completely mat free version.

/*
__device__
void mulK(float x[12], float B[3][3][BLOCK_SIZE], float c1[BLOCK_SIZE], float c2[BLOCK_SIZE])
{
	int ltid = threadIdx.x;
	float temp[6];
	float temp2[6];
	float b[3][3];

	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			b[i][j] = B[i][j][ltid];
	
	float con1 = c1[ltid];
	float con2 = c2[ltid];
	float con3 = (con1 - con2)/2.0;

	float b4 = -b[0][0] -b[1][0] -b[2][0]; 
	float c4 = -b[0][1] -b[1][1] -b[2][1]; 
	float d4 = -b[0][2] -b[1][2] -b[2][2]; 

	temp[0] = b[0][0] * x[0] + b[1][0] * x[3] + b[2][0] * x[6] + b4 * x[9];
	temp[1] = b[0][1] * x[1] + b[1][1] * x[4] + b[2][1] * x[7] + c4 * x[10];
	temp[2] = b[0][2] * x[2] + b[1][2] * x[5] + b[2][2] * x[8] + d4 * x[11];
	temp[3] = b[0][1] * x[0] + b[0][0] * x[1] + b[1][1] * x[3] + b[1][0] * x[4] + b[2][1] * x[6] + b[2][0] * x[7] + c4 * x[9] + b4 * x[10];
	temp[4] = b[0][2] * x[1] + b[0][1] * x[2] + b[1][2] * x[4] + b[1][1] * x[5] + b[2][2] * x[7] + b[2][1] * x[8] + d4 * x[10] + c4 * x[11];
	temp[5] = b[0][2] * x[0] + b[0][0] * x[2] + b[1][2] * x[3] + b[1][0] * x[5] + b[2][2] * x[6] + b[2][0] * x[8] + d4 * x[9] + b4 * x[11];

	temp2[0] = temp[0] * con1 + temp[1] * con2 + temp[2] * con2;
	temp2[1] = temp[0] * con2 + temp[1] * con1 + temp[2] * con2;
	temp2[2] = temp[0] * con2 + temp[1] * con2 + temp[2] * con1;
	temp2[3] = temp[3] * con3;
	temp2[4] = temp[4] * con3;
	temp2[5] = temp[5] * con3;

	x[0] = b[0][0] * temp2[0] + b[0][1] * temp2[3] + b[0][2] * temp2[5];
	x[1] = b[0][1] * temp2[1] + b[0][0] * temp2[3] + b[0][2] * temp2[4];
	x[2] = b[0][2] * temp2[2] + b[0][1] * temp2[4] + b[0][0] * temp2[5];

	x[3] = b[1][0] * temp2[0] + b[1][1] * temp2[3] + b[1][2] * temp2[5];
	x[4] = b[1][1] * temp2[1] + b[1][0] * temp2[3] + b[1][2] * temp2[4];
	x[5] = b[1][2] * temp2[2] + b[1][1] * temp2[4] + b[1][0] * temp2[5];

	x[6] = b[2][0] * temp2[0] + b[2][1] * temp2[3] + b[2][2] * temp2[5];
	x[7] = b[2][1] * temp2[1] + b[2][0] * temp2[3] + b[2][2] * temp2[4];
	x[8] = b[2][2] * temp2[2] + b[2][1] * temp2[4] + b[2][0] * temp2[5];

	x[9] = b4 * temp2[0] + c4 * temp2[3] + d4 * temp2[5];
	x[10] = c4 * temp2[1] + b4 * temp2[3] + d4 * temp2[4];
	x[11] = d4 * temp2[2] + c4 * temp2[4] + b4 * temp2[5];
}
*/