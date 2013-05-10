#pragma once
#include <cuda.h>
#include <stdio.h>

// conjugate gradient solver. by SKH.

#define GPU_SOLVER_MAX_ITER 1000
#define GPU_SOLVER_EPS 0.001
#define BLOCK_SIZE 256

__device__ float dot(float* a, float *b, int n)
{
	float ret = 0;
	for(int i=0;i<n;i++)
	{
		ret += a[i] * b[i];
	}
	return ret;
}

__global__ void gpuCGSolve(float* A, float* x, float* b,
							float* d, float* r, float* q,
							int n)
{

	int id = blockIdx.x *blockDim.x + threadIdx.x;
	
	if(id < n)
	{
		int i = 0;
		float alpha, beta, deltaOld, delta0, deltaNew;

		r[id] = b[id] - dot(&A[id * n], x, n);
		d[id] = r[id];
	
		__syncthreads();

		deltaNew = dot(r,r,n);
		delta0 = deltaNew;

		while(i<GPU_SOLVER_MAX_ITER && deltaNew > GPU_SOLVER_EPS * GPU_SOLVER_EPS * delta0)
		{
			q[id] = dot(&A[id * n],d,n);

			__syncthreads();
		
			alpha = (deltaNew)/dot(d,q,n);
			x[id] += alpha * d[id];

			__syncthreads();

			if(i%50)
			{
				r[id] = b[id] - dot(&A[id * n], x,n);
			}
			else
			{
				r[id] -= alpha * q[id];
			}

			deltaOld = deltaNew;
			deltaNew = dot(r,r,n);
			beta = deltaNew/deltaOld;

			d[id] = r[id] + beta * d[id];

			i = i+1;

			__syncthreads();

		}
	}
}

void CGSolverGPU(float* A, float* x, float* b, int n)
{
	float *gpu_A, *gpu_x, *gpu_b,
			*gpu_d, *gpu_r, *gpu_q;

	int ARR_SIZE = sizeof(float) * n;

	//performance testing/////////////////////////
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start,0);
	//end///////////////////////////////////////////

	cudaMalloc(&gpu_A, sizeof(float) * n * n);

	cudaMemcpy(gpu_A, A, sizeof(float) * n * n, cudaMemcpyHostToDevice);


	cudaMalloc(&gpu_x, ARR_SIZE);
	cudaMalloc(&gpu_b, ARR_SIZE);
	cudaMemcpy(gpu_x, x, ARR_SIZE, cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_b, b, ARR_SIZE, cudaMemcpyHostToDevice);

	cudaMalloc(&gpu_d, ARR_SIZE);
	cudaMalloc(&gpu_r, ARR_SIZE);
	cudaMalloc(&gpu_q, ARR_SIZE);

	int no_blocks = (n + BLOCK_SIZE - 1)/BLOCK_SIZE;

	gpuCGSolve <<<no_blocks, BLOCK_SIZE>>> (gpu_A, gpu_x, gpu_b, gpu_d, gpu_r, gpu_q, n);

	cudaMemcpy(x, gpu_x, ARR_SIZE, cudaMemcpyDeviceToHost);

	cudaFree(gpu_A);
	cudaFree(gpu_x);
	cudaFree(gpu_b);
	cudaFree(gpu_d);
	cudaFree(gpu_r);

	//performance testing///////////////////////
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);

	float timing;
	cudaEventElapsedTime( &timing, start, stop );

	printf("Time taken %.4f ms\n",timing);

	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	//end//////////////////////////////////////////
}
