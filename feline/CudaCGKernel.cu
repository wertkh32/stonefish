#pragma once
#include <cuda.h>

// conjugate gradient solver. by SKH.

#define GPU_SOLVER_MAX_ITER 1000
#define GPU_SOLVER_EPS 0.001
#define BLOCK_SIZE 256

__device__ double dot(double* a, double *b, int n)
{
	double ret = 0;
	for(int i=0;i<n;i++)
	{
		ret += a[i] * b[i];
	}
	return ret;
}

__global__ void gpuCGSolve(double** A, double* x, double* b,
							double* d, double* r, double* q,
							int n)
{

	int id = blockIdx.x *blockDim.x + threadIdx.x;
	
	if(id < n)
	{
		int i = 0;
		double alpha, beta, deltaOld, delta0, deltaNew;

		r[id] = b[id] - dot(A[id], x, n);
		d[id] = r[id];
	
		__syncthreads();

		deltaNew = dot(r,r,n);
		delta0 = deltaNew;

		while(i<GPU_SOLVER_MAX_ITER && deltaNew > GPU_SOLVER_EPS * GPU_SOLVER_EPS * delta0)
		{
			q[id] = dot(A[id],d,n);

			__syncthreads();
		
			alpha = (deltaNew)/dot(d,q,n);
			x[id] += alpha * d[id];

			__syncthreads();

			if(i%50)
			{
				r[id] = b[id] - dot(A[id], x,n);
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

void CGSolverGPU(double** A, double* x, double* b, int n)
{
	double **gpu_A, *gpu_x, *gpu_b,
			*gpu_d, *gpu_r, *gpu_q;

	int ARR_SIZE = sizeof(double) * n;

	cudaMalloc(&gpu_A, sizeof(double*) * n);
	for(int i=0;i<n;i++)
	{
		cudaMalloc(&(gpu_A[i]), ARR_SIZE);
		cudaMemcpy(gpu_A[i], A[i], ARR_SIZE, cudaMemcpyHostToDevice);
	}

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
	
	for(int i=0;i<n;i++)
	{
		cudaFree(gpu_A[i]);
	}

	cudaFree(gpu_A);
	cudaFree(gpu_x);
	cudaFree(gpu_b);
	cudaFree(gpu_d);
	cudaFree(gpu_r);

}