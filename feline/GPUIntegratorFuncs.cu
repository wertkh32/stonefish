#pragma once
#include <cuda.h>
#include <stdio.h>
#include "GPUDataStructs.cuh"

#define BLOCK_SIZE 256

GPUElement* gpuptrElements;
mat12d* gpuptrWarpK;
float*   gpuptr_x0;//const
float*   gpuptr_xt;//dynamic
float*   gpuptr_vt;//dynamic
float*	 gpuptr_extforces;//dynamic

void
gpuInitVars(int numele, int numnodes)
{
	cudaMalloc(&gpuptrElements, numele * sizeof(GPUElement));
	cudaMalloc(&gpuptrWarpK, numele * sizeof(mat12d));
	cudaMalloc(&gpuptr_x0, numnodes * 3 * sizeof(float));
	cudaMalloc(&gpuptr_xt, numnodes * 3 * sizeof(float));
	cudaMalloc(&gpuptr_vt, numnodes * 3 * sizeof(float));
	cudaMalloc(&gpuptr_extforces, numnodes * 3 * sizeof(float));
}

void
gpuUploadExtForces(float* extforces, int numnodes)
{
	cudaMemcpy(gpuptr_extforces, extforces, numnodes*3*sizeof(float),cudaMemcpyHostToDevice);
}

__device__
void zeroOut(float* f, int n)
{
	for(int i=0;i<n;i++)
	{
		f[i]=0;
	}
}

__global__
void
timestep(GPUElement* elements, float* x0, float* xt, float*vt,
		float* extforces, int numeles, int numnodes)
{
	
}

