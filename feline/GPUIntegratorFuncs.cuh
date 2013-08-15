#pragma once
#include "GPUDataStructs.cuh"

void gpuUploadVars(GPUElement* gpuElements, GPUNode* gpuNodes,float* xt, 
			  float* vt, float* extforces, float* mass, int numnodes, int numelements);

void gpuUploadExtForces(float* extforces, int numnodes);

void gpuDownloadVars(float* xt, int numnodes);

void gpuInitVars(int numele, int numnodes);

void gpuTimeStep(int numelements, int numnodes);

void gpuDestroyVars();