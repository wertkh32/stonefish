#pragma once

#define BLOCK_SIZE 128				//element block size
#define DOT_BLOCK_SIZE 512			//dot product block size

#define ELEMENTS_PER_THREAD 8
#define NODE_THREADS 4
#define NODE_BLOCK_SIZE 128			//node block size
#define GATHER_THREAD_NO 512

#define VECTOR_BLOCK_SIZE 512


struct mulData
{
	float system[12][12][BLOCK_SIZE];
	float product[12][BLOCK_SIZE];
	float b[12][BLOCK_SIZE];
};

struct CGVars
{
	float delta0, deltaNew, deltaOld, alpha, beta;
};

struct GPUElement
{
	int nodeindex[4][BLOCK_SIZE];
	float unwarpK[12][12][BLOCK_SIZE];
	float x0[12][BLOCK_SIZE];
	float undefShapeMatInv[3][3][BLOCK_SIZE];
	float nodalmass[BLOCK_SIZE];
};


struct GPUNode
{
	int n[NODE_BLOCK_SIZE][4];
	//{tet_index,tet_node_index}
	int elementindex[ELEMENTS_PER_THREAD][2][NODE_BLOCK_SIZE][4];
};