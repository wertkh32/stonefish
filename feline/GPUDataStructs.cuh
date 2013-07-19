#pragma once

#define MAX_ELEMENTS_PER_NODE 32
#define BLOCK_SIZE 128				//element block size
#define DOT_BLOCK_SIZE 512			//dot product block size
#define NODE_BLOCK_SIZE 256			//node block size


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
	int n[NODE_BLOCK_SIZE];
	//{tet_index,tet_node_index}
	int elementindex[MAX_ELEMENTS_PER_NODE][2][NODE_BLOCK_SIZE];
};