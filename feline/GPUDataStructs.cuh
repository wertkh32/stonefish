#pragma once

#define MAX_ELEMENTS_PER_NODE 32

struct mulData
{
	float system[12][12];
	float product[12];
};

struct GPUElement
{
	int nodeindex[4];
	float unwarpK[12][12];
	float x0[12];
	float b[12];
	float nodalmass;
};

struct GPUNode
{
	int n;
	//{tet_index,tet_node_index}
	int elementindex[MAX_ELEMENTS_PER_NODE][2];
};