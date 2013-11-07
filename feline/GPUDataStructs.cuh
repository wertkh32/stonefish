#pragma once
#include "defines.h"

#ifdef  _GAUSSIAN_QUADRATURE_
	#define BLOCK_SIZE 32				//element block size
#else
	#define BLOCK_SIZE 32
#endif


#define DOT_BLOCK_SIZE 512			//dot product block size

#define ELEMENTS_PER_THREAD 64
#define NODE_THREADS 2
#define NODE_BLOCK_SIZE 128			//node block size
#define GATHER_THREAD_NO 256

#define SHARED_MEM_OFFSET 0

#define VECTOR_BLOCK_SIZE 512

#if defined(_BERSTEIN_POLY_)

	#define THREADS_PER_ELE 5
	#define THREADS_PER_BLOCK (BLOCK_SIZE * THREADS_PER_ELE)

#elif defined(_GAUSSIAN_QUADRATURE_)

	#define THREADS_PER_ELE 4
	#define THREADS_PER_BLOCK (BLOCK_SIZE * THREADS_PER_ELE)

#endif


//#define INDEX(x,y) (( ( x * (x+1) ) >> 1 ) + y)


struct mulData
{
	float R[3][3][BLOCK_SIZE];
	float product[NUM_NODES_PER_ELE * 3][BLOCK_SIZE];
	float b[NUM_NODES_PER_ELE * 3][BLOCK_SIZE];
};

struct CGVars
{
	float delta0, deltaNew, deltaOld, alpha, beta;
};

struct GPUElement
{
	int nodeindex[NUM_NODES_PER_ELE][BLOCK_SIZE];
	
	#if defined(_QUAD_TET_) && defined(_BERSTEIN_POLY_)
			float system[30][30][BLOCK_SIZE];
	#endif

	float B[3][3][BLOCK_SIZE]; //undefShapeMatInv ({1,2,3},{b,c,d}), ({4},{b,c,d}) = SUM(-({1,2,3},{b,c,d}))
	float c1[BLOCK_SIZE], c2[BLOCK_SIZE];

	//#ifdef _LINEAR_TET_
	//float diag[NUM_NODES_PER_ELE][3][3][BLOCK_SIZE];
	//#endif

	float f0[NUM_NODES_PER_ELE * 3][BLOCK_SIZE];

};


struct GPUNode
{
	int n[NODE_THREADS][NODE_BLOCK_SIZE];
	//{tet_index,tet_node_index}
	int elementindex[ELEMENTS_PER_THREAD][2][NODE_THREADS][NODE_BLOCK_SIZE];
};



//EXTRAS

//prototype
/*
struct GPUTinyElement
{
	int nodeindex[4][BLOCK_SIZE];
	float B[3][3][BLOCK_SIZE]; //({1,2,3},{b,c,d}), ({4},{b,c,d}) = SUM(-({1,2,3},{b,c,d}))
	float c1[BLOCK_SIZE], c2[BLOCK_SIZE];


	float x0[12][BLOCK_SIZE];
	float undefShapeMatInv[3][3][BLOCK_SIZE];
	float nodalmass[BLOCK_SIZE];
};
*/

/*
#if defined(__CUDACC__) // NVCC
   #define MY_ALIGN(n) __align__(n)
#elif defined(__GNUC__) // GCC
  #define MY_ALIGN(n) __attribute__((aligned(n)))
#elif defined(_MSC_VER) // MSVC
  #define MY_ALIGN(n) __declspec(align(n))
#else
  #error "Please provide a definition for MY_ALIGN macro for your host compiler!"
#endif
*/