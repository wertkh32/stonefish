#pragma once

#define BLOCK_SIZE 128				//element block size
#define DOT_BLOCK_SIZE 512			//dot product block size

#define ELEMENTS_PER_THREAD 16
#define NODE_THREADS 4
#define NODE_BLOCK_SIZE 128			//node block size
#define GATHER_THREAD_NO 512

#define VECTOR_BLOCK_SIZE 512

#define INDEX(x,y) (( ( x * (x+1) ) >> 1 ) + y)


struct mulData
{
	//float system[12][12][BLOCK_SIZE];
	float R[3][3][BLOCK_SIZE];
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
	//float unwarpK[12][12][BLOCK_SIZE];
	float B[3][3][BLOCK_SIZE]; //({1,2,3},{b,c,d}), ({4},{b,c,d}) = SUM(-({1,2,3},{b,c,d}))
	float c1[BLOCK_SIZE], c2[BLOCK_SIZE];
	
	float f0[12][BLOCK_SIZE];
	float undefShapeMatInv[3][3][BLOCK_SIZE];
	float nodalmass[BLOCK_SIZE];
};


struct GPUNode
{
	int n[NODE_BLOCK_SIZE][NODE_THREADS];
	//{tet_index,tet_node_index}
	int elementindex[ELEMENTS_PER_THREAD][2][NODE_BLOCK_SIZE][NODE_THREADS];
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