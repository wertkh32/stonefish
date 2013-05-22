#include "SymSparseMatrix.h"


SparseMatrix::SparseMatrix(int _no_blocks)
{
	no_blocks = _no_blocks;
	matCollection = new Matrix3d[no_blocks * no_blocks];
	matmap = new bool[no_blocks * no_blocks];
}


SparseMatrix::~SparseMatrix(void)
{
	delete[] matCollection;
	delete[] matmap;
}

void
SparseMatrix::matProduct(float* in, float* out)
{
	for(int i=0;i<n;i++)
	{
			out[i] = 0;
	}

	for(int i=0;i<no_blocks;i++)
		for(int j=0;j<no_blocks;j++)
		{
			if(isBlockFilled(i,j))
			{
				Matrix3d& mat = getBlock(i,j);
				out[i * 3] += mat(0,0) * in[j * 3] + mat(0,1) * in[j * 3 + 1] + mat(0,2) * in[j * 3 + 2];
				out[i * 3 + 1] += mat(1,0) * in[j * 3] + mat(1,1) * in[j * 3 + 1] + mat(1,2) * in[j * 3 + 2];
				out[i * 3 + 2] += mat(2,0) * in[j * 3] + mat(2,1) * in[j * 3 + 1] + mat(2,2) * in[j * 3 + 2];
			}
		}
}
