#pragma once
#include "Matrix3d.h"

class SparseMatrix
{
	int n;
	int no_blocks;
	bool* matmap;
	Matrix3d* matCollection;

	//void swap(int& x, int& y){int temp = x; x = y; y = temp;}
public:
	SparseMatrix(int _no_blocks);
	~SparseMatrix(void);
	void matProduct(float* in, float* out);

	Matrix3d& getBlock(int x, int y)
	{
		return matCollection[x * no_blocks + y];
	}

	void setBlock(int x, int y, Matrix3d& mat)
	{
		matCollection[x * no_blocks + y] = mat;
	}

	void setBlock(int x, int y, Matrix3d& mat, bool filled)
	{
		matCollection[x * no_blocks + y] = mat;
		setBlockFilled(x,y,filled);
	}

	void addBlock(int x, int y, Matrix3d& mat)
	{
		matCollection[x * no_blocks + y] = matCollection[x * no_blocks + y] + mat;
	}

	bool isBlockFilled(int x, int y)
	{
		return matmap[x * no_blocks + y];
	}

	void setBlockFilled(int x, int y, bool val)
	{
		//skip the swap
		matmap[x * no_blocks + y] = val;
	}

	float getValue(int globalx, int globaly)
	{
		int blockx = globalx / no_blocks;
		int blocky = globaly / no_blocks;
		return getBlock(blockx, blocky)(globalx % no_blocks, globaly % no_blocks);
	}

	void setValue(int globalx, int globaly, float val)
	{
		int blockx = globalx / no_blocks;
		int blocky = globaly / no_blocks;
		getBlock(blockx, blocky)(globalx % no_blocks, globaly % no_blocks) = val;
	}

	float operator()(int x, int y)
	{
		return getValue(x,y);
	}

	float scalarMul(float val)
	{
		for(int i=0;i<no_blocks;i++)
			for(int j=0;j<no_blocks;j++)
			{
				if(isBlockFilled(i,j))
					getBlock(i,j) = getBlock(i,j) * val;
			}
	}
};

