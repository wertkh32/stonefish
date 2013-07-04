#pragma once

class gpuMat3d
{
	public:
	float mat[3][3];
	__device__
	gpuMat3d(void)
	{
		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
			{
				mat[i][j] = 0;
			}
	}
	
	__device__
	gpuMat3d(float m[3][3])
	{
		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
				mat[i][j]=m[i][j];
	}

	__device__
	gpuMat3d(float a00,float a01,float a02,
			 float a10,float a11,float a12,
			 float a20,float a21,float a22)
	{
		mat[0][0]=a00;mat[0][1]=a01;mat[0][2]=a02;
		mat[1][0]=a10;mat[1][1]=a11;mat[1][2]=a12;
		mat[2][0]=a20;mat[2][1]=a21;mat[2][2]=a22;
	}
	
	__device__
	float operator()(int x, int y)
	{
		return mat[x][y];
	}

	__device__ 
	gpuMat3d operator*(const gpuMat3d& mm){
		float m[3][3]={0};
		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
				for(int k=0;k<3;k++)
					m[i][j]+= mat[i][k]*mm.mat[k][j];
		return gpuMat3d(m);
	}

	__device__
	void vmul(float* in, float* out)
	{
		for(int i=0;i<3;i++)
		{
			out[i]=0;
			for(int j=0;j<3;j++)
				out[i]+=mat[i][j] * in[j];
		}
	}

	__device__
	gpuMat3d transpose()
	{
		return gpuMat3d(mat[0][0],mat[1][0],mat[2][0],
						mat[0][1],mat[1][1],mat[2][1],
						mat[0][2],mat[1][2],mat[2][2]);
	}
};
