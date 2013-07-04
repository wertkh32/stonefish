#pragma once

struct mat12d
{
	float mat[12][12];
	
	__device__
	void makeRK(float R[3][3])
	{
		float RK[12][12];
		for(int i=0;i<4;i++)
		{
			for(int j=0;j<4;j++)
			{
				for(int a=0;a<3;a++)
				{
					for(int b=0;b<3;b++)
					{
						RK[a + i * 3][b + j * 3]=0;
						
						for(int c=0;c<3;c++)
							RK[a + i * 3][b + j * 3] += R[a][c] * mat[c + i * 3][b + j * 3];
					}
				}
			}
		}
	}

};

struct GPUElement
{
	mat12d unwarpK;
	float x0[12];
	float xt[12];
	float nodalmass;
};