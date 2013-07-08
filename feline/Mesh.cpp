#include "Mesh.h"


Mesh::Mesh(Node node_list[], int n)
{
	for(int i=0;i<n;i++)
		nodes.push(node_list+i);
	numnodes = n;

	globalStiffness = (float**)malloc(sizeof(float*) * n * 3);
	globalMass = (float**)malloc(sizeof(float*) * n * 3);
	for(int i=0;i<n * 3;i++)
	{
		globalStiffness[i] = (float*)malloc(sizeof(float) * n * 3);
		globalMass[i] = (float*)malloc(sizeof(float) * n * 3);
	}

	dt = 1.0/FPS;

	resetGlobalStiffness();
	resetGlobalMass();
}

float**
Mesh::assembleGlobalMass()
{
	resetGlobalMass();

	for(int i=0;i<elements.size();i++)
	{
		GenMatrix<float,12,12>* masse = elements[i]->getMassMat();

		for(int a=0;a<4;a++)
			for(int b=0;b<4;b++)
				for(int j=0;j<3;j++)
					for(int k=0;k<3;k++)
					{
						globalMass[nodeIndices[i][a] * 3 + j][nodeIndices[i][b] * 3 + k] += (*masse)(a + j, b + k);
					}
			
	}

	return globalMass;
}

float**
Mesh::assembleGlobalStiffness()
{
	resetGlobalStiffness();

	for(int i=0;i<elements.size();i++)
	{
		GenMatrix<float,12,12>* stiffe = elements[i]->getStiffnessMat();

		for(int a=0;a<4;a++)
			for(int b=0;b<4;b++)
				for(int j=0;j<3;j++)
					for(int k=0;k<3;k++)
					{
						globalStiffness[nodeIndices[i][a] * 3+j][nodeIndices[i][b] * 3+k] += (*stiffe)(a + j, b + k);
						//printf("%lf ",(*stiffe)(a+j,b+k));
					}
			
	}

	return globalStiffness;
}

void
Mesh::resetGlobalStiffness()
{
	for(int i=0;i<numnodes * 3;i++)
		for(int j=0;j<numnodes * 3;j++)
			globalStiffness[i][j] = 0.0;
}

void
Mesh::resetGlobalMass()
{
	for(int i=0;i<numnodes * 3;i++)
		for(int j=0;j<numnodes * 3;j++)
			globalMass[i][j] = 0.0;
}

void
Mesh::addElement(int ind0, int ind1, int ind2, int ind3, float _E, float _v, float density)
{
	nodeIndices[elements.size()][0] = ind0;
	nodeIndices[elements.size()][1] = ind1;
	nodeIndices[elements.size()][2] = ind2;
	nodeIndices[elements.size()][3] = ind3;
	elements.push(new Element(nodes[ind0], nodes[ind1], nodes[ind2], nodes[ind3], _E, _v, density));
}

/*
void
Mesh::getElastic(float* out)
{
	//element computeMatFreeVars must be called first.
	//out[i]=0;
	for(int i=0;i<elements.size();i++)
	{
		int i0 = nodeIndices[i][0] * 3,i1 = nodeIndices[i][1] * 3,i2 = nodeIndices[i][2] * 3,i3 = nodeIndices[i][3] * 3;

		Matrix3d P = elements[i]->computePiolaStressTensor();
		Matrix3d dH = P * elements[i]->getUndeformShapeMatInvT() * (-elements[i]->getVolume());
		
		out[i0] += dH(0,0); 
		out[i0 + 1] += dH(1,0); 
		out[i0 + 2] += dH(2,0);
		
		out[i1] += dH(0,1); 
		out[i1 + 1] += dH(1,1); 
		out[i1 + 2] += dH(2,1);
		
		out[i2] += dH(0,2); 
		out[i2 + 1] += dH(1,2); 
		out[i2 + 2] += dH(2,2);

		out[i3] +=     -dH(0,0)-dH(0,1)-dH(0,2);
		out[i3 + 1] += -dH(1,0)-dH(1,1)-dH(1,2);
		out[i3 + 2] += -dH(2,0)-dH(2,1)-dH(2,2);
	}
}

void
Mesh::mulRK(float* in, float* out)
{
	for(int i=0;i<elements.size();i++)
	{
		int i0 = nodeIndices[i][0] * 3,i1 = nodeIndices[i][1] * 3,i2 = nodeIndices[i][2] * 3,i3 = nodeIndices[i][3] * 3;
		
		Matrix3d dF = Matrix3d( in[i0    ] - in[i3    ],in[i1    ] - in[i3    ],in[i2    ] - in[i3    ],
								in[i0 + 1] - in[i3 + 1],in[i1 + 1] - in[i3 + 1],in[i2 + 1] - in[i3 + 1],
								in[i0 + 2] - in[i3 + 2],in[i1 + 2] - in[i3 + 2],in[i2 + 2] - in[i3 + 2]) * elements[i]->getUndeformShapeMatInv();

		Matrix3d dP = elements[i]->computeRKtensor(dF);
		Matrix3d dH = dP * elements[i]->getUndeformShapeMatInvT() * (-elements[i]->getVolume());

		out[i0] += dH(0,0); 
		out[i0 + 1] += dH(1,0); 
		out[i0 + 2] += dH(2,0);
		
		out[i1] += dH(0,1); 
		out[i1 + 1] += dH(1,1); 
		out[i1 + 2] += dH(2,1);
		
		out[i2] += dH(0,2); 
		out[i2 + 1] += dH(1,2); 
		out[i2 + 2] += dH(2,2);

		out[i3] +=     -dH(0,0)-dH(0,1)-dH(0,2);
		out[i3 + 1] += -dH(1,0)-dH(1,1)-dH(1,2);
		out[i3 + 2] += -dH(2,0)-dH(2,1)-dH(2,2);
	}
}

void
Mesh::mulRKRT(float* in, float* out)
{
	for(int i=0;i<elements.size();i++)
	{
		int i0 = nodeIndices[i][0] * 3,i1 = nodeIndices[i][1] * 3,i2 = nodeIndices[i][2] * 3,i3 = nodeIndices[i][3] * 3;
		
		Matrix3d dF = Matrix3d( in[i0    ] - in[i3    ],in[i1    ] - in[i3    ],in[i2    ] - in[i3    ],
								in[i0 + 1] - in[i3 + 1],in[i1 + 1] - in[i3 + 1],in[i2 + 1] - in[i3 + 1],
								in[i0 + 2] - in[i3 + 2],in[i1 + 2] - in[i3 + 2],in[i2 + 2] - in[i3 + 2]) * elements[i]->getUndeformShapeMatInv();
		Matrix3d dP = elements[i]->computeDiffPiolaStressTensor(dF);
		//dP.print();
		Matrix3d dH = dP * elements[i]->getUndeformShapeMatInvT() * (-elements[i]->getVolume());

		out[i0] += dH(0,0); 
		out[i0 + 1] += dH(1,0); 
		out[i0 + 2] += dH(2,0);
		
		out[i1] += dH(0,1); 
		out[i1 + 1] += dH(1,1); 
		out[i1 + 2] += dH(2,1);
		
		out[i2] += dH(0,2); 
		out[i2 + 1] += dH(1,2); 
		out[i2 + 2] += dH(2,2);

		out[i3] +=     -dH(0,0)-dH(0,1)-dH(0,2);
		out[i3 + 1] += -dH(1,0)-dH(1,1)-dH(1,2);
		out[i3 + 2] += -dH(2,0)-dH(2,1)-dH(2,2);

		//printf("%f %f %f %f\n",out[i0],out[i1],out[i2],out[i3]);
		//system("pause");
		//printf("%f %f %f\n %f %f %f \n %f %f %f\n %f %f %f\n",out[i0],out[i0+1],out[i0+2],out[i1],out[i1+1],out[i1+2],out[i2],out[i2+1],out[i2+2],out[i3],out[i3+1],out[i3+2]);
		//system("pause");
	}
	
}

void Mesh::mulA(float* in, float* out)
{
	static const float alpha = 0.1, beta = 0.3;
	static const float coeffK = dt * beta + dt * dt, coeffM = 1 + dt * alpha;
	//element computeMatFreeVars must be called first.
	//out[i]=0;
	for(int i=0;i<elements.size();i++)
	{
		int i0 = nodeIndices[i][0] * 3,i1 = nodeIndices[i][1] * 3,i2 = nodeIndices[i][2] * 3,i3 = nodeIndices[i][3] * 3;
		
		//Matrix3d dF = Matrix3d( in[i3    ] - in[i0    ],in[i3    ] - in[i1    ],in[i3    ] - in[i2    ],
		//						in[i3 + 1] - in[i0 + 1],in[i3 + 1] - in[i1 + 1],in[i3 + 1] - in[i2 + 1],
		//						in[i3 + 2] - in[i0 + 2],in[i3 + 2] - in[i1 + 2],in[i3 + 2] - in[i2 + 2]) * elements[i]->getUndeformShapeMatInv();
	
		
		Matrix3d dF = Matrix3d( in[i0    ] - in[i3    ],in[i1    ] - in[i3    ],in[i2    ] - in[i3    ],
								in[i0 + 1] - in[i3 + 1],in[i1 + 1] - in[i3 + 1],in[i2 + 1] - in[i3 + 1],
								in[i0 + 2] - in[i3 + 2],in[i1 + 2] - in[i3 + 2],in[i2 + 2] - in[i3 + 2]) * elements[i]->getUndeformShapeMatInv();
		
		Matrix3d dP = elements[i]->computeDiffPiolaStressTensor(dF);
		Matrix3d dH = dP * elements[i]->getUndeformShapeMatInvT() * (-elements[i]->getVolume());
		Matrix3d A = dH * coeffK + dF * (lumpedmass[i0] * coeffM); //assume all nodalmass of element are the same.
		
		out[i0] += A(0,0); 
		out[i0 + 1] += A(1,0); 
		out[i0 + 2] += A(2,0);
		
		out[i1] += A(0,1);            
		out[i1 + 1] += A(1,1); 
		out[i1 + 2] += A(2,1);
		
		out[i2] += A(0,2); 
		out[i2 + 1] += A(1,2); 
		out[i2 + 2] += A(2,2);

		out[i3] +=     -A(0,0)-A(0,1)-A(0,2);
		out[i3 + 1] += -A(1,0)-A(1,1)-A(1,2);
		out[i3 + 2] += -A(2,0)-A(2,1)-A(2,2);
		
		//printf("%f %f %f\n %f %f %f \n %f %f %f\n %f %f %f\n",out[i0],out[i0+1],out[i0+2],out[i1],out[i1+1],out[i1+2],out[i2],out[i2+1],out[i2+2],out[i3],out[i3+1],out[i3+2]);
	}
	//system("pause");
	//printf("\nAx\n");
	//for(int i=0;i<numnodes;i++)
		//printf("%f %f %f\n",out[i],out[i+1],out[i+2]);
	//system("pause");
}
*/
Mesh::~Mesh(void)
{
}
