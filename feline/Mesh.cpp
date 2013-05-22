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
	
	//uncomment to activate
	//F_t.push(Matrix3d());
	//dR_t.push(Matrix3d());
}

void
Mesh::MatFree_TimestepPrecomp()
{
	for(int i=0;i<elements.size();i++)
	{
		Matrix3d R,S;
		F_t[i] = elements[i]->computeDeformationMat();
		PolarDecompose::compute(F_t[i],R,S);
		dR_t[i] = R - dR_t[i];
	}
}

void
Mesh::MatFree_stressTensor(Matrix3d& dF)
{
	//to be continued
}

void
Mesh::MatFree_stiffnessProduct(float* x, float* out)
{
	for(int i=0;i<numnodes * 3;i++)
	{
		out[i] = 0;
	}

	for(int i=0;i<elements.size();i++)
	{
		Matrix3d F = elements[i]->computeDeformationMat();
		int i0 = nodeIndices[i][0], i1 = nodeIndices[i][1], i2 = nodeIndices[i][2], i3 = nodeIndices[i][3];
		
		Matrix3d dD(x[i0 * 3]	  - x[i3 * 3],	   x[i1 * 3]	 - x[i3 * 3],	  x[i2 * 3]		- x[i3 * 3],
					x[i0 * 3 + 1] - x[i3 * 3 + 1], x[i1 * 3 + 1] - x[i3 * 3 + 1], x[i2 * 3 + 1] - x[i3 * 3 + 1],
					x[i0 * 3 + 2] - x[i3 * 3 + 2], x[i1 * 3 + 2] - x[i3 * 3 + 2], x[i2 * 3 + 2] - x[i3 * 3 + 2]);
		
		Matrix3d shapeInv = elements[i]->getUndeformShapeMatInv();
		
		Matrix3d dF = dD * shapeInv;

		float miu = elements[i]->getYoungMod(), lambda = elements[i]->getPoissonRatio();

		//Matrix3d dP = 2 * miu * (dF - dR
	}
}

Mesh::~Mesh(void)
{
}
