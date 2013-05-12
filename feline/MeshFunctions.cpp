#include "MeshFunctions.h"

void MeshFunctions::makeLever(Mesh** m, int n)
{
	Node* list = (Node*)malloc(sizeof(Node) * 4 * (n+1));
	for(int i=0;i<n+1;i++)
	{
		list[i*4 + 0] = Node(vector3<float>(i,0,0),vector3<float>(),vector3<float>(0,0,0),100);
		list[i*4 + 1] = Node(vector3<float>(i,0,1),vector3<float>(),vector3<float>(0,0,0),100);
		list[i*4 + 2] = Node(vector3<float>(i,1,0),vector3<float>(),vector3<float>(0,0,0),100);
		list[i*4 + 3] = Node(vector3<float>(i,1,1),vector3<float>(),vector3<float>(0,0,0),100);
	}

	*m = new Mesh(list, 4 * (n+1));

	for(int i=0;i<n;i++)
	{
		int next = i * 4;
			(*m)->addElement(1 + next,0 + next,4 + next,2 + next,50,0.1,100);
			(*m)->addElement(1 + next,5 + next,4 + next,7 + next,50,0.1,100);
			(*m)->addElement(2 + next,3 + next,7 + next,1 + next,50,0.1,100);
			(*m)->addElement(2 + next,6 + next,7 + next,4 + next,50,0.1,100);
			(*m)->addElement(1 + next,2 + next,7 + next,4 + next,50,0.1,100);
	}
}

void MeshFunctions::makeSheet(Mesh** mesh, int n, int m)
{
	n++;
	m++;

	Node* list = (Node*)malloc(sizeof(Node) * n * m * 2);
	
	for(int k=0; k<2;k++)
		for(int i=0;i<n;i++)
			for(int j=0;j<m;j++)	
			{
				list[k * m * n + i * m + j] = Node(vector3<float>(i,k,j),vector3<float>(),vector3<float>(0,0,0),100);
			}

	*mesh = new Mesh(list, n*m*2);

	for(int i=0;i<n-1;i++)
		for(int j=0;j<m-1;j++)
		{
			//int nextbase = i * m + j;
			//int nexttop = m * n + i * m + j;
			(*mesh)->addElement(i * m + j, (i+1) * m + j, i * m + (j + 1), m * n + i * m + j, 50,0.1,100);
			(*mesh)->addElement((i+1) * m + (j+1), (i+1) * m + j, i * m + (j + 1), m * n + (i + 1) * m + (j + 1), 50,0.1,100);
			(*mesh)->addElement(m * n + i * m + j, m * n + (i+1) * m + j, m * n + (i+1) * m + (j+1), (i+1) * m + j, 50,0.1,100);
			(*mesh)->addElement(m * n + i * m + j, m * n + i * m + (j+1), m * n + (i+1) * m + (j+1), i * m + (j+1), 50,0.1,100);
			(*mesh)->addElement(m * n + i * m + j, m * n + (i+1) * m + (j+1), i * m + (j + 1),(i+1) * m + j, 50,0.1,100);
		}
}