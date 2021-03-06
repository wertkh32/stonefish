#include "MeshFunctions.h"

void MeshFunctions::makeLever(Mesh** m, int n)
{
	Node* list = (Node*)malloc(sizeof(Node) * 4 * (n+1));
	for(int i=0;i<n+1;i++)
	{
		list[i*4 + 0] = Node(vector3<float>(i,0,0),vector3<float>(),vector3<float>(0,0,0));
		list[i*4 + 1] = Node(vector3<float>(i,0,1),vector3<float>(),vector3<float>(0,0,0));
		list[i*4 + 2] = Node(vector3<float>(i,1,0),vector3<float>(),vector3<float>(0,0,0));
		list[i*4 + 3] = Node(vector3<float>(i,1,1),vector3<float>(),vector3<float>(0,0,0));
	}

	*m = new Mesh(list, 4 * (n+1));

	for(int i=0;i<n;i++)
	{
		int next = i * 4;
		int tet1[4] = {1 + next,0 + next,4 + next,2 + next};
		int tet2[4] = {1 + next,5 + next,4 + next,7 + next};
		int tet3[4] = {2 + next,3 + next,7 + next,1 + next};
		int tet4[4] = {2 + next,6 + next,7 + next,4 + next};
		int tet5[4] = {1 + next,2 + next,7 + next,4 + next};

			(*m)->addElement(tet1,50,0.1,100);
			(*m)->addElement(tet2,50,0.1,100);
			(*m)->addElement(tet3,50,0.1,100);
			(*m)->addElement(tet4,50,0.1,100);
			(*m)->addElement(tet5,50,0.1,100);
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
				list[k * m * n + i * m + j] = Node(vector3<float>(i,k,j),vector3<float>(),vector3<float>(0,0,0));
			}

	*mesh = new Mesh(list, n*m*2);

	for(int i=0;i<n-1;i++)
		for(int j=0;j<m-1;j++)
		{
			//int nextbase = i * m + j;
			//int nexttop = m * n + i * m + j;
			int tet1[4] = {i * m + j, (i+1) * m + j, i * m + (j + 1), m * n + i * m + j};
			int tet2[4] = {(i+1) * m + (j+1), (i+1) * m + j, i * m + (j + 1), m * n + (i + 1) * m + (j + 1)};
			int tet3[4] = {m * n + i * m + j, m * n + (i+1) * m + j, m * n + (i+1) * m + (j+1), (i+1) * m + j};
			int tet4[4] = {m * n + i * m + j, m * n + i * m + (j+1), m * n + (i+1) * m + (j+1), i * m + (j+1)};
			int tet5[4] = {m * n + i * m + j, m * n + (i+1) * m + (j+1), i * m + (j + 1),(i+1) * m + j};

			(*mesh)->addElement(tet1, 50,0.1,10);
			(*mesh)->addElement(tet2, 50,0.1,10);
			(*mesh)->addElement(tet3, 50,0.1,10);
			(*mesh)->addElement(tet4, 50,0.1,10);
			(*mesh)->addElement(tet5, 50,0.1,10);
		}
}


void MeshFunctions::makeCube(Mesh** mesh, int n, int m, int d)
{
	n++;
	m++;
	d++;

	Node* list = (Node*)malloc(sizeof(Node) * n * m * d);
	
	for(int k=0; k<d;k++)
		for(int i=0;i<n;i++)
			for(int j=0;j<m;j++)	
			{
				list[k * m * n + i * m + j] = Node(vector3<float>(i,k,j),vector3<float>(),vector3<float>(0,0,0));
			}

	*mesh = new Mesh(list, n*m*d);

	for(int k=0; k<d-1;k++)
		for(int i=0;i<n-1;i++)
			for(int j=0;j<m-1;j++)
			{
				int tet1[4] = {k * m * n + i * m + j, 
								k * m * n + (i+1) * m + j, 
								k * m * n + i * m + (j + 1), 
								k * m * n + m * n + i * m + j};
				int tet2[4] = {k * m * n + (i+1) * m + (j+1), 
								k * m * n + (i+1) * m + j, 
								k * m * n + i * m + (j + 1), 
								k * m * n + m * n + (i + 1) * m + (j + 1)};
				int tet3[4] = {k * m * n + m * n + i * m + j, 
								k * m * n + m * n + (i+1) * m + j, 
								k * m * n + m * n + (i+1) * m + (j+1), 
								k * m * n + (i+1) * m + j};
				int tet4[4] = {k * m * n + m * n + i * m + j, 
								k * m * n + m * n + i * m + (j+1), 
								k * m * n + m * n + (i+1) * m + (j+1), 
								k * m * n + i * m + (j+1)};
				int tet5[4] = {k * m * n + m * n + i * m + j, 
								k * m * n + m * n + (i+1) * m + (j+1), 
								k * m * n + i * m + (j + 1),
								k * m * n + (i+1) * m + j};
				//int nextbase = i * m + j;
				//int nexttop = m * n + i * m + j;
				(*mesh)->addElement(tet1, 50,0.1,10);

				(*mesh)->addElement(tet2, 50,0.1,10);

				(*mesh)->addElement(tet3, 50,0.1,10);

				(*mesh)->addElement(tet4, 50,0.1,10);

				(*mesh)->addElement(tet5, 50,0.1,10);
			}
}