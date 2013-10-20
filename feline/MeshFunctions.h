#pragma once
#include "Mesh.h"
#include "QuadTetMesh.h"
#define C ((N+1) * (M+1) * 2)

class MeshFunctions
{
public:

static void makeLever(Mesh** m, int n);
static void makeSheet(Mesh** mesh, int n, int m);
static void makeCube(Mesh** mesh, int n, int m, int d);


template<int N, int M>
static void makeQuadTetSheet(MESH** mesh, int edgemap[C][C]);
};


template<int N, int M>
void MeshFunctions::makeQuadTetSheet(MESH** mesh, int edgemap[C][C])
{
	int n = N + 1;
	int m = M + 1;

	int numcorners = n * m * 2;

	int c = (N+1) * (M+1) * 2;

	printf("%d\n", c * 10);

	Node* list = (Node*)malloc(sizeof(Node) * c * 10);

	//int edgemap[C][C] = {0};

	for(int k=0; k<2;k++)
		for(int i=0;i<n;i++)
			for(int j=0;j<m;j++)	
			{
				list[k * m * n + i * m + j] = Node(vector3<float>(i,k,j),vector3<float>(),vector3<float>(0,0,0));
			}

	int iter = C;

	for(int i=0;i<n-1;i++)
		for(int j=0;j<m-1;j++)
		{
			//int nextbase = i * m + j;
			//int nexttop = m * n + i * m + j;
			int tets[5][4] = {{i * m + j, (i+1) * m + j, i * m + (j + 1), m * n + i * m + j},
							{(i+1) * m + (j+1), (i+1) * m + j, i * m + (j + 1), m * n + (i + 1) * m + (j + 1)},
							{m * n + i * m + j, m * n + (i+1) * m + j, m * n + (i+1) * m + (j+1), (i+1) * m + j},
							{m * n + i * m + j, m * n + i * m + (j+1), m * n + (i+1) * m + (j+1), i * m + (j+1)},
							{m * n + i * m + j, m * n + (i+1) * m + (j+1), i * m + (j + 1),(i+1) * m + j}};

			for(int k=0;k<5;k++)
				for(int a=0;a<4;a++)
					for(int b=a+1;b<4;b++)
					{
						int s = tets[k][a], t = tets[k][b];
						if(edgemap[s][t] == 0)
						{
							
							
							list[iter] = Node((list[s].pos + list[t].pos) * 0.5,vector3<float>(),vector3<float>());
							edgemap[s][t] = iter;
							edgemap[t][s] = iter;
							iter++;
						}
					}
		}

	*mesh = new MESH(list, iter);

	for(int i=0;i<n-1;i++)
		for(int j=0;j<m-1;j++)
		{
			//int nextbase = i * m + j;
			//int nexttop = m * n + i * m + j;
			int tets[5][10] = {{i * m + j, (i+1) * m + j, i * m + (j + 1), m * n + i * m + j},
							{(i+1) * m + (j+1), (i+1) * m + j, i * m + (j + 1), m * n + (i + 1) * m + (j + 1)},
							{m * n + i * m + j, m * n + (i+1) * m + j, m * n + (i+1) * m + (j+1), (i+1) * m + j},
							{m * n + i * m + j, m * n + i * m + (j+1), m * n + (i+1) * m + (j+1), i * m + (j+1)},
							{m * n + i * m + j, m * n + (i+1) * m + (j+1), i * m + (j + 1),(i+1) * m + j}};

			//1,2,3,4, 5 = 12, 6 = 23, 7 = 31, 8 = 41, 9 = 42, 10 = 43
			for(int a=0;a<5;a++)
			{
				tets[a][4] = edgemap[tets[a][0]][tets[a][1]];
				tets[a][5] = edgemap[tets[a][1]][tets[a][2]];
				tets[a][6] = edgemap[tets[a][2]][tets[a][0]];
				tets[a][7] = edgemap[tets[a][3]][tets[a][0]];
				tets[a][8] = edgemap[tets[a][3]][tets[a][1]];
				tets[a][9] = edgemap[tets[a][3]][tets[a][2]];

				(*mesh)->addElement(tets[a], 50,0.1,10);
			}
		}
}