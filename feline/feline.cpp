#include "includes.h"
#include "Integrator.h"
#include "GPUIntegrator.h"
#include "PolarDecompose.h"
#include "Model.h"
#include "MeshFunctions.h"
#include "ModelFunctions.h"
#include "perfmon.h"

//extern void CGSolverGPU(float* A, float* x, float* b, int n);
//#define DIM 100 50k
#define DIM 10
int edgemap[(DIM+1) * (DIM+1) * 2][(DIM+1) * (DIM+1) * 2] = {0};

/* GLUT callback Handlers */
static INTEGRATOR *inte;
static Mesh* tet;
static MESH* quadtet;
static int iter = 10;
Model * mod;

float tx=0, ty=0, tz=0;

static void 
resize(int width, int height)
{
    const float ar = (float) width / (float) height;
    
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-ar, ar, -1.0, 1.0, 2.0, 200.0);
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity() ;
}

bool first = true;

static void 
display(void)
{
static float rot = 0.0;

//rot += 0.3;
  iter++;
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glColor3d(1,0,0);
	glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
    glPushMatrix();
    glColor3f(1.0,0,0);
    //glTranslatef(-40 + tx,-1 + ty,-140 + tz);
	glTranslatef(-5 + tx,-3 + ty,-15 + tz);
	glRotatef(30,1,0,0);
    
	perfmon p;
	p.startTimer();
	inte->timeStep();
	p.stopTimer();
	p.print();

	//for(int i=0;i<quadtet->getNoOfElements();i++)
	//	quadtet->elements[i]->renderElement();
	mod->interpolateVerts();
	mod->render();
	//for(int i=0;i<tet->getNoOfElements();i++)
	//	tet->elements[i]->renderElement();
	//}
	
	int end = (DIM+1) * (DIM+1);
	int start = end - (DIM+1) - 1;
	int end2 = end * 2;
	int start2 = end2 - (DIM+1) - 1;
	
	if(iter>=10)
	{
		//#ifdef _LINEAR_TET_
			
			for(int i=start;i<end;i++)
				quadtet->nodes[i]->force = vector3<float>();
		
			for(int i=start2;i<end2;i++)
				quadtet->nodes[i]->force = vector3<float>();
		//#endif

		#ifdef _QUAD_TET_

			for(int i=start;i<end-1;i++)
				quadtet->nodes[edgemap[i][i+1]]->force = vector3<float>();

			for(int i=start2;i<end2-1;i++)
				quadtet->nodes[edgemap[i][i+1]]->force = vector3<float>();

			for(int i=start;i<end;i++)
				for(int j=start2;j<end2;j++)
					quadtet->nodes[edgemap[i][j]]->force = vector3<float>();

			for(int i=start;i<end-1;i++)
				for(int j=start2+1;j<end2;j++)
					quadtet->nodes[edgemap[i][j]]->force = vector3<float>();
		#endif

	}
	else
	{

		//#ifdef _LINEAR_TET_
			
			for(int i=start;i<end;i++)
				quadtet->nodes[i]->force = vector3<float>(0.001,0.001,0);
		
			for(int i=start2;i<end2;i++)
				quadtet->nodes[i]->force = vector3<float>(0.001,0.001,0);
		//#endif

		#ifdef _QUAD_TET_

			for(int i=start;i<end-1;i++)
				quadtet->nodes[edgemap[i][i+1]]->force = vector3<float>(10,10,0);

			for(int i=start2;i<end2-1;i++)
				quadtet->nodes[edgemap[i][i+1]]->force = vector3<float>(10,10,0);

			for(int i=start;i<end;i++)
				for(int j=start2;j<end2;j++)
					quadtet->nodes[edgemap[i][j]]->force = vector3<float>(10,10,0);

			for(int i=start;i<end-1;i++)
				for(int j=start2+1;j<end2;j++)
					quadtet->nodes[edgemap[i][j]]->force = vector3<float>(10,10,0);
		#endif
	
	}
	
	
	//mod->interpolateVerts();
	//mod->render();

	glPopMatrix();

    glutSwapBuffers();
}


static void 
key(unsigned char key, int x, int y)
{
    switch (key) 
    {
        case 27 : 
        case 'q':
			inte->~INTEGRATOR();
            //inte->~Integrator();
			exit(0);
            break;

		case 'w':
			tz += 1;
			break;
		case 's':
			tz -= 1;
			break;
		case 'd':
			tx -= 1;
			break;
		case 'a':
			tx += 1;
			break;
		case 'e':
			ty += 1;
			break;
		case 'r':
			ty -= 1;
			break;
		case ' ':
			iter = 0;
			break;

    }

    glutPostRedisplay();
}

static void 
idle(void)
{
    glutPostRedisplay();
}

void init_light(){
     const GLfloat light_ambient[]  = { 0.0f, 0.0f, 0.0f, 1.0f };
     const GLfloat light_diffuse[]  = { 1.0f, 1.0f, 1.0f, 1.0f };
     const GLfloat light_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
     const GLfloat light_position[] = { -2.0f, 10.0f, 5.0f, 0.0f };

     glEnable(GL_NORMALIZE);
     glEnable(GL_LIGHTING);

     glEnable(GL_LIGHT0);
     glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient);
     glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
     glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
     glLightfv(GL_LIGHT0, GL_POSITION, light_position);     
}

void init_material(){
    const GLfloat mat_ambient[]    = { 0.7f, 0.7f, 0.7f, 1.0f };
    const GLfloat mat_diffuse[]    = { 0.8f, 0.8f, 0.8f, 1.0f };
    const GLfloat mat_specular[]   = { 1.0f, 1.0f, 1.0f, 1.0f };
    const GLfloat high_shininess[] = { 100.0f };
    
    glEnable(GL_COLOR_MATERIAL);
    
    glMaterialfv(GL_FRONT, GL_AMBIENT,   mat_ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,   mat_diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR,  mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, high_shininess);
}

void init_general(){
    glClearColor(1,1,1,1);
    //glEnable(GL_CULL_FACE);
    //glCullFace(GL_BACK);
	//glFrontFace(GL_CCW);


    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);     
}

void timer(int fps) {
	
    glutPostRedisplay();
    glutTimerFunc(1000/FPS, timer, FPS);
}

/* Program entry point */

void makebox(Mesh** m)
{
	MeshFunctions::makeLever(m, 1);
}

void makerod(Mesh** m)
{
	MeshFunctions::makeLever(m, 5);
}

Mesh* loadMesh(char* nodefile, char* tetfile)
{
	FILE* nodef = fopen(nodefile,"r");
	FILE* tetf = fopen(tetfile,"r");

	int numnodes;

	fscanf(nodef,"%d %*d %*d %*d",&numnodes);

	Node* list = (Node*)malloc(numnodes * sizeof(Node));

	for(int i=0;i<numnodes;i++)
	{
		double x,y,z;
		fscanf(nodef,"%*d %lf %lf %lf",&x, &y, &z);
		list[i] = Node(vector3<float>(x,y,z),vector3<float>(),vector3<float>());
	}

	Mesh* mesh = new Mesh(list,numnodes);

	int numele;
	fscanf(tetf,"%d %*d %*d",&numele);

	for(int i=0;i<numele;i++)
	{
		int node[4];
		fscanf(tetf,"%*d %d %d %d %d",node,node+1,node+2,node+3);
		//node[0]--;
		//node[1]--;
		//node[2]--;
		//node[3]--;
		for(int j=0;j<4;j++)
			if(node[j] >= numnodes)
			{
				printf("ele %d, node: %d",i,node[j]);
				system("pause");
			}

		mesh->addElement(node,1,0.001,0.1);
	}

	fclose(nodef);
	fclose(tetf);

	return mesh;
}

QuadTetMesh* loadQuadMesh(char* nodefile, char* tetfile)
{
	FILE* nodef = fopen(nodefile,"r");
	FILE* tetf = fopen(tetfile,"r");

	int numnodes;

	fscanf(nodef,"%d %*d %*d %*d",&numnodes);

	Node* list = (Node*)malloc(numnodes * 10 * sizeof(Node));

	int** edgemat = (int**)malloc(numnodes * sizeof(int));
	for(int i=0;i<numnodes;i++)
	{
		edgemat[i] = (int*)malloc(numnodes * sizeof(int));
		
		for(int j=0;j<numnodes;j++)
			edgemat[i][j] = -1;
	}

	int realcount = numnodes;

	for(int i=0;i<numnodes;i++)
	{
		double x,y,z;
		fscanf(nodef,"%*d %lf %lf %lf",&x, &y, &z);
		list[i] = Node(vector3<float>(x,y,z),vector3<float>(),vector3<float>());
	}


	int numele;
	fscanf(tetf,"%d %*d %*d",&numele);
	for(int i=0;i<numele;i++)
	{
		int node[4];
		fscanf(tetf,"%*d %d %d %d %d",node,node+1,node+2,node+3);
		
		for(int a=0;a<4;a++)
			for(int b=a+1;b<4;b++)
			{
				if(edgemat[node[a]][node[b]] == -1)
				{
					list[realcount] = Node((list[node[a]].pos + list[node[b]].pos) * 0.5,vector3<float>(),vector3<float>());
					edgemat[node[a]][node[b]] = edgemat[node[b]][node[a]] = realcount;
					realcount++;
				}
			}
	}

	QuadTetMesh* mesh = new QuadTetMesh(list,realcount);

	fclose(tetf);
	tetf = fopen(tetfile,"r");

	fscanf(tetf,"%d %*d %*d",&numele);

	for(int i=0;i<numele;i++)
	{
		int node[10];
		fscanf(tetf,"%*d %d %d %d %d",node,node+1,node+2,node+3);
		//node[0]--;
		//node[1]--;
		//node[2]--;
		//node[3]--;
		for(int j=0;j<4;j++)
			if(node[j] >= numnodes)
			{
				printf("ele %d, node: %d",i,node[j]);
				system("pause");
			}

		node[4] = edgemat[node[0]][node[1]];
		node[5] = edgemat[node[1]][node[2]];
		node[6] = edgemat[node[2]][node[0]];
		node[7] = edgemat[node[0]][node[3]];
		node[8] = edgemat[node[1]][node[3]];
		node[9] = edgemat[node[2]][node[3]];

		mesh->addElement(node,1,0.001,0.1);
	}

	fclose(nodef);
	fclose(tetf);

	return mesh;
}

int 
main(int argc, char *argv[])
{
    glutInit(&argc, argv);
    glutInitWindowSize(640,480);
    glutInitWindowPosition(10,10);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);

    glutCreateWindow("Wertkh32's Feline Physics!");

    glutReshapeFunc(resize);
    glutDisplayFunc(display);
    glutKeyboardFunc(key);
    glutTimerFunc(1000./FPS, timer, FPS);

	//mod = new Model(ModelFunctions::rodFunc,MeshFunctions::makeQuadTetSheet);
	//mod = new Model(ModelFunctions::sphereFunc,makebox);

	ConstrainedRows rows;
	
	rows.add(0);
	rows.add(1);
	rows.add(4);
	rows.add(5);

	rows.add(0);
	rows.add(1);
	rows.add(2);
	rows.add(3);
	rows.add(16);
	rows.add(17);
	rows.add(18);
	rows.add(19);

	rows.add(20);
	rows.add(21);
	rows.add(22);
	rows.add(23);
	
	//tet = mod->mesh;

	//sheet//
	/*
	#ifdef _QUAD_TET_
		MeshFunctions::makeQuadTetSheet<DIM,DIM>(&quadtet,edgemap);
		mod = new Model(ModelFunctions::rodFunc,quadtet);
	#endif
	
	#ifdef _LINEAR_TET_
		MeshFunctions::makeSheet(&quadtet,DIM,DIM);	
		mod = new Model(ModelFunctions::rodFunc,quadtet);
	#endif
	
	for(int i=0;i<DIM+1;i++)
		rows.add(i);

	for(int i=(DIM+1) * (DIM+1);i<(DIM+1) * (DIM+1) + DIM+1;i++)
		rows.add(i);

	#ifdef _QUAD_TET_
	for(int i=0;i<DIM;i++)
		rows.add(edgemap[i][i+1]);

	for(int i=(DIM+1) * (DIM+1);i<(DIM+1) * (DIM+1) + DIM;i++)
		rows.add(edgemap[i][i+1]);
	#endif
	*/
	//MESH();
	//sheet//
	//quad tet ele stiffness test
	/*
	Node nodeset[14];
	nodeset[0] = Node(vector3<float>(0,0,0),vector3<float>(),vector3<float>());
	nodeset[1] = Node(vector3<float>(0,0,1),vector3<float>(),vector3<float>());
	nodeset[2] = Node(vector3<float>(1,0,0),vector3<float>(),vector3<float>());
	nodeset[3] = Node(vector3<float>(0,1,0),vector3<float>(),vector3<float>());
	nodeset[4] = Node((nodeset[0].pos + nodeset[1].pos) * 0.5,vector3<float>(),vector3<float>());
	nodeset[5] = Node((nodeset[1].pos + nodeset[2].pos) * 0.5,vector3<float>(),vector3<float>());
	nodeset[6] = Node((nodeset[2].pos + nodeset[0].pos) * 0.5,vector3<float>(),vector3<float>());
	nodeset[7] = Node((nodeset[0].pos + nodeset[3].pos) * 0.5,vector3<float>(),vector3<float>());
	nodeset[8] = Node((nodeset[1].pos + nodeset[3].pos) * 0.5,vector3<float>(),vector3<float>());
	nodeset[9] = Node((nodeset[2].pos + nodeset[3].pos) * 0.5,vector3<float>(),vector3<float>());

	nodeset[10] = Node(vector3<float>(0,-1,0),vector3<float>(),vector3<float>());
	nodeset[11] = Node((nodeset[0].pos + nodeset[10].pos) * 0.5,vector3<float>(),vector3<float>());
	nodeset[12] = Node((nodeset[1].pos + nodeset[10].pos) * 0.5,vector3<float>(),vector3<float>());
	nodeset[13] = Node((nodeset[2].pos + nodeset[10].pos) * 0.5,vector3<float>(),vector3<float>());
	*/
	//ELEMENT q = ELEMENT(nodeset,480,(1.0/3.0),0.5);
	//for(int i=0;i<30;i++, putchar('\n'),putchar('\n'))
	//	for(int j=0;j<30;j++)
	//		printf("%f ",q.K(i,j));
	//system("pause");
	/*
	int aaa[10] = {0,1,2,3,4,5,6,7,8,9};
	int bbb[10] = {0,1,2,10,4,5,6,11,12,13};
	quadtet = new MESH(nodeset,14);
	quadtet->addElement(aaa,5,0.1,10);
	quadtet->addElement(bbb,5,0.1,10);
	*/
	//ginte = new GPUIntegrator(tet,&rows);
	
	//quadtet = loadMesh("C:\\Users\\wertkh32\\Desktop\\felineforever\\smalldragon_nodes.txt","C:\\Users\\wertkh32\\Desktop\\felineforever\\smalldragon_tets.txt");
	
	quadtet = loadMesh("michelin07.1.node","michelin07.1.ele");
	//quadtet = loadMesh("dragon.node","dragon.ele");
	mod = new Model("michelin04_fine.ply",quadtet);
	inte = new INTEGRATOR(quadtet,&rows);
	
	//conjugate gradient test
	// works
	/*
	float **t = (float**)malloc(sizeof(float*) * 3);
	t[0] = (float*)malloc(sizeof(float) * 3);
	t[1] = (float*)malloc(sizeof(float) * 3);
	t[2] = (float*)malloc(sizeof(float) * 3);

	t[0][0] = 2.0;t[0][1] = 1.0;t[0][2] = -1.0;
	t[1][0] = -3.0;t[1][1] = -1.0;t[1][2] = 2.0;
	t[2][0] = -2.0;t[2][1] = 1.0;t[2][2] = 2.0;

	float b[3] = {8,-11,-3};
	float x[3] = {0,0,0};

	ConjugateGradientSolver cg(3,t);
	cg.solve(x,b);

	printf("%f %f %f\n", x[0],x[1],x[2]);
	*/

	//matrix3d test
	//test ok
	/*
	Matrix3d in = Matrix3d(2,0,0,
							0,2,0,
							0,0,2).inverse();

	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
			printf("%f ",in(i,j));
		printf("\n");
	}
	*/
	//matrix4d test
	//test ok
	/*
	Matrix4d in = Matrix4d(2,0,0,0,
			0,2,0,0,
			0,0,2,0,
			0,0,0,2).inverse();

	for(int i=0;i<4;i++)
	{
		for(int j=0;j<4;j++)
			printf("%f ",in(i,j));
		printf("\n");
	}*/
	//polar decompose test
	///*test ok
	//Matrix3d h(1,2,3,
	//		3,2,1,
	//		1,3,2);
	//Matrix3d R, S;

	//PolarDecompose::compute(h,R,S);

	//for(int i=0;i<3;i++)
	//{
	//	for(int j=0;j<3;j++)
	//		printf("%lf ",R(i,j));
	//	printf("\n");
	//}
	//system("pause");
	//*/

	//gpu cg solve test
	/*
	float t[] = {2,0,0,
				 0,-1,0,
				 0,0,3};

	float b[3] = {8,-11,-3};
	float x[3] = {0,0,0};

	CGSolverGPU(t,x,b,3);

	printf("%f %f %f\n", x[0],x[1],x[2]);
	*/
	


    init_general();
    init_light();
    init_material();

    glutMainLoop();

    return EXIT_SUCCESS;
}

