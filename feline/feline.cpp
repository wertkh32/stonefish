#include "includes.h"
#include "Integrator.h"
#include "PolarDecompose.h"
#include "Model.h"
#include "MeshFunctions.h"
#include "ModelFunctions.h"

extern void CGSolverGPU(float* A, float* x, float* b, int n);

/* GLUT callback Handlers */
static Integrator* inte;
static Mesh* tet;
static int iter = 10;
Model * mod;
/*
Node nodelist[] =
{
	Node(vector3<float>(1,1,0),vector3<float>(),vector3<float>(),10),
	Node(vector3<float>(1,0,0),vector3<float>(),vector3<float>(),10),
	Node(vector3<float>(1,0,1),vector3<float>(),vector3<float>(),10),
	Node(vector3<float>(0,0,0),vector3<float>(),vector3<float>(),10),
	//Node(vector3<float>(1,-1,0),vector3<float>(),vector3<float>(),10)
};
*/

static void 
resize(int width, int height)
{
    const float ar = (float) width / (float) height;
    
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-ar, ar, -1.0, 1.0, 2.0, 100.0);
    
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
    glTranslatef(-3,-1,-10);
	glRotatef(20 + rot,1,1,0);
    //glutSolidSphere(3,30,30);
	//inte->debug();
	//inte->debug();

	inte->timeStep();

	for(int i=0;i<tet->getNoOfElements();i++)
		tet->elements[i]->renderElement();

	if(iter>=10)
	{
		for(int i=90;i<98;i++)
		tet->nodes[i]->force = vector3<float>();
		//for(int i=2;i<8;i++)
		//tet->nodes[i]->force = vector3<float>();
	}
	else
	{
		for(int i=90;i<98;i++)
		tet->nodes[i]->force = vector3<float>(0,100,0);
		//for(int i=2;i<8;i++)
		//tet->nodes[i]->force = vector3<float>(0,40,0);
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
            exit(0);
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

	//mod = new Model(ModelFunctions::rodFunc,makerod);
	//mod = new Model(ModelFunctions::sphereFunc,makebox);

	ConstrainedRows rows;

	//rows.add(0);
	//rows.add(1);
	//rows.add(4);
	//rows.add(5);

	//rows.add(0);
	//rows.add(1);
	//rows.add(2);
	//rows.add(3);
	//rows.add(16);
	//rows.add(17);
	//rows.add(18);
	//rows.add(19);

	//rows.add(20);
	//rows.add(21);
	//rows.add(22);
	//rows.add(23);
	
	//tet = mod->mesh;

	//sheet//
	rows.add(0);
	rows.add(1);
	rows.add(2);
	rows.add(3);
	rows.add(4);
	rows.add(5);
	rows.add(6);

	rows.add(49);
	rows.add(50);
	rows.add(51);
	rows.add(52);
	rows.add(53);
	rows.add(54);
	rows.add(55);


	MeshFunctions::makeSheet(&tet,6,6);
	//sheet//


	inte = new Integrator(tet,&rows);
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
	/*test ok
	Matrix3d h(4,0,0,
			0,4,0,
			0,0,4);
	Matrix3d R, S;

	PolarDecompose::compute(h,R,S);

	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
			printf("%lf ",R(i,j));
		printf("\n");
	}
	*/

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

