#include "includes.h"
#include "Integrator.h"
#include "PolarDecompose.h"

/* GLUT callback Handlers */
static Integrator* inte;
static Mesh* tet;
static int iter = 10;
/*
Node nodelist[] =
{
	Node(vector3<double>(1,1,0),vector3<double>(),vector3<double>(),10),
	Node(vector3<double>(1,0,0),vector3<double>(),vector3<double>(),10),
	Node(vector3<double>(1,0,1),vector3<double>(),vector3<double>(),10),
	Node(vector3<double>(0,0,0),vector3<double>(),vector3<double>(),10),
	//Node(vector3<double>(1,-1,0),vector3<double>(),vector3<double>(),10)
};
*/

Node nodelist[] =
{
	Node(vector3<double>(0,0,0),vector3<double>(),vector3<double>(0,0,0),100),
	Node(vector3<double>(0,0,1),vector3<double>(),vector3<double>(0,0,0),100),
	Node(vector3<double>(0,1,0),vector3<double>(),vector3<double>(0,0,0),100),
	Node(vector3<double>(0,1,1),vector3<double>(),vector3<double>(0,0,0),100),

	Node(vector3<double>(1,0,0),vector3<double>(),vector3<double>(0,0,0),100),
	Node(vector3<double>(1,0,1),vector3<double>(),vector3<double>(0,0,0),100),
	Node(vector3<double>(1,1,0),vector3<double>(),vector3<double>(0,0,0),100),
	Node(vector3<double>(1,1,1),vector3<double>(),vector3<double>(0,0,0),100),

	Node(vector3<double>(2,0,0),vector3<double>(),vector3<double>(0,0,0),100),
	Node(vector3<double>(2,0,1),vector3<double>(),vector3<double>(0,0,0),100),
	Node(vector3<double>(2,1,0),vector3<double>(),vector3<double>(0,0,0),100),
	Node(vector3<double>(2,1,1),vector3<double>(),vector3<double>(0,0,0),100),

	Node(vector3<double>(3,0,0),vector3<double>(),vector3<double>(0,0,0),100),
	Node(vector3<double>(3,0,1),vector3<double>(),vector3<double>(0,0,0),100),
	Node(vector3<double>(3,1,0),vector3<double>(),vector3<double>(0,0,0),100),
	Node(vector3<double>(3,1,1),vector3<double>(),vector3<double>(0,0,0),100),
}; 

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

static void 
display(void)
{

  iter++;
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glColor3d(1,0,0);
	glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
    glPushMatrix();
    glColor3f(1.0,0,0);
    glTranslatef(-3,-1,-10);
	glRotatef(20,1,1,0);
    //glutSolidSphere(3,30,30);
	//inte->debug();

	//inte->debug();
	inte->timeStep();

	for(int i=0;i<tet->getNoOfElements();i++)
		tet->elements[i]->renderElement();

	if(iter>=10)
	{
		for(int i=0;i<8;i++)
		tet->nodes[i]->force = vector3<double>();
	}
	else
	{
		for(int i=0;i<8;i++)
		tet->nodes[i]->force = vector3<double>(-50,50,0);
	}

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
     const GLfloat light_position[] = { 2.0f, 5.0f, 5.0f, 0.0f };

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
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);     
}

void timer(int fps) {
	
    glutPostRedisplay();
    glutTimerFunc(1000/FPS, timer, FPS);
}

/* Program entry point */

void makelever(Mesh** m, int n)
{
	Node* list = (Node*)malloc(sizeof(Node) * 4 * (n+1));
	for(int i=0;i<n+1;i++)
	{
		list[i*4 + 0] = Node(vector3<double>(i,0,0),vector3<double>(),vector3<double>(0,0,0),100);
		list[i*4 + 1] = Node(vector3<double>(i,0,1),vector3<double>(),vector3<double>(0,0,0),100);
		list[i*4 + 2] = Node(vector3<double>(i,1,0),vector3<double>(),vector3<double>(0,0,0),100);
		list[i*4 + 3] = Node(vector3<double>(i,1,1),vector3<double>(),vector3<double>(0,0,0),100);
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


	makelever(&tet,10);

	ConstrainedRows rows;
	rows.add(40);
	rows.add(41);
	rows.add(42);
	rows.add(43);

	inte = new Integrator(tet, &rows);
	
	//conjugate gradient test
	// works
	/*
	double **t = (double**)malloc(sizeof(double*) * 3);
	t[0] = (double*)malloc(sizeof(double) * 3);
	t[1] = (double*)malloc(sizeof(double) * 3);
	t[2] = (double*)malloc(sizeof(double) * 3);

	t[0][0] = 2.0;t[0][1] = 1.0;t[0][2] = -1.0;
	t[1][0] = -3.0;t[1][1] = -1.0;t[1][2] = 2.0;
	t[2][0] = -2.0;t[2][1] = 1.0;t[2][2] = 2.0;

	double b[3] = {8,-11,-3};
	double x[3] = {0,0,0};

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
    init_general();
    init_light();
    init_material();

    glutMainLoop();

    return EXIT_SUCCESS;
}

