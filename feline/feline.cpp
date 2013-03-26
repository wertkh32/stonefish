#include "includes.h"
#include "Integrator.h"

/* GLUT callback Handlers */
static Integrator* inte;
static Mesh* tet;
Node nodelist[] =
{
	Node(vector3<float>(1,1,0),vector3<float>(),vector3<float>(),10),
	Node(vector3<float>(1,0,0),vector3<float>(),vector3<float>(0,2,0),10),
	Node(vector3<float>(1,0,1),vector3<float>(),vector3<float>(),10),
	Node(vector3<float>(0,0,0),vector3<float>(),vector3<float>(),10),
	Node(vector3<float>(1,-1,0),vector3<float>(),vector3<float>(),10)
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
  
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glColor3d(1,0,0);
	glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
    glPushMatrix();
    glColor3f(1.0,0,0);
    glTranslatef(0,0,-5);
	glRotatef(10,1,1,0);
    //glutSolidSphere(3,30,30);
	inte->timeStep();
	
	for(int i=0;i<tet->getNoOfElements();i++)
		tet->elements[i]->renderElement();

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
    glutTimerFunc(1000/FPS, timer, FPS);

	tet = new Mesh(nodelist,5);
	tet->addElement(0,1,2,3,0.1,0.3,10);
	tet->addElement(4,1,2,3,0.1,0.3,10);
	inte = new Integrator(tet);
	
	//conjugate gradient test
	/* works
	float **t = (float**)malloc(sizeof(float*) * 2);
	t[0] = (float*)malloc(sizeof(float) * 2);
	t[1] = (float*)malloc(sizeof(float) * 2);

	t[0][0] = 1.0;
	t[0][1] = 2.0;

	t[1][0] = 2.0;
	t[1][1] = 1.0;

	float b[2] = {8,7};
	float x[2] = {0,0};

	ConjugateGradientSolver cg(2,t);
	cg.solve(x,b);

	printf("%f %f\n", x[0],x[1]);
	*/

    init_general();
    init_light();
    init_material();

    glutMainLoop();

    return EXIT_SUCCESS;
}

