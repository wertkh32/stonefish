#include "Model.h"

Model::Model(void (*Modelfunc)(vertArray*,	edgeArray*,	faceArray*), void (*Meshfunc)(MESH**))
	:verts(),faces(),edges()
{
	Modelfunc(&verts, &edges, &faces);
	Meshfunc(&mesh);
	computeBarycentricCoords();

}

Model::Model(void (*Modelfunc)(vertArray*,	edgeArray*,	faceArray*), MESH* _mesh)
	:verts(),faces(),edges()
{
	Modelfunc(&verts, &edges, &faces);
	mesh = _mesh;
	computeBarycentricCoords();
}

Model::Model(char* filename, MESH* _mesh)
	:verts(),faces(),edges()
{
	FILE*f=fopen(filename,"r");
    vector3<float> v;
	int i=0;
	int numverts;
	int numfaces;

	fscanf(f,"%d",&numverts);
	fscanf(f,"%d",&numfaces);
	//vector3<float> v;
	for(int i=0;i<numverts;i++)
	{
		fscanf(f,"%f %f %f",&v.x, &v.y, &v.z);
		vertex vvv;
		vvv.coords = v;
		vvv.norm = vector3<float>();
		verts.push(vvv);
	}

	for(int i=0;i<numfaces;i++)
	{
		int u,v,w;
		fscanf(f,"%*d %d %d %d",&u, &v, &w);
		face f;
		f.type = TRIANGLE;
		f.vindex[0] = u;
		f.vindex[1] = v;
		f.vindex[2] = w;

		vector3<float> v1 = verts[v].coords - verts[u].coords;
		vector3<float> v2 = verts[w].coords - verts[u].coords;
		vector3<float> norm = v1.cross(v2);

		if(verts[u].norm.mag() < 0.0001)
			verts[u].norm = norm;
		else
			verts[u].norm = (verts[u].norm + norm)/2.0;

		if(verts[v].norm.mag() < 0.0001)
			verts[v].norm = norm;
		else
			verts[v].norm = (verts[v].norm + norm)/2.0;

		if(verts[w].norm.mag() < 0.0001)
			verts[w].norm = norm;
		else
			verts[w].norm = (verts[w].norm + norm)/2.0;


		faces.push(f);
	}

      fclose(f);
	  mesh = _mesh;
	  computeBarycentricCoords();
}

void
Model::computeBarycentricCoords()
{
	for(int i=0;i<verts.size();i++)
	{
		int j=0;
		vector3<float> barytest;
		bool valid = false;

		for(j=0;j<mesh->elements.size();j++)
		{
			barytest = mesh->elements[j]->getUndeformShapeMatInv() * (verts[i].coords - mesh->elements[j]->nodes[3]->pos);
			if(barytest.x >= 0 && barytest.y >= 0 && barytest.z >= 0 && (barytest.x + barytest.y + barytest.z) <= 1.0)
			{
					valid = true;
					break;
			}
		}
		bary temp = {j,barytest};
		
		if(!valid)
		{
			//printf("oh thats a problem:%d\n",i);
			//system("pause");
			temp.element_no = -1;
		}
		
		barys.push(temp);
	}

}

float Model::fact(int n)
{
	int r = 1;
	while(n) r *= (n--);
	return r;
}

void
Model::interpolateVerts()
{
	for(int i=0;i<verts.size();i++)
	{
		if(barys[i].element_no == -1) continue;
		#ifdef _LINEAR_TET_
		verts[i].coords = (mesh->elements[barys[i].element_no]->computeDeformShapeMat() * barys[i].barycoords) + mesh->elements[barys[i].element_no]->nodes[3]->pos_t;
		#endif
		#ifdef _QUAD_TET_
		{
			
			float s1 = barys[i].barycoords.x,s2 = barys[i].barycoords.y,s3 = barys[i].barycoords.z;
			float s4 = 1.0-s1-s2-s3;
			Node** nds =  (mesh->elements[barys[i].element_no]->nodes);
			#if defined(_GAUSSIAN_QUADRATURE_)
				verts[i].coords = nds[0]->pos_t * s1*(2 * s1 - 1) 
								+ nds[1]->pos_t * s2*(2 * s2 - 1)
								+ nds[2]->pos_t * s3*(2 * s3 - 1)
								+ nds[3]->pos_t * s4*(2 * s4 - 1)
								+ nds[4]->pos_t * 4 * s1 * s2
								+ nds[5]->pos_t * 4 * s2 * s3
								+ nds[6]->pos_t * 4 * s3 * s1
								+ nds[7]->pos_t * 4 * s1 * s4
								+ nds[8]->pos_t * 4 * s2 * s4
								+ nds[9]->pos_t * 4 * s3 * s4;
			#elif defined(_BERSTEIN_POLY_)
				verts[i].coords = nds[0]->pos_t * s1 * s1
								+ nds[1]->pos_t * s2 * s2
								+ nds[2]->pos_t * s3 * s3 
								+ nds[3]->pos_t * s4 * s4
								+ nds[4]->pos_t * 2 * s1 * s2
								+ nds[5]->pos_t * 2 * s2 * s3
								+ nds[6]->pos_t * 2 * s3 * s1
								+ nds[7]->pos_t * 2 * s1 * s4
								+ nds[8]->pos_t * 2 * s2 * s4
								+ nds[9]->pos_t * 2 * s3 * s4;
			#endif

		}
		#endif
	}
}

void
Model::render()
{
	glDisable(GL_CULL_FACE);
	glPushMatrix();
	for(int i=0;i<faces.size();i++)
	{
		if(faces[i].type == QUAD)
		{
			glBegin(GL_QUADS);

			glNormal3fv( verts[faces[i].vindex[0]].norm.coords );
			glVertex3fv( verts[faces[i].vindex[0]].coords.coords );

			glNormal3fv( verts[faces[i].vindex[1]].norm.coords );
			glVertex3fv( verts[faces[i].vindex[1]].coords.coords );

			glNormal3fv( verts[faces[i].vindex[2]].norm.coords );
			glVertex3fv( verts[faces[i].vindex[2]].coords.coords );

			glNormal3fv( verts[faces[i].vindex[3]].norm.coords );
			glVertex3fv( verts[faces[i].vindex[3]].coords.coords );
	
			glEnd();
		}
		else if(faces[i].type == TRIANGLE)
		{
			//glPointSize(5);
			glBegin(GL_TRIANGLES);

			glNormal3fv( verts[faces[i].vindex[0]].norm.coords );
			glVertex3fv( verts[faces[i].vindex[0]].coords.coords );

			glNormal3fv( verts[faces[i].vindex[1]].norm.coords );
			glVertex3fv( verts[faces[i].vindex[1]].coords.coords );

			glNormal3fv( verts[faces[i].vindex[2]].norm.coords );
			glVertex3fv( verts[faces[i].vindex[2]].coords.coords );
	
			glEnd();
		}
	}
	glPopMatrix();
	glEnable(GL_CULL_FACE);
}

Model::~Model(void)
{
}
