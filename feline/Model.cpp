#include "Model.h"

Model::Model(void (*Modelfunc)(vertArray*,	edgeArray*,	faceArray*), void (*Meshfunc)(MESH**))
	:verts(),faces(),edges()
{
	Modelfunc(&verts, &edges, &faces);
	Meshfunc(&mesh);
	computeBarycentricCoords();

}

Model::Model(void (*Modelfunc)(vertArray*,	edgeArray*,	faceArray*), MESH* _mesh)
{
	Modelfunc(&verts, &edges, &faces);
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
			barytest = mesh->elements[j]->getUndeformShapeMatInv() * (verts[i].coords - mesh->elements[j]->nodes[3]->pos_t);
			if(barytest.x >= 0 && barytest.y >= 0 && barytest.z >= 0)
			{
					valid = true;
					break;
			}
		}
		bary temp = {j,barytest};
		
		if(!valid) 
			temp.element_no = -1;
		
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
		verts[i].coords = (mesh->elements[barys[i].element_no]->computeDeformShapeMat() * barys[i].barycoords) + mesh->elements[barys[i].element_no]->getNodes()[3]->pos_t;
		#endif
		#ifdef _QUAD_TET_
		{
			
			float s1 = barys[i].barycoords.x,s2 = barys[i].barycoords.y,s3 = barys[i].barycoords.z;
			float s4 = 1-s1-s2-s3;
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
