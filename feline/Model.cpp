#include "Model.h"

Model::Model(void (*Modelfunc)(vertArray*,	edgeArray*,	faceArray*), void (*Meshfunc)(Mesh**))
	:verts(),faces(),edges()
{
	Modelfunc(&verts, &edges, &faces);
	Meshfunc(&mesh);
	computeBarycentricCoords();

}

void
Model::computeBarycentricCoords()
{
	for(int i=0;i<verts.size();i++)
	{
		int j=0;
		vector3<double> barytest;
		bool valid = false;

		for(j=0;j<mesh->elements.size();j++)
		{
			barytest = mesh->elements[j]->getUndeformShapeMatInv() * (verts[i].coords - mesh->elements[j]->getNodes()[3]->pos_t);
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

void
Model::interpolateVerts()
{
	for(int i=0;i<verts.size();i++)
	{
		if(barys[i].element_no == -1) continue;
		verts[i].coords = (mesh->elements[barys[i].element_no]->computeDeformShapeMat() * barys[i].barycoords) + mesh->elements[barys[i].element_no]->getNodes()[3]->pos_t;
	}
}

void
Model::render()
{
	glDisable(GL_CULL_FACE);
	glPushMatrix();
	glBegin(GL_QUADS);

	for(int i=0;i<faces.size();i++)
	{
		glNormal3dv( (*(faces[i].verts))[0].norm.coords );
		glVertex3dv( (*(faces[i].verts))[0].coords.coords );

		glNormal3dv( (*(faces[i].verts))[1].norm.coords );
		glVertex3dv( (*(faces[i].verts))[1].coords.coords );

		glNormal3dv( (*(faces[i].verts))[2].norm.coords );
		glVertex3dv( (*(faces[i].verts))[2].coords.coords );

		glNormal3dv( (*(faces[i].verts))[3].norm.coords );
		glVertex3dv( (*(faces[i].verts))[3].coords.coords );
	}

	glEnd();
	glPopMatrix();
	glEnable(GL_CULL_FACE);
}

Model::~Model(void)
{
}