#include "ModelFunctions.h"

void ModelFunctions::sphereFunc(vertArray* v,	edgeArray* e,	faceArray* f)
{
	double r = 0.499;
	int n = 8;
	double tx = 0.5, ty = 0.5, tz = 0.5; 
	for(int i=0;i<n;i++)
		for(int j=0;j<2*n;j++)
		{

                vector3<double> v1 = vector3<double>(r*sin(i*M_PI/n)*cos(j*M_PI/n) + tx,              r*cos(i*M_PI/n)*cos(j*M_PI/n) + ty,       r*sin(j*M_PI/n) + tz);

                vector3<double> v2 = vector3<double>(r*sin((i+1)*M_PI/n)*cos(j*M_PI/n) + tx,          r*cos((i+1)*M_PI/n)*cos(j*M_PI/n) + ty,    r*sin(j*M_PI/n) + tz);

                vector3<double> v3 = vector3<double>(r*sin((i+1)*M_PI/n)*cos((j+1)*M_PI/n) + tx,      r*cos((i+1)*M_PI/n)*cos((j+1)*M_PI/n) + ty,    r*sin((j+1)*M_PI/n) + tz);

                vector3<double> v4 = vector3<double>(r*sin(i*M_PI/n)*cos((j+1)*M_PI/n)+ tx,           r*cos(i*M_PI/n)*cos((j+1)*M_PI/n) + ty,       r*sin((j+1)*M_PI/n) + tz);

				vertex vv1 = {v1,v1};
				vertex vv2 = {v2,v2};
				vertex vv3 = {v3,v3};
				vertex vv4 = {v4,v4};

				v->push(vv1);
				v->push(vv2);
				v->push(vv3);
				v->push(vv4);
		}

	for(int i=0;i<v->size();i+=4)
	{
		face f1 = {QUAD,i,i+1,i+2,i+3};
		f->push(f1);
	}
}


void ModelFunctions::rodFunc(vertArray* v,	edgeArray* e,	faceArray* f)
{
	int len = 5;
	int len_quad_no = len * 3;
	double quad_len = (double)len/(double)len_quad_no;
	double EPS = 0.01;
	double r = 0.499;
	int n = 10;

	for(int i=0;i<n;i++)
	{
		double angle = ((M_PI * 2.0)/n) * i;
		double next_angle = ((M_PI * 2.0)/n) * (i+1);
		vector3<double> v1 = vector3<double>(EPS,0.5,0.5);
		vector3<double> v2 = vector3<double>(EPS,0.5 + r * cos(angle),0.5 + r * sin(angle));
		vector3<double> v3 = vector3<double>(EPS,0.5 + r * cos(next_angle),0.5 + r * sin(next_angle));
		vector3<double> norm = vector3<double>(-1,0,0);
		vertex vv1 = {v1,norm};
		vertex vv2 = {v2,norm};
		vertex vv3 = {v3,norm};
		v->push(vv1);
		v->push(vv2);
		v->push(vv3);

		
	}

	for(int i=0;i<n;i++)
	{
		double angle = ((M_PI * 2.0)/n) * i;
		double next_angle = ((M_PI * 2.0)/n) * (i+1);
		vector3<double> v1 = vector3<double>(len-EPS,0.5,0.5);
		vector3<double> v2 = vector3<double>(len-EPS,0.5 + r * cos(angle),0.5 + r * sin(angle));
		vector3<double> v3 = vector3<double>(len-EPS,0.5 + r * cos(next_angle),0.5 + r * sin(next_angle));
		vector3<double> norm = vector3<double>(1,0,0);
		vertex vv1 = {v1,norm};
		vertex vv2 = {v2,norm};
		vertex vv3 = {v3,norm};
		v->push(vv1);
		v->push(vv2);
		v->push(vv3);
	}
	
	for(int i=0;i<v->size();i+=3)
	{
		face f1 = {TRIANGLE,i,i+1,i+2,0};
		f->push(f1);
	}

	int body_start = v->size();

	for(int i=0;i<len_quad_no;i++)
	{
		for(int j=0;j<n;j++)
		{
			double angle = ((M_PI * 2.0)/n) * j;
			double next_angle = ((M_PI * 2.0)/n) * (j+1);

			vector3<double> norm = vector3<double>(0,0.5 + r * cos(angle),0.5 + r * sin(angle));
			vector3<double> next_norm = vector3<double>(0,0.5 + r * cos(next_angle),0.5 + r * sin(next_angle));

			vector3<double> v1 = vector3<double>(EPS + quad_len * i,0.5 + r * cos(angle),0.5 + r * sin(angle));
			vector3<double> v2 = vector3<double>(EPS + quad_len * i,0.5 + r * cos(next_angle),0.5 + r * sin(next_angle));
			vector3<double> v3 = vector3<double>(EPS + quad_len * (i+1),0.5 + r * cos(next_angle),0.5 + r * sin(next_angle));
			vector3<double> v4 = vector3<double>(EPS + quad_len * (i+1),0.5 + r * cos(angle),0.5 + r * sin(angle));

			vertex vv1 = {v1,norm};
			vertex vv2 = {v2,next_norm};
			vertex vv3 = {v3,next_norm};
			vertex vv4 = {v4,norm};
		
			v->push(vv1);
			v->push(vv2);
			v->push(vv3);
			v->push(vv4);
		}
	}

	for(int i=body_start;i<v->size();i+=4)
	{
		face f1 = {QUAD,i, i+1, i+2, i+3};
		f->push(f1);
	}
}