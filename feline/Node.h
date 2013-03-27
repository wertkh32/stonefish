#pragma once
#include "vector3.h"
#include "GenMatrix.h"

class Node
{
	//implements a node
public:
	double mass;
	vector3<double> pos, pos_t, vec_t,force;
	Node(vector3<double> _pos, vector3<double> _vel, vector3<double> _force, double _mass);
	Node(void);
	~Node(void);
};

