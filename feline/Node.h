#pragma once
#include "vector3.h"
#include "GenMatrix.h"

class Node
{
	//implements a node
public:
	float mass;
	vector3<float> pos, pos_t, vec_t,force;
	Node(vector3<float> _pos, vector3<float> _vel, vector3<float> _force, float _mass);
	Node(void);
	~Node(void);
};

