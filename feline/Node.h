#pragma once
#include "vector3.h"
#include "GenMatrix.h"

class Node
{
	//implements a node
public:
	float mass;
	vector3<float> pos, pos_t, force;
	Node(void);
	~Node(void);
};

