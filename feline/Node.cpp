#include "Node.h"


Node::Node(void)
:pos(0,0,0), pos_t(0,0,0), force(0,0,0), mass(0.0), vec_t(0,0,0)
{
}

Node::Node(vector3<double> _pos, vector3<double> _vel, vector3<double> _force, double _mass)
:pos(_pos), pos_t(_pos), force(_force), mass(_mass), vec_t(_vel)
{
}


Node::~Node(void)
{
}
