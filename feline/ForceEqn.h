#pragma once
#include "TetElement.h"

//implements linear elasticity. For now.

class ForceEqn
{
	TetElement* element;
public:
	ForceEqn(TetElement* e);
	~ForceEqn(void);
};

