#pragma once
#include "Element.h"

//implements linear elasticity. For now.

class ForceEqn
{
	Element* element;
public:
	ForceEqn(Element* e);
	~ForceEqn(void);
};

