#pragma once
#include "includes.h"
#include "Element.h"
#include "QuickArray.h"

class Mesh
{
	
public:
	QuickArray<Element*,MAX_ELEMENT> elements;
	Mesh(void);
	~Mesh(void);
};

