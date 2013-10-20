#pragma once
#include "Model.h"

class ModelFunctions
{
public:
	static void sphereFunc(vertArray* v,	edgeArray* e,	faceArray* f);
	static void rodFunc(vertArray* v,	edgeArray* e,	faceArray* f);
	static void sheetFunc(vertArray* v, edgeArray* e, faceArray* f);
};

