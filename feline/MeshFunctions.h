#pragma once
#include "Mesh.h"

class MeshFunctions
{
public:

static void makeLever(Mesh** m, int n);
static void makeSheet(Mesh** mesh, int n, int m);
static void makeCube(Mesh** mesh, int n, int m, int d);
};

