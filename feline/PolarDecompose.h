#pragma once
#include "includes.h"
#include "Matrix3d.h"
#include "vector3.h"
#define TOLERANCE 0.0001
class PolarDecompose
{
public:
	static double oneNorm(Matrix3d& F);
	static double infNorm(Matrix3d& F);
	static void compute(Matrix3d& F, Matrix3d& R, Matrix3d& S);
};
