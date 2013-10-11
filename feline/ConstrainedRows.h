#pragma once
#include "QuickArray.h"
#define MAX_ROWS_CONSTRAINED 10000

enum DOF
{
	X = 0, Y = 1, Z = 2
};

struct ConstrainedRows
{
	QuickArray<int,MAX_ROWS_CONSTRAINED> list;
	ConstrainedRows()
	{
	}
	
	void add(int node, DOF dof)
	{
		list.push(node * 3 + dof);
	}
	
	void add(int node)
	{
		list.push(node * 3);
		list.push(node * 3 + 1);
		list.push(node * 3 + 2);
	}


};