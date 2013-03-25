#pragma once
class Integrator
{
	int n;
	float **globalK, **globalM;
public:
	Integrator(int no_of_nodes);
	~Integrator(void);
};

