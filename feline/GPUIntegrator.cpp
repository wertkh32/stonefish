#include "GPUIntegrator.h"


GPUIntegrator::GPUIntegrator(Mesh* _mesh, ConstrainedRows* r)
{
	mesh = _mesh;
	n = mesh->getNoOfNodes();
	initVars();
}

void
GPUIntegrator::initVars()
{

}

GPUIntegrator::~GPUIntegrator(void)
{
}
