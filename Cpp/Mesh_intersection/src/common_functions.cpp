#include "common_functions.h"

libMesh::Real kronecker_delta(unsigned int i,
				   unsigned int j)
{
	return i == j ? 1. : 0.;
};
