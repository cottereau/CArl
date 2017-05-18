/*
 * non_linear_FEMSystem.h
 *
 *  Created on: Jul 21, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef EXECS_ASSEMBLE_FUNCTIONS_NON_LINEAR_FEMSYSTEM_H_
#define EXECS_ASSEMBLE_FUNCTIONS_NON_LINEAR_FEMSYSTEM_H_

#include "carl_headers.h"

#include "weight_parameter_function.h"

#include "neohooke_elasticity.h"

// Code adapted from libMesh's FEM example 2

namespace carl
{

class NonLinearSystem: public libMesh::FEMSystem {
public:
	// Constructor
	NonLinearSystem(libMesh::EquationSystems& es, const std::string& name,
			  	  	  const unsigned int number);

	// System initialization
	virtual void init_data();

	// Context initialization
	virtual void init_context(libMesh::DiffContext &context);

	// Element residual and jacobian calculations
	virtual bool element_time_derivative(bool request_jacobian,
			libMesh::DiffContext& context);

	// Contributions for adding boundary conditions
	// virtual bool side_time_derivative(bool request_jacobian,
	// 		libMesh::DiffContext& context);

	virtual bool eulerian_residual(bool, libMesh::DiffContext &) {
	return false;
	}

	// Simulation parameters
	GetPot args;

	// Custom Identifier
	virtual std::string system_type() const {
		return "NonLinearSolid";
	}

	// override method to update mesh also
	virtual void update();

	// save the undeformed mesh to an auxiliary system
	void save_initial_mesh();

	// variable numbers of primary variables in the current system
	unsigned int var[3];

	// variable numbers of primary variables in auxiliary system (for accessing the
	// undeformed configuration)
	unsigned int undefo_var[3];
};

}


#endif /* EXECS_ASSEMBLE_FUNCTIONS_NON_LINEAR_FEMSYSTEM_H_ */
