/*
 * neohooke_elasticity.h
 *
 *  Created on: Jul 21, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef EXECS_ASSEMBLE_FUNCTIONS_NEOHOOKE_ELASTICITY_H_
#define EXECS_ASSEMBLE_FUNCTIONS_NEOHOOKE_ELASTICITY_H_

// Code adapted from libMesh's FEM example 2

#include "carl_headers.h"

#include "weight_parameter_function.h"

namespace carl
{

class NonlinearNeoHookeCurrentConfig
{
public:
	NonlinearNeoHookeCurrentConfig(const std::vector<std::vector<libMesh::RealGradient> >& dphi_in,
								 GetPot& args,
								 bool calculate_linearized_stiffness_in) :
	calculate_linearized_stiffness(calculate_linearized_stiffness_in),
	dphi(dphi_in)
	{
	E = args("material/neohooke/e_modulus", 10000.0);
	nu = args("material/neohooke/nu", 0.3);
	}

	/**
	* Initialize the class for the given displacement gradient at the
	* specified quadrature point.
	*/
	void init_for_qp(libMesh::VectorValue<libMesh::Gradient> & grad_u, unsigned int qp);

	/**
	* Return the residual vector for the current state.
	*/
	void get_residual(libMesh::DenseVector<libMesh::Real> & residuum, unsigned int & i);

	/**
	* Return the stiffness matrix for the current state.
	*/
	void get_linearized_stiffness(libMesh::DenseMatrix<libMesh::Real> & stiffness,
								unsigned int & i, unsigned int & j);

	/**
	* Flag to indicate if it is necessary to calculate values for stiffness
	* matrix during initialization.
	*/
	bool calculate_linearized_stiffness;
private:
	void build_b_0_mat(int i, libMesh::DenseMatrix<libMesh::Real>& b_l_mat);
	void calculate_stress();
	void calculate_tangent();
	static void tensor_to_voigt(const libMesh::RealTensor &tensor, libMesh::DenseVector<libMesh::Real> &vec);

	unsigned int current_qp;
	const std::vector<std::vector<libMesh::RealGradient> >& dphi;

	libMesh::DenseMatrix<libMesh::Real> C_mat;
	libMesh::Real E;
	libMesh::Real nu;
	libMesh::RealTensor F, S, tau, sigma;
	libMesh::DenseMatrix<libMesh::Real> B_L;
	libMesh::DenseMatrix<libMesh::Real> B_K;
};

}


#endif /* EXECS_ASSEMBLE_FUNCTIONS_NEOHOOKE_ELASTICITY_H_ */
