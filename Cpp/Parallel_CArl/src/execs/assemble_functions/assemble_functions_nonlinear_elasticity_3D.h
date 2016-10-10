/*
 * assemble_functions_nonlinear_elasticity_3D.h
 *
 *  Created on: Aug 28, 2016
 *      Author: Thiago Milanetto Schlittler
 *
 *  This file is based on libMesh's systems_of_equations_ex7.C
 */

#ifndef EXECS_ASSEMBLE_FUNCTIONS_ASSEMBLE_FUNCTIONS_NONLINEAR_ELASTICITY_3D_H_
#define EXECS_ASSEMBLE_FUNCTIONS_ASSEMBLE_FUNCTIONS_NONLINEAR_ELASTICITY_3D_H_

#include "carl_headers.h"

#include "weight_parameter_function.h"

libMesh::NonlinearImplicitSystem& add_nonlinear_elasticity(	libMesh::EquationSystems& input_systems,
														libMesh::Order order = libMesh::FIRST,
														libMesh::FEFamily family = libMesh::LAGRANGE);

class LargeDeformationElasticity : public libMesh::NonlinearImplicitSystem::ComputeResidual,
                                   public libMesh::NonlinearImplicitSystem::ComputeJacobian
{
private:
	libMesh::EquationSystems &es;
	carl::weight_parameter_function& weight_mask;
	libMesh::SparseMatrix<libMesh::Number> * m_Coupling_Term_Matrix;
	libMesh::NumericVector<libMesh::Number> * m_Coupling_Term_Vector;

public:

	LargeDeformationElasticity (libMesh::EquationSystems &es_in, carl::weight_parameter_function& weight_mask_in) :
		es(es_in),
		weight_mask(weight_mask_in),
		m_Coupling_Term_Matrix(NULL),
		m_Coupling_Term_Vector(NULL)
	{}

	libMesh::Real eval_lambda_1(libMesh::Real E, libMesh::Real mu)
	{
	return mu*(E - 2*mu)/(3*mu-E);
	}

	/**
	* Evaluate the fourth order tensor (C_ijkl) that relates stress to strain.
	*/
	libMesh::Real eval_elasticity_tensor(unsigned int i,
			  unsigned int j,
			  unsigned int k,
			  unsigned int l,
			  libMesh::Number E,
			  libMesh::Number mu)
	{
	const libMesh::Real lambda_1 = eval_lambda_1(E,mu);
	const libMesh::Real lambda_2 = mu;

	return lambda_1 * kronecker_delta(i,j) * kronecker_delta(k,l)
		+ lambda_2 * (kronecker_delta(i,k) * kronecker_delta(j,l)
		+ kronecker_delta(i,l) * kronecker_delta(j,k));
	}

	void set_coupling_term_matrix(libMesh::SparseMatrix<libMesh::Number> * input_matrix)
	{
		m_Coupling_Term_Matrix = input_matrix;
	}

	void set_coupling_term_vector(libMesh::NumericVector<libMesh::Number> * input_vector)
	{
		m_Coupling_Term_Vector = input_vector;
	}
	/**
	* Evaluate the Jacobian of the nonlinear system.
	*/
	virtual void jacobian (	const libMesh::NumericVector<libMesh::Number>& soln,
							libMesh::SparseMatrix<libMesh::Number>&  jacobian,
							libMesh::NonlinearImplicitSystem& /*sys*/);

	/**
	* Evaluate the residual of the nonlinear system.
	*/
	virtual void residual (const libMesh::NumericVector<libMesh::Number>& soln,
							libMesh::NumericVector<libMesh::Number>& residual,
							libMesh::NonlinearImplicitSystem& /*sys*/);

	/**
	* Compute the Cauchy stress for the current solution.
	*/
	void compute_stresses();
};



#endif /* EXECS_ASSEMBLE_FUNCTIONS_ASSEMBLE_FUNCTIONS_NONLINEAR_ELASTICITY_3D_H_ */
