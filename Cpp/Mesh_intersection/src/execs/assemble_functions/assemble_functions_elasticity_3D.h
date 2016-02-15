/*
 * assemble_functions_elasticity_3D.h
 *
 *  Created on: Nov 30, 2015
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef ELASTICITY_3D_ASSEMBLE_FUNCTIONS_ELASTICITY_3D_H_
#define ELASTICITY_3D_ASSEMBLE_FUNCTIONS_ELASTICITY_3D_H_

#include "common_header_libmesh.h"
#include "common_functions.h"

#include "weight_parameter_function.h"

const bool MASTER_bPerfLog_assemble_fem = false;

// Some classes and functions dealing with the boundary conditions
// ---> Border displacememnt
class border_displacement : public libMesh::FunctionBase<libMesh::Number>
{
	private:
		const unsigned int _u_var, _v_var, _w_var;
		const libMesh::Real _x_displ, _y_displ, _z_displ;

	public:
		border_displacement (	unsigned int u_var,
								unsigned int v_var,
								unsigned int w_var,
								libMesh::Real x_displ = 0,
								libMesh::Real y_displ = 0,
								libMesh::Real z_displ = 0)
								: _u_var(u_var),
								  _v_var(v_var),
								  _w_var(w_var),
								  _x_displ(x_displ),
								  _y_displ(y_displ),
								  _z_displ(z_displ)
		{
			this->_initialized = true;
		}

		virtual libMesh::Number operator() (const libMesh::Point&, const libMesh::Real = 0)
		{
			libmesh_not_implemented();
		}

		virtual void operator() (	const libMesh::Point& p,
									const libMesh::Real,
									libMesh::DenseVector<libMesh::Number>& output)
		{
			output.resize(3);
			output.zero();

			output(_u_var) = _x_displ;
			output(_v_var) = _y_displ;
			output(_w_var) = _z_displ;
		}

		virtual libMesh::UniquePtr<FunctionBase<libMesh::Number> > clone() const
		{
			return libMesh::UniquePtr<FunctionBase<libMesh::Number> >
				(new border_displacement(	_u_var, _v_var, _w_var,
											_x_displ, _y_displ, _z_displ));
		}
};

class boundary_id_cube
{
public:
	int MIN_Z;
	int MIN_Y;
	int MAX_X;
	int MAX_Y;
	int MIN_X;
	int MAX_Z;

	boundary_id_cube()
	{
		MIN_Z = 1;
		MIN_Y = 2;
		MAX_X = 3;
		MAX_Y = 4;
		MIN_X = 5;
		MAX_Z = 6;
	}
};

class boundary_displacement
{
public:
	double x_displ;
	double y_displ;
	double z_displ;

	boundary_displacement(): x_displ { 0.5 }, y_displ { 0 }, z_displ { 0 }
	{
	}

	boundary_displacement(double inputX, double inputY, double inputZ): x_displ { inputX }, y_displ { inputY }, z_displ { inputZ }
	{
	}
};

void set_x_displacement(libMesh::ImplicitSystem& elasticity_system, boundary_displacement& displ, boundary_id_cube& boundary_ids);

void set_displaced_border_translation(libMesh::ImplicitSystem& elasticity_system, boundary_displacement& displ, int boundary_id);

void set_clamped_border(libMesh::ImplicitSystem& elasticity_system, int boundary_id);

libMesh::ExplicitSystem& add_stress(libMesh::EquationSystems& input_systems);

libMesh::LinearImplicitSystem& add_elasticity(	libMesh::EquationSystems& input_systems,
												libMesh::Order order = libMesh::FIRST,
												libMesh::FEFamily family = libMesh::LAGRANGE);

libMesh::LinearImplicitSystem& add_elasticity_with_assemble(	libMesh::EquationSystems& input_systems,
												void fptr(	libMesh::EquationSystems& es,
															const std::string& name),
												libMesh::Order order = libMesh::FIRST,
												libMesh::FEFamily family = libMesh::LAGRANGE);

void Update_SubK(	libMesh::DenseSubMatrix<libMesh::Number>& SubK,
					unsigned int qp,
					unsigned int C_i,
					unsigned int C_k,
					const std::vector<std::vector<libMesh::RealGradient> >& dphi,

					const unsigned int n_components,
					const unsigned int n_u_dofs,
					const std::vector<libMesh::Real>& JxW,
					libMesh::Number E = 1.0,
					libMesh::Number mu = 0.4,
					double cte = 1
					);

void set_physical_properties(libMesh::EquationSystems& es, std::string& physicalParamsFile, double& meanE, double& meanMu);

void set_constant_physical_properties(libMesh::EquationSystems& es, double meanE, double meanMu);

void assemble_elasticity(libMesh::EquationSystems& es,
					   const std::string& system_name);

void assemble_elasticity_heterogeneous(libMesh::EquationSystems& es,
					   const std::string& system_name);

void assemble_elasticity_with_weight(libMesh::EquationSystems& es,
					   const std::string& system_name, carl::weight_parameter_function& weight_mask);

void assemble_elasticity_heterogeneous_with_weight(libMesh::EquationSystems& es,
					   const std::string& system_name, carl::weight_parameter_function& weight_mask);


void compute_stresses(libMesh::EquationSystems& es);

libMesh::Real eval_lambda_1(libMesh::Real E, libMesh::Real mu);

libMesh::Real eval_elasticity_tensor(unsigned int i,
						  unsigned int j,
						  unsigned int k,
						  unsigned int l,
						  libMesh::Number E = 1.,
						  libMesh::Number mu = 0.4);

#endif /* ELASTICITY_3D_ASSEMBLE_FUNCTIONS_ELASTICITY_3D_H_ */
