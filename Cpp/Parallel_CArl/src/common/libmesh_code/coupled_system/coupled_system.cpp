#include "coupled_system.h"

// Class members
void carl::coupled_system::clear()
{
	// Clean up all systems
	libMesh::EquationSystems *EqBIGSys = m_BIG_EquationSystem.second;
	EqBIGSys->clear();
	delete EqBIGSys;
	EqBIGSys = NULL;
	m_BIG_EquationSystem.second = NULL;

	while(!m_micro_EquationSystemMap.empty())
	{
		EqSystem_iterator toClean = m_micro_EquationSystemMap.begin();

		libMesh::EquationSystems *EqSys = toClean->second;
		EqSys->clear();
		delete EqSys;
		EqSys = NULL;

		m_micro_EquationSystemMap.erase(toClean);
	}

	if(m_bHasDefinedMeshRestrictions)
	{
		libMesh::EquationSystems *EqRBIGSys = m_R_BIG_EquationSystem.second;
		EqRBIGSys->clear();
		delete EqRBIGSys;
		EqRBIGSys = NULL;
		m_R_BIG_EquationSystem.second = NULL;

		while(!m_R_micro_EquationSystemMap.empty())
		{
			EqSystem_iterator toClean = m_R_micro_EquationSystemMap.begin();

			libMesh::EquationSystems *EqSys = toClean->second;
			EqSys->clear();
			delete EqSys;
			EqSys = NULL;

			m_R_micro_EquationSystemMap.erase(toClean);
		}
	}

	if(m_bHasDefinedCoordVector_BIG)
	{
		libMesh::PetscVector<libMesh::Number> *VectBIG = m_coord_vect_BIG.second;
		VectBIG->clear();
		delete VectBIG;
		VectBIG = NULL;
		 m_coord_vect_BIG.second = NULL;
	}

	if(m_bHasDefinedCoordVector_micro)
	{
		while(!m_coord_vect_microMap.empty())
		{
			Vector_iterator toClean = m_coord_vect_microMap.begin();

			libMesh::PetscVector<libMesh::Number> *Vect = toClean->second;
			Vect->clear();
			delete Vect;
			Vect = NULL;

			m_coord_vect_microMap.erase(toClean);
		}
	}
	while(!m_inter_EquationSystemMap.empty())
	{
		EqSystem_iterator toClean = m_inter_EquationSystemMap.begin();

		libMesh::EquationSystems *EqSys = toClean->second;
		EqSys->clear();
		delete EqSys;
		EqSys = NULL;

		m_inter_EquationSystemMap.erase(toClean);
	}

	while(!m_mediator_EquationSystemMap.empty())
	{
		EqSystem_iterator toClean = m_mediator_EquationSystemMap.begin();

		libMesh::EquationSystems *EqSys = toClean->second;
		EqSys->clear();
		delete EqSys;
		EqSys = NULL;

		m_mediator_EquationSystemMap.erase(toClean);
	}

	while(!m_couplingMatrixMap_mediator_micro.empty())
	{
		Matrix_iterator toClean = m_couplingMatrixMap_mediator_micro.begin();

		libMesh::PetscMatrix<libMesh::Number> *Mat = toClean->second;
		Mat->clear();
		delete Mat;
		Mat = NULL;

		m_couplingMatrixMap_mediator_micro.erase(toClean);
	}

	while(!m_couplingMatrixMap_mediator_BIG.empty())
	{
		Matrix_iterator toClean = m_couplingMatrixMap_mediator_BIG.begin();

		libMesh::PetscMatrix<libMesh::Number> *Mat = toClean->second;
		Mat->clear();
		delete Mat;
		Mat = NULL;

		m_couplingMatrixMap_mediator_BIG.erase(toClean);
	}

	while(!m_couplingMatrixMap_mediator_mediator.empty())
	{
		Matrix_iterator toClean = m_couplingMatrixMap_mediator_mediator.begin();

		libMesh::PetscMatrix<libMesh::Number> *Mat = toClean->second;
		Mat->clear();
		delete Mat;
		Mat = NULL;

		m_couplingMatrixMap_mediator_mediator.erase(toClean);
	}

	while(!m_alpha_masks.empty())
	{
		alpha_mask_iterator toClean = m_alpha_masks.begin();

		weight_parameter_function *alpha = toClean->second;
		alpha->clear();
		delete alpha;

		m_alpha_masks.erase(toClean);
	}
};

void carl::coupled_system::set_corrected_shapes(	const std::vector<std::vector<libMesh::Real> >& 	lambda_weights,
							const std::vector<std::vector<libMesh::Real> >& 	phi_inter,
							std::vector<std::vector<libMesh::Real> >& 			phi_corrected)
{
	unsigned int lambda_corr_size  	= lambda_weights.size();
	unsigned int lambda_inter_size 	= lambda_weights[0].size();

	unsigned int phi_inter_dof		= phi_inter.size();
	unsigned int phi_inter_qp		= phi_inter[0].size();

	unsigned int phi_corrected_dof		= phi_corrected.size();
	unsigned int phi_corrected_qp		= phi_corrected[0].size();

	/*
	 * 		lambda 	: [ n_dofs_corr ] x [ n_dofs_inter ]
	 * 		inter  	: [ n_dofs_inter ] x [ n_qp ]
	 * 		corr	: [ n_dofs_corr ] x [ n_qp ]
	 */

	homemade_assert_msg(phi_inter_qp == phi_corrected_qp,
					" Different numbers of quadrature points!");
	homemade_assert_msg(lambda_corr_size == phi_corrected_dof,
					" Incompatible corrected shape table and barycentric coordinates!");
	homemade_assert_msg(lambda_inter_size == phi_inter_dof,
					" Incompatible intersection shape table and barycentric coordinates!");

	// phi_A,i (qp) = l_i,1 * phi_I,1 (qp) + l_i,2 * phi_I,2 (qp) ...
	//			    = sum [ l_i,j * phi_I,j (qp) ]
	for(unsigned int qp = 0; qp < phi_inter_qp; ++qp)
	{
		for(unsigned int iii = 0; iii < phi_corrected_dof; ++iii)
		{
			phi_corrected[iii][qp] = 0;
			for(unsigned int jjj = 0; jjj < phi_inter_dof; ++jjj)
			{
				phi_corrected[iii][qp] += lambda_weights[iii][jjj]*phi_inter[jjj][qp];
			}
		}
	}
}

void carl::coupled_system::set_corrected_shape_gradients(	const std::vector<std::vector<libMesh::Real> >& 	lambda_weights,
		const std::vector<std::vector<libMesh::RealGradient> >& 	dphi_inter,
		std::vector<std::vector<libMesh::RealGradient> >& 			dphi_corrected)
{
	unsigned int lambda_corr_size  	= lambda_weights.size();
	unsigned int lambda_inter_size 	= lambda_weights[0].size();

	unsigned int phi_inter_dof		= dphi_inter.size();
	unsigned int phi_inter_qp		= dphi_inter[0].size();

	unsigned int phi_corrected_dof		= dphi_corrected.size();
	unsigned int phi_corrected_qp		= dphi_corrected[0].size();

	/*
	 * 		lambda 	: [ n_dofs_corr ] x [ n_dofs_inter ]
	 * 		inter  	: [ n_dofs_inter ] x [ n_qp ]
	 * 		corr	: [ n_dofs_corr ] x [ n_qp ]
	 */

	homemade_assert_msg(phi_inter_qp == phi_corrected_qp,
					" Different numbers of quadrature points!");
	homemade_assert_msg(lambda_corr_size == phi_corrected_dof,
					" Incompatible corrected shape table and barycentric coordinates!");
	homemade_assert_msg(lambda_inter_size == phi_inter_dof,
					" Incompatible intersection shape table and barycentric coordinates!");

	// phi_A,i (qp) = l_i,1 * phi_I,1 (qp) + l_i,2 * phi_I,2 (qp) ...
	//			    = sum [ l_i,j * phi_I,j (qp) ]
	for(unsigned int qp = 0; qp < phi_inter_qp; ++qp)
	{
		for(unsigned int iii = 0; iii < phi_corrected_dof; ++iii)
		{
			dphi_corrected[iii][qp] = 0;
			for(unsigned int jjj = 0; jjj < phi_inter_dof; ++jjj)
			{
				dphi_corrected[iii][qp] += lambda_weights[iii][jjj]*dphi_inter[jjj][qp];
			}
		}
	}
};

void carl::coupled_system::get_lambdas(	const unsigned int 							dim,
					const libMesh::FEType& 						fe_t,
					const libMesh::Elem* 						base_elem,
					const std::vector<libMesh::Point>& 			phys_points,
					std::vector<libMesh::Point>& 				ref_points,
					std::vector<std::vector<libMesh::Real> >& 	lambda_weights)
{
	// Test if we are using the correct type!
	homemade_assert_msg(base_elem->type() == libMesh::TET4 || base_elem->type() == libMesh::HEX8,
			" Only implemented for TET4 or HEX8 elements!");

	homemade_assert_msg(dim == 3,
			" Only implemented for 3D!");

	unsigned int lambda_inter_size = lambda_weights[0].size();
	unsigned int lambda_base_size  = lambda_weights.size();

	// Set vectors to have same sizes
	if(phys_points.size() != ref_points.size())
	{
		ref_points.resize( phys_points.size() );
	}

	// Convert the points
	libMesh::FEInterface::inverse_map(dim, fe_t, base_elem, phys_points, ref_points);

//	for(unsigned int iii = 0; iii < ref_points.size(); ++iii)
//	{
//		test_barycentric = 	ref_points[iii](0) < 1 + 1E-16 &&
//							ref_points[iii](0) > - 1E-16 &&
//							ref_points[iii](1) < 1 + 1E-16 &&
//							ref_points[iii](1) > - 1E-16 &&
//							ref_points[iii](2) < 1 + 1E-16 &&
//							ref_points[iii](2) > - 1E-16;
//
//			homemade_assert_msg(test_barycentric,
//						" One of the barycentric coordinates is wrong!");
//	}

	// Calculate the lambdas
	// -> lines : DoF from base
	// -> cols  : DoF from inter
	for(unsigned int iii = 0; iii < lambda_base_size; ++iii)
	{
		for(unsigned int jjj = 0; jjj < lambda_inter_size; ++jjj)
		{
			lambda_weights[iii][jjj] =
					libMesh::FE<3,libMesh::LAGRANGE>::shape(
							base_elem->type(),
							libMesh::FIRST,
							iii,
							ref_points[jjj]);
		}
	}
};

libMesh::PetscMatrix<libMesh::Number>& carl::coupled_system::get_micro_coupling_matrix(const std::string& name)
{
	return * m_couplingMatrixMap_mediator_micro[name];
}

libMesh::PetscMatrix<libMesh::Number>& carl::coupled_system::get_BIG_coupling_matrix(const std::string& name)
{
	return * m_couplingMatrixMap_mediator_BIG[name];
}

libMesh::PetscMatrix<libMesh::Number>& carl::coupled_system::get_mediator_coupling_matrix(const std::string& name)
{
	return * m_couplingMatrixMap_mediator_mediator[name];
}

void carl::coupled_system::print_matrix_micro_info(const std::string& name)
{
	libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix =
						* m_couplingMatrixMap_mediator_micro[name];
	std::cout << "| Restrict - Micro matrix -> " << name << std::endl;
	print_matrix(CouplingTestMatrix);
}

void carl::coupled_system::set_restart(		bool bUseRestart,
											bool bPrintRestart,
											const std::string restart_base_filename,
											bool bPrintMatrix)
{
	// Set the solver parameters
	switch(m_solver_type)
	{
	case LATIN_MODIFIED_STIFFNESS:
	case LATIN_ORIGINAL_STIFFNESS:
		{
			std::shared_ptr<PETSC_LATIN_coupled_solver> cast_LATIN_solver = std::dynamic_pointer_cast<PETSC_LATIN_coupled_solver>(m_coupled_solver);

			cast_LATIN_solver->set_restart(bUseRestart,bPrintRestart,restart_base_filename);
			cast_LATIN_solver->set_info(bPrintMatrix,restart_base_filename);
		}
		break;
	case CG:
		{
			std::shared_ptr<PETSC_CG_coupled_solver> cast_CG_solver = std::dynamic_pointer_cast<PETSC_CG_coupled_solver>(m_coupled_solver);

			cast_CG_solver->set_restart(bUseRestart,bPrintRestart,restart_base_filename);
			cast_CG_solver->set_info(bPrintMatrix,restart_base_filename);
		}
		break;
	}
};

void carl::coupled_system::set_macro_system(
		const std::string micro_name,
		const std::string type_name,
		void fptr_assemble(		libMesh::EquationSystems& es,
								const std::string& name, weight_parameter_function& weight_mask))
{
	// Get the system
	libMesh::EquationSystems& EqSystems_BIG = * m_BIG_EquationSystem.second;
	libMesh::ImplicitSystem& Sys_BIG =   libMesh::cast_ref<libMesh::ImplicitSystem&>(EqSystems_BIG.get_system(type_name));

	// Assemble the systems
	if(fptr_assemble)
	{
		// Get the weight functions
		weight_parameter_function& alpha_mask = * m_alpha_masks[micro_name];

		// Assemble
		fptr_assemble(EqSystems_BIG,type_name,alpha_mask);
	}
	else
	{
		Sys_BIG.assemble();
	}

	// Close matrix and vector
	Sys_BIG.matrix->close();
	Sys_BIG.rhs->close();
	m_bHasAssembled_BIG = true;
}

void carl::coupled_system::set_micro_system(
		const std::string micro_name,
		const std::string type_name,
		void fptr_assemble(		libMesh::EquationSystems& es,
										const std::string& name, weight_parameter_function& weight_mask))
{
	// Get the system
	libMesh::EquationSystems& EqSystems_micro = * m_micro_EquationSystemMap[micro_name];
	libMesh::ImplicitSystem& Sys_micro = libMesh::cast_ref<libMesh::ImplicitSystem&>(EqSystems_micro.get_system(type_name));

	// Assemble the systems
	if(fptr_assemble)
	{
		// Get the weight functions
		weight_parameter_function& alpha_mask = * m_alpha_masks[micro_name];

		// Assemble
		fptr_assemble(EqSystems_micro,type_name,alpha_mask);
	}
	else
	{
		Sys_micro.assemble();
	}

	// Close matrix and vector
	Sys_micro.matrix->close();
	Sys_micro.rhs->close();
	m_bHasAssembled_micro[micro_name] = true;
}

void carl::coupled_system::set_micro_system(
		const std::string micro_name,
		const std::string type_name,
		void fptr_assemble(		libMesh::EquationSystems& es,
								const std::string& name, weight_parameter_function& weight_mask,
								anisotropic_elasticity_tensor_cubic_sym& anisotropy_obj_input),
								anisotropic_elasticity_tensor_cubic_sym& anisotropy_obj)
{
	// Get the system
	libMesh::EquationSystems& EqSystems_micro = * m_micro_EquationSystemMap[micro_name];
	libMesh::ImplicitSystem& Sys_micro = libMesh::cast_ref<libMesh::ImplicitSystem&>(EqSystems_micro.get_system(type_name));

	// Assemble the systems
	if(fptr_assemble)
	{
		// Get the weight functions
		weight_parameter_function& alpha_mask = * m_alpha_masks[micro_name];

		// Assemble
		fptr_assemble(EqSystems_micro,type_name,alpha_mask,anisotropy_obj);
	}
	else
	{
		Sys_micro.assemble();
	}

	// Close matrix and vector
	Sys_micro.matrix->close();
	Sys_micro.rhs->close();
	m_bHasAssembled_micro[micro_name] = true;
}

void carl::coupled_system::set_LATIN_solver(const std::string micro_name,
		const std::string type_name,
		double k_dA, double k_dB, double k_cA, double k_cB,
		double eps, int convIter, double relax)
{
	this->set_macro_system(micro_name,type_name);
	this->set_micro_system(micro_name,type_name);

	libMesh::PetscMatrix<libMesh::Number>& C_RA = * m_couplingMatrixMap_mediator_BIG[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& C_RB = * m_couplingMatrixMap_mediator_micro[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& C_RR = * m_couplingMatrixMap_mediator_mediator[micro_name];

	libMesh::PetscMatrix<libMesh::Number>& M_A = libMesh::cast_ref<libMesh::PetscMatrix<libMesh::Number>& >(
			* libMesh::cast_ref<libMesh::ImplicitSystem&>(m_BIG_EquationSystem.second->get_system(type_name)).matrix);
	libMesh::PetscVector<libMesh::Number>& F_A =libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>& >(
			* libMesh::cast_ref<libMesh::ImplicitSystem&>(m_BIG_EquationSystem.second->get_system(type_name)).rhs);

	libMesh::PetscMatrix<libMesh::Number>& M_B = libMesh::cast_ref<libMesh::PetscMatrix<libMesh::Number>& >(
			* libMesh::cast_ref<libMesh::ImplicitSystem&>(m_micro_EquationSystemMap[micro_name]->get_system(type_name)).matrix);
	libMesh::PetscVector<libMesh::Number>& F_B =libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>& >(
			* libMesh::cast_ref<libMesh::ImplicitSystem&>(m_micro_EquationSystemMap[micro_name]->get_system(type_name)).rhs);

	if(m_bUseNullSpace_BIG)
	{
		std::cout << "| -> Using null space for macro system!!!" << std::endl;
		this->set_rigid_body_modes_BIG(type_name);
	}
	if(m_bUseNullSpace_micro[micro_name])
	{
		std::cout << "| -> Using null space for micro system " << micro_name << "!!!" << std::endl;
		this->set_rigid_body_modes_micro(micro_name,type_name);
	}

	// Set the solver parameters
	switch(m_solver_type)
	{
	case LATIN_MODIFIED_STIFFNESS:
	case LATIN_ORIGINAL_STIFFNESS:
		{
			std::shared_ptr<PETSC_LATIN_coupled_solver> cast_LATIN_solver = std::dynamic_pointer_cast<PETSC_LATIN_coupled_solver>(m_coupled_solver);

			cast_LATIN_solver->set_params(k_dA,k_dB,k_cA,k_cB);

			// Set the solver matrices
			cast_LATIN_solver->set_matrices(M_A,M_B,C_RA,C_RB,C_RR);

			// Set the solver matrices
			cast_LATIN_solver->set_forces(F_A,F_B);

			// Set LATIN parameters (convergence, relaxation ... )
			cast_LATIN_solver->set_convergence_limits(eps,convIter);
			cast_LATIN_solver->set_relaxation(relax);
			break;
		}
	case CG:
		break;
	}
};

void carl::coupled_system::set_LATIN_solver(const std::string micro_name, const std::string type_name,
												void fptr_BIG(		libMesh::EquationSystems& es,
																	const std::string& name, weight_parameter_function& weight_mask),
												void fptr_micro(	libMesh::EquationSystems& es,
																	const std::string& name, weight_parameter_function& weight_mask),
												double k_dA, double k_dB, double k_cA, double k_cB,
												double eps, int convIter, double relax)
{
	this->set_macro_system(micro_name,type_name,fptr_BIG);
	this->set_micro_system(micro_name,type_name,fptr_micro);

	// Get the matrices
	libMesh::PetscMatrix<libMesh::Number>& C_RA = * m_couplingMatrixMap_mediator_BIG[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& C_RB = * m_couplingMatrixMap_mediator_micro[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& C_RR = * m_couplingMatrixMap_mediator_mediator[micro_name];

	libMesh::PetscMatrix<libMesh::Number>& M_A = libMesh::cast_ref<libMesh::PetscMatrix<libMesh::Number>& >(
			* libMesh::cast_ref<libMesh::ImplicitSystem&>(m_BIG_EquationSystem.second->get_system(type_name)).matrix);
	libMesh::PetscVector<libMesh::Number>& F_A =libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>& >(
			* libMesh::cast_ref<libMesh::ImplicitSystem&>(m_BIG_EquationSystem.second->get_system(type_name)).rhs);

	libMesh::PetscMatrix<libMesh::Number>& M_B = libMesh::cast_ref<libMesh::PetscMatrix<libMesh::Number>& >(
			* libMesh::cast_ref<libMesh::ImplicitSystem&>(m_micro_EquationSystemMap[micro_name]->get_system(type_name)).matrix);
	libMesh::PetscVector<libMesh::Number>& F_B =libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>& >(
			* libMesh::cast_ref<libMesh::ImplicitSystem&>(m_micro_EquationSystemMap[micro_name]->get_system(type_name)).rhs);

	if(m_bUseNullSpace_BIG)
	{
		std::cout << "| -> Using null space for macro system!!!" << std::endl;
		this->set_rigid_body_modes_BIG(type_name);
	}
	if(m_bUseNullSpace_micro[micro_name])
	{
		std::cout << "| -> Using null space for micro system " << micro_name << "!!!" << std::endl;
		this->set_rigid_body_modes_micro(micro_name,type_name);
	}

	switch(m_solver_type)
	{
	case LATIN_MODIFIED_STIFFNESS:
	case LATIN_ORIGINAL_STIFFNESS:
		{
			std::shared_ptr<PETSC_LATIN_coupled_solver> cast_LATIN_solver = std::dynamic_pointer_cast<PETSC_LATIN_coupled_solver>(m_coupled_solver);

			// Set the solver parameters
			cast_LATIN_solver->set_params(k_dA,k_dB,k_cA,k_cB);

			// Set the solver matrices
			cast_LATIN_solver->set_matrices(M_A,M_B,C_RA,C_RB,C_RR);

			// Set the solver matrices
			cast_LATIN_solver->set_forces(F_A,F_B);

			// Set LATIN parameters (convergence, relaxation ... )
			cast_LATIN_solver->set_convergence_limits(eps,convIter);
			cast_LATIN_solver->set_relaxation(relax);
			break;
		}
	case CG:
		break;
	}
};

void carl::coupled_system::set_LATIN_solver(const std::string micro_name, const std::string type_name,
												void fptr_BIG(		libMesh::EquationSystems& es,
																	const std::string& name, weight_parameter_function& weight_mask),
												void fptr_micro(	libMesh::EquationSystems& es,
																	const std::string& name, weight_parameter_function& weight_mask, anisotropic_elasticity_tensor_cubic_sym& anisotropy_obj_input),
												anisotropic_elasticity_tensor_cubic_sym& anisotropy_obj,
												double k_dA, double k_dB, double k_cA, double k_cB,
												double eps, int convIter, double relax)
{
	this->set_macro_system(micro_name,type_name,fptr_BIG);
	this->set_micro_system(micro_name,type_name,fptr_micro,anisotropy_obj);

	// Get the matrices
	libMesh::PetscMatrix<libMesh::Number>& C_RA = * m_couplingMatrixMap_mediator_BIG[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& C_RB = * m_couplingMatrixMap_mediator_micro[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& C_RR = * m_couplingMatrixMap_mediator_mediator[micro_name];

	libMesh::PetscMatrix<libMesh::Number>& M_A = libMesh::cast_ref<libMesh::PetscMatrix<libMesh::Number>& >(
			* libMesh::cast_ref<libMesh::ImplicitSystem&>(m_BIG_EquationSystem.second->get_system(type_name)).matrix);
	libMesh::PetscVector<libMesh::Number>& F_A =libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>& >(
			* libMesh::cast_ref<libMesh::ImplicitSystem&>(m_BIG_EquationSystem.second->get_system(type_name)).rhs);

	libMesh::PetscMatrix<libMesh::Number>& M_B = libMesh::cast_ref<libMesh::PetscMatrix<libMesh::Number>& >(
			* libMesh::cast_ref<libMesh::ImplicitSystem&>(m_micro_EquationSystemMap[micro_name]->get_system(type_name)).matrix);
	libMesh::PetscVector<libMesh::Number>& F_B =libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>& >(
			* libMesh::cast_ref<libMesh::ImplicitSystem&>(m_micro_EquationSystemMap[micro_name]->get_system(type_name)).rhs);

	if(m_bUseNullSpace_BIG)
	{
		std::cout << "| -> Using null space for macro system!!!" << std::endl;
		this->set_rigid_body_modes_BIG(type_name);
	}
	if(m_bUseNullSpace_micro[micro_name])
	{
		std::cout << "| -> Using null space for micro system " << micro_name << "!!!" << std::endl;
		this->set_rigid_body_modes_micro(micro_name,type_name);
	}

	switch(m_solver_type)
	{
	case LATIN_MODIFIED_STIFFNESS:
	case LATIN_ORIGINAL_STIFFNESS:
		{
			std::shared_ptr<PETSC_LATIN_coupled_solver> cast_LATIN_solver = std::dynamic_pointer_cast<PETSC_LATIN_coupled_solver>(m_coupled_solver);

			// Set the solver parameters
			cast_LATIN_solver->set_params(k_dA,k_dB,k_cA,k_cB);

			// Set the solver matrices
			cast_LATIN_solver->set_matrices(M_A,M_B,C_RA,C_RB,C_RR);

			// Set the solver matrices
			cast_LATIN_solver->set_forces(F_A,F_B);

			// Set LATIN parameters (convergence, relaxation ... )
			cast_LATIN_solver->set_convergence_limits(eps,convIter);
			cast_LATIN_solver->set_relaxation(relax);
			break;
		}
	case CG:
		break;
	}
};

void carl::coupled_system::set_CG_solver(const std::string micro_name,
		const std::string type_name,
		double eps_abs, double eps_rel, int convIter, double div_tol)
{
	this->set_macro_system(micro_name,type_name);
	this->set_micro_system(micro_name,type_name);

	// Get the matrices
	libMesh::PetscMatrix<libMesh::Number>& C_RA = * m_couplingMatrixMap_mediator_BIG[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& C_RB = * m_couplingMatrixMap_mediator_micro[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& C_RR = * m_couplingMatrixMap_mediator_mediator[micro_name];

	libMesh::PetscMatrix<libMesh::Number>& M_A = libMesh::cast_ref<libMesh::PetscMatrix<libMesh::Number>& >(
			* libMesh::cast_ref<libMesh::ImplicitSystem&>(m_BIG_EquationSystem.second->get_system(type_name)).matrix);
	libMesh::PetscVector<libMesh::Number>& F_A =libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>& >(
			* libMesh::cast_ref<libMesh::ImplicitSystem&>(m_BIG_EquationSystem.second->get_system(type_name)).rhs);

	libMesh::PetscMatrix<libMesh::Number>& M_B = libMesh::cast_ref<libMesh::PetscMatrix<libMesh::Number>& >(
			* libMesh::cast_ref<libMesh::ImplicitSystem&>(m_micro_EquationSystemMap[micro_name]->get_system(type_name)).matrix);
	libMesh::PetscVector<libMesh::Number>& F_B =libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>& >(
			* libMesh::cast_ref<libMesh::ImplicitSystem&>(m_micro_EquationSystemMap[micro_name]->get_system(type_name)).rhs);

	if(m_bUseNullSpace_BIG)
	{
		std::cout << "| -> Using null space for macro system!!!" << std::endl;
		this->set_rigid_body_modes_BIG(type_name);
	}
	if(m_bUseNullSpace_micro[micro_name])
	{
		std::cout << "| -> Using null space for micro system " << micro_name << "!!!" << std::endl;
		this->set_rigid_body_modes_micro(micro_name,type_name);
	}

	// Set the solver parameters
	switch(m_solver_type)
	{
	case CG:
		{
			std::shared_ptr<PETSC_CG_coupled_solver> cast_CG_solver = std::dynamic_pointer_cast<PETSC_CG_coupled_solver>(m_coupled_solver);
			std::shared_ptr<generic_solver_interface> cast_sys_A_solver = std::dynamic_pointer_cast<KSP_linear_solver>(m_sys_A_solver);
			std::shared_ptr<generic_solver_interface> cast_sys_B_solver = std::dynamic_pointer_cast<KSP_linear_solver>(m_sys_B_solver);

			// Use preconditioner
			cast_CG_solver->set_preconditioner_type(m_CG_precond_type);

			// Set the solver matrices
			cast_CG_solver->set_matrices(M_A,M_B,C_RA,C_RB,C_RR);

			// Set the solver matrices
			cast_CG_solver->set_forces(F_A,F_B);

			if(m_bUseNullSpace_micro[micro_name])
			{
				cast_CG_solver->build_null_space_projection_matrices(M_B,C_RB);
			}

			// Set CG parameters (convergence )
			cast_CG_solver->set_convergence_limits(eps_abs,eps_rel,convIter,div_tol);

			// Set the system solvers
			cast_CG_solver->set_solvers(cast_sys_A_solver.get(),cast_sys_B_solver.get());
			break;
		}
	case LATIN_MODIFIED_STIFFNESS:
	case LATIN_ORIGINAL_STIFFNESS:
		break;
	}
};

void carl::coupled_system::set_CG_solver(const std::string micro_name, const std::string type_name,
												void fptr_BIG(		libMesh::EquationSystems& es,
																	const std::string& name, weight_parameter_function& weight_mask),
												void fptr_micro(	libMesh::EquationSystems& es,
																	const std::string& name, weight_parameter_function& weight_mask),
												double eps_abs, double eps_rel, int convIter, double div_tol)
{
	this->set_macro_system(micro_name,type_name,fptr_BIG);
	this->set_micro_system(micro_name,type_name,fptr_micro);

	// Get the matrices
	libMesh::PetscMatrix<libMesh::Number>& C_RA = * m_couplingMatrixMap_mediator_BIG[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& C_RB = * m_couplingMatrixMap_mediator_micro[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& C_RR = * m_couplingMatrixMap_mediator_mediator[micro_name];

	libMesh::PetscMatrix<libMesh::Number>& M_A = libMesh::cast_ref<libMesh::PetscMatrix<libMesh::Number>& >(
			* libMesh::cast_ref<libMesh::ImplicitSystem&>(m_BIG_EquationSystem.second->get_system(type_name)).matrix);
	libMesh::PetscVector<libMesh::Number>& F_A =libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>& >(
			* libMesh::cast_ref<libMesh::ImplicitSystem&>(m_BIG_EquationSystem.second->get_system(type_name)).rhs);

	libMesh::PetscMatrix<libMesh::Number>& M_B = libMesh::cast_ref<libMesh::PetscMatrix<libMesh::Number>& >(
			* libMesh::cast_ref<libMesh::ImplicitSystem&>(m_micro_EquationSystemMap[micro_name]->get_system(type_name)).matrix);
	libMesh::PetscVector<libMesh::Number>& F_B =libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>& >(
			* libMesh::cast_ref<libMesh::ImplicitSystem&>(m_micro_EquationSystemMap[micro_name]->get_system(type_name)).rhs);

	if(m_bUseNullSpace_BIG)
	{
		std::cout << "| -> Using null space for macro system!!!" << std::endl;
		this->set_rigid_body_modes_BIG(type_name);
	}
	if(m_bUseNullSpace_micro[micro_name])
	{
		std::cout << "| -> Using null space for micro system " << micro_name << "!!!" << std::endl;
		this->set_rigid_body_modes_micro(micro_name,type_name);
	}

	// Set the solver parameters
	switch(m_solver_type)
	{
	case CG:
		{
			std::shared_ptr<PETSC_CG_coupled_solver> cast_CG_solver = std::dynamic_pointer_cast<PETSC_CG_coupled_solver>(m_coupled_solver);
			std::shared_ptr<generic_solver_interface> cast_sys_A_solver = std::dynamic_pointer_cast<KSP_linear_solver>(m_sys_A_solver);
			std::shared_ptr<generic_solver_interface> cast_sys_B_solver = std::dynamic_pointer_cast<KSP_linear_solver>(m_sys_B_solver);

			// Use preconditioner
			cast_CG_solver->set_preconditioner_type(m_CG_precond_type);

			// Set the solver matrices
			cast_CG_solver->set_matrices(M_A,M_B,C_RA,C_RB,C_RR);

			// Set the solver matrices
			cast_CG_solver->set_forces(F_A,F_B);

			if(m_bUseNullSpace_micro[micro_name])
			{
				cast_CG_solver->build_null_space_projection_matrices(M_B,C_RB);
			}

			// Set CG parameters (convergence )
			cast_CG_solver->set_convergence_limits(eps_abs,eps_rel,convIter,div_tol);

			// Set the system solvers
			cast_CG_solver->set_solvers(cast_sys_A_solver.get(),cast_sys_B_solver.get());
			break;
		}
	case LATIN_MODIFIED_STIFFNESS:
	case LATIN_ORIGINAL_STIFFNESS:
		break;
	}
};

void carl::coupled_system::set_CG_solver(const std::string micro_name, const std::string type_name,
												void fptr_BIG(		libMesh::EquationSystems& es,
																	const std::string& name, weight_parameter_function& weight_mask),
												void fptr_micro(	libMesh::EquationSystems& es,
																	const std::string& name, weight_parameter_function& weight_mask, anisotropic_elasticity_tensor_cubic_sym& anisotropy_obj_input),
												anisotropic_elasticity_tensor_cubic_sym& anisotropy_obj,
												double eps_abs, double eps_rel, int convIter, double div_tol)
{
	this->set_macro_system(micro_name,type_name,fptr_BIG);
	this->set_micro_system(micro_name,type_name,fptr_micro,anisotropy_obj);

	// Get the matrices
	libMesh::PetscMatrix<libMesh::Number>& C_RA = * m_couplingMatrixMap_mediator_BIG[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& C_RB = * m_couplingMatrixMap_mediator_micro[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& C_RR = * m_couplingMatrixMap_mediator_mediator[micro_name];

	libMesh::PetscMatrix<libMesh::Number>& M_A = libMesh::cast_ref<libMesh::PetscMatrix<libMesh::Number>& >(
			* libMesh::cast_ref<libMesh::ImplicitSystem&>(m_BIG_EquationSystem.second->get_system(type_name)).matrix);
	libMesh::PetscVector<libMesh::Number>& F_A =libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>& >(
			* libMesh::cast_ref<libMesh::ImplicitSystem&>(m_BIG_EquationSystem.second->get_system(type_name)).rhs);

	libMesh::PetscMatrix<libMesh::Number>& M_B = libMesh::cast_ref<libMesh::PetscMatrix<libMesh::Number>& >(
			* libMesh::cast_ref<libMesh::ImplicitSystem&>(m_micro_EquationSystemMap[micro_name]->get_system(type_name)).matrix);
	libMesh::PetscVector<libMesh::Number>& F_B =libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>& >(
			* libMesh::cast_ref<libMesh::ImplicitSystem&>(m_micro_EquationSystemMap[micro_name]->get_system(type_name)).rhs);

	if(m_bUseNullSpace_BIG)
	{
		std::cout << "| -> Using null space for macro system!!!" << std::endl;
		this->set_rigid_body_modes_BIG(type_name);
	}
	if(m_bUseNullSpace_micro[micro_name])
	{
		std::cout << "| -> Using null space for micro system " << micro_name << "!!!" << std::endl;
		this->set_rigid_body_modes_micro(micro_name,type_name);
	}

	// Set the solver parameters
	switch(m_solver_type)
	{
	case CG:
		{
			std::shared_ptr<PETSC_CG_coupled_solver> cast_CG_solver = std::dynamic_pointer_cast<PETSC_CG_coupled_solver>(m_coupled_solver);
			std::shared_ptr<generic_solver_interface> cast_sys_A_solver = std::dynamic_pointer_cast<KSP_linear_solver>(m_sys_A_solver);
			std::shared_ptr<generic_solver_interface> cast_sys_B_solver = std::dynamic_pointer_cast<KSP_linear_solver>(m_sys_B_solver);

			// Use preconditioner
			cast_CG_solver->set_preconditioner_type(m_CG_precond_type);

			// Set the solver matrices
			cast_CG_solver->set_matrices(M_A,M_B,C_RA,C_RB,C_RR);

			// Set the solver matrices
			cast_CG_solver->set_forces(F_A,F_B);

			if(m_bUseNullSpace_micro[micro_name])
			{
				cast_CG_solver->build_null_space_projection_matrices(M_B,C_RB);
			}

			// Set CG parameters (convergence )
			cast_CG_solver->set_convergence_limits(eps_abs,eps_rel,convIter,div_tol);

			// Set the system solvers
			cast_CG_solver->set_solvers(cast_sys_A_solver.get(),cast_sys_B_solver.get());
			break;
		}
	case LATIN_MODIFIED_STIFFNESS:
	case LATIN_ORIGINAL_STIFFNESS:
		break;
	}
};

void carl::coupled_system::set_LATIN_nonlinear_solver(const std::string micro_name, const std::string type_name_BIG,
		 	 	 	 	 	 	 	 	 	 	 const std::string type_name_micro,
												void fptr_BIG(		libMesh::EquationSystems& es,
																	const std::string& name, weight_parameter_function& weight_mask),
												libMesh::EquationSystems& eq_sys_nonlinear,
												double k_dA, double k_dB, double k_cA, double k_cB,
												double eps, int convIter, double relax)
{
	this->set_macro_system(micro_name,type_name_BIG,fptr_BIG);

	// Get the systems
	libMesh::EquationSystems& EqSystems_micro = * m_micro_EquationSystemMap[micro_name];

	// Since the micro system is non-linear, the assemble steps must be done during the solves

	// Get the matrices
	libMesh::PetscMatrix<libMesh::Number>& C_RA = * m_couplingMatrixMap_mediator_BIG[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& C_RB = * m_couplingMatrixMap_mediator_micro[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& C_RR = * m_couplingMatrixMap_mediator_mediator[micro_name];

	libMesh::PetscMatrix<libMesh::Number>& M_A = libMesh::cast_ref<libMesh::PetscMatrix<libMesh::Number>& >(
			* libMesh::cast_ref<libMesh::ImplicitSystem&>(m_BIG_EquationSystem.second->get_system(type_name_BIG)).matrix);
	libMesh::PetscVector<libMesh::Number>& F_A =libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>& >(
			* libMesh::cast_ref<libMesh::ImplicitSystem&>(m_BIG_EquationSystem.second->get_system(type_name_BIG)).rhs);


	switch(m_solver_type)
	{
	case LATIN_MODIFIED_STIFFNESS:
	case LATIN_ORIGINAL_STIFFNESS:
		{
			std::shared_ptr<PETSC_LATIN_coupled_solver> cast_LATIN_solver = std::dynamic_pointer_cast<PETSC_LATIN_coupled_solver>(m_coupled_solver);

			// Set the solver parameters
			cast_LATIN_solver->set_params(k_dA,k_dB,k_cA,k_cB);

			// Set the solver matrices
			cast_LATIN_solver->set_matrices_nonlinear(M_A,C_RA,C_RB,C_RR);

			// Set the solver matrices
			cast_LATIN_solver->set_forces_nonlinear(F_A);

			// Set LATIN parameters (convergence, relaxation ... )
			cast_LATIN_solver->set_convergence_limits(eps,convIter);
			cast_LATIN_solver->set_relaxation(relax);
			break;
		}
	case CG:
		break;
	}
};

void carl::coupled_system::solve(const std::string micro_name, const std::string type_name, const std::string conv_name)
{
	// Solve!
	m_coupled_solver->solve();

	// Get the solutions
	libMesh::NumericVector<libMesh::Number>& sol_BIG =
			libMesh::cast_ref<libMesh::NumericVector<libMesh::Number>& >(m_coupled_solver->get_solution_BIG());
	libMesh::NumericVector<libMesh::Number>& sol_micro =
			libMesh::cast_ref<libMesh::NumericVector<libMesh::Number>& >(m_coupled_solver->get_solution_micro());

	// Get the systems
	libMesh::EquationSystems& EqSystems_BIG = * m_BIG_EquationSystem.second;
	libMesh::EquationSystems& EqSystems_micro = * m_micro_EquationSystemMap[micro_name];

	libMesh::System& Sys_BIG =   libMesh::cast_ref<libMesh::System&>(EqSystems_BIG.get_system(type_name));
	libMesh::System& Sys_micro = libMesh::cast_ref<libMesh::System&>(EqSystems_micro.get_system(type_name));

	// Set the solutions!
	*(Sys_BIG.solution) = sol_BIG;
	Sys_BIG.solution->close();
	Sys_BIG.update();
	*(Sys_micro.solution) = sol_micro;
	Sys_micro.solution->close();
	Sys_micro.update();

	this->print_convergence(conv_name);
}

void carl::coupled_system::solve_LATIN_nonlinear(const std::string micro_name, const std::string type_name_micro, const std::string type_name_BIG, const std::string conv_name)
{
	libMesh::EquationSystems& EqSystems_micro = * m_micro_EquationSystemMap[micro_name];
	libMesh::NonlinearImplicitSystem& Sys_micro = libMesh::cast_ref<libMesh::NonlinearImplicitSystem&>(EqSystems_micro.get_system(type_name_micro));

	// Solve!
	switch(m_solver_type)
	{
	case LATIN_MODIFIED_STIFFNESS:
	case LATIN_ORIGINAL_STIFFNESS:
		{
			std::shared_ptr<PETSC_LATIN_coupled_solver> cast_LATIN_solver = std::dynamic_pointer_cast<PETSC_LATIN_coupled_solver>(m_coupled_solver);

			cast_LATIN_solver->solve_nonlinear(EqSystems_micro,type_name_micro);
			break;
		}
	case CG:
		break;
	}


	// Get the solutions
	libMesh::NumericVector<libMesh::Number>& sol_BIG =
			libMesh::cast_ref<libMesh::NumericVector<libMesh::Number>& >(m_coupled_solver->get_solution_BIG());
	libMesh::NumericVector<libMesh::Number>& sol_micro =
			libMesh::cast_ref<libMesh::NumericVector<libMesh::Number>& >(m_coupled_solver->get_solution_micro());

	// Get the systems
	libMesh::EquationSystems& EqSystems_BIG = * m_BIG_EquationSystem.second;

	libMesh::System& Sys_BIG =   libMesh::cast_ref<libMesh::System&>(EqSystems_BIG.get_system(type_name_BIG));

	// Set the solutions!
	*(Sys_BIG.solution) = sol_BIG;
	Sys_BIG.solution->close();
	Sys_BIG.update();
	*(Sys_micro.solution) = sol_micro;
	Sys_micro.solution->close();
	Sys_micro.update();

	this->print_convergence(conv_name);

}

void carl::coupled_system::print_matrix_BIG_info(const std::string& name)
{
	libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix =
						* m_couplingMatrixMap_mediator_BIG[name];
	std::cout << "| Restrict - Macro matrix -> " << name << std::endl;
	print_matrix(CouplingTestMatrix);
}

void carl::coupled_system::print_matrices_matlab(const std::string& name, const std::string& outputRoot)
{
	libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix_BIG =
							* m_couplingMatrixMap_mediator_BIG[name];
	libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix_micro=
							* m_couplingMatrixMap_mediator_micro[name];
	libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix_mediator=
							* m_couplingMatrixMap_mediator_mediator[name];

	CouplingTestMatrix_BIG.print_matlab(outputRoot + "_BIG.m");
	CouplingTestMatrix_micro.print_matlab(outputRoot + "_micro.m");
	CouplingTestMatrix_mediator.print_matlab(outputRoot + "_mediator.m");
}

void carl::coupled_system::print_matrix_mediator_info(const std::string& name)
{
	libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix =
						* m_couplingMatrixMap_mediator_mediator[name];
	std::cout << "| Restrict - Restrict matrix -> " << name << std::endl;
	print_matrix(CouplingTestMatrix);
}

void carl::coupled_system::print_convergence(const std::string& filename)
{
	if(m_comm->rank() == 0)
	{
		std::ofstream convergenceOut(filename);

		switch(m_solver_type)
		{
		case LATIN_MODIFIED_STIFFNESS:
		case LATIN_ORIGINAL_STIFFNESS:
			{
				std::shared_ptr<PETSC_LATIN_coupled_solver> cast_LATIN_solver = std::dynamic_pointer_cast<PETSC_LATIN_coupled_solver>(m_coupled_solver);

				cast_LATIN_solver->print_convergence(convergenceOut);
				break;
			}
		case CG:
			{
				std::shared_ptr<PETSC_CG_coupled_solver> cast_CG_solver = std::dynamic_pointer_cast<PETSC_CG_coupled_solver>(m_coupled_solver);

				cast_CG_solver->print_convergence(convergenceOut);
				break;
			}
		}

		convergenceOut.close();
	}
}

void carl::coupled_system::set_rigid_body_mode(libMesh::ImplicitSystem&  input_system,
												libMesh::PetscVector<libMesh::Number>* coord_vec,
												const std::string& sys_type
												)
{
	// Set up some temporary variables to simplify code
	libMesh::PetscMatrix<libMesh::Number> * mat_sys = libMesh::cast_ptr<libMesh::PetscMatrix<libMesh::Number>* >(input_system.matrix);
	const libMesh::MeshBase& mesh_sys = input_system.get_mesh();
	unsigned int sys_sys_number = input_system.number();

	// Set vector structure as *almost* the same of rigidity matrix
	Vec PETSc_vec_sys;
	PetscInt local_N;
	MatGetLocalSize(mat_sys->mat(),NULL,&local_N);
	VecCreate(mesh_sys.comm().get(),&PETSc_vec_sys);
	VecSetSizes(PETSc_vec_sys,local_N,mat_sys->n());
	VecSetBlockSize(PETSc_vec_sys,mesh_sys.mesh_dimension());
	VecSetFromOptions(PETSc_vec_sys);

	libMesh::PetscVector<libMesh::Number> vec_sys(PETSc_vec_sys,mesh_sys.comm());
	m_coord_vect_BIG.second->init(vec_sys);

	VecDestroy(&PETSc_vec_sys);

	std::cout.flush();

	auto node_it = mesh_sys.local_nodes_begin();
	auto node_it_end = mesh_sys.local_nodes_end();

	unsigned int dof_number_BIG = 0;

	for( ; node_it != node_it_end; ++node_it)
	{
		const libMesh::Node* node_BIG = *node_it;

		for(unsigned int var=0; var<node_BIG->n_dofs(sys_sys_number); var++)
		{
			dof_number_BIG = node_BIG->dof_number(sys_sys_number,var,0);
			m_coord_vect_BIG.second->set(dof_number_BIG,node_BIG->operator ()(var));
		}
	}

	MatNullSpace nullsp_sys;
	MatNullSpaceCreateRigidBody(m_coord_vect_BIG.second->vec(),&nullsp_sys);
	MatSetNullSpace(mat_sys->mat(),nullsp_sys);
	MatNullSpaceDestroy(&nullsp_sys);
}

void carl::coupled_system::set_rigid_body_modes_BIG(const std::string& sys_name)
{

	homemade_assert_msg(m_bHasAssembled_BIG,"Macro system not assembled yet!\n");

	// Set up some temporary variables to simplify code
	libMesh::EquationSystems& EqSystems_BIG = * m_BIG_EquationSystem.second;
	libMesh::ImplicitSystem& Sys_BIG =
			libMesh::cast_ref<libMesh::ImplicitSystem&>(EqSystems_BIG.get_system(sys_name));

	set_rigid_body_mode(Sys_BIG,m_coord_vect_BIG.second,m_BIG_EquationSystem.first);
}

void carl::coupled_system::set_rigid_body_modes_micro(const std::string micro_name, const std::string& sys_name)
{

	homemade_assert_msg(m_bHasAssembled_micro[micro_name],"Micro system not assembled yet!\n");

	// Set up some temporary variables to simplify code
	libMesh::EquationSystems& EqSystems_micro = * m_micro_EquationSystemMap[micro_name];
	libMesh::ImplicitSystem& Sys_micro =
			libMesh::cast_ref<libMesh::ImplicitSystem&>(EqSystems_micro.get_system(sys_name));

	set_rigid_body_mode(Sys_micro,m_coord_vect_microMap[micro_name],micro_name);
}


void carl::coupled_system::print_perf_log(std::string filename_input)
{
	m_coupled_solver->print_perf_log(filename_input);
}
