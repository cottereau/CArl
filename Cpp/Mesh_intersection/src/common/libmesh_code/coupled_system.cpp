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
	homemade_assert_msg(base_elem->type() == libMesh::TET4 || base_elem->type() == libMesh::TET10,
			" Only implemented for TET4 or TET10 elements!");

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
	for(int iii = 0; iii < lambda_base_size; ++iii)
	{
		for(int jjj = 0; jjj < lambda_inter_size; ++jjj)
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

void carl::coupled_system::set_LATIN_solver(const std::string micro_name,
		const std::string type_name, double k_dA, double k_dB, double k_cA, double k_cB,
		double eps, int convIter, double relax)
{
	// Get the systems
	libMesh::EquationSystems& EqSystems_BIG = * m_BIG_EquationSystem.second;
	libMesh::EquationSystems& EqSystems_micro = * m_micro_EquationSystemMap[micro_name];

	libMesh::ImplicitSystem& Sys_BIG =   libMesh::cast_ref<libMesh::ImplicitSystem&>(EqSystems_BIG.get_system(type_name));
	libMesh::ImplicitSystem& Sys_micro = libMesh::cast_ref<libMesh::ImplicitSystem&>(EqSystems_micro.get_system(type_name));

	// Assemble the systems
	Sys_BIG.assemble();
	Sys_BIG.matrix->close();
	Sys_BIG.rhs->close();
	m_bHasAssembled_BIG = true;

	Sys_micro.assemble();
	Sys_micro.matrix->close();
	Sys_micro.rhs->close();
	m_bHasAssembled_micro[micro_name] = true;

	// Get the matrices
	libMesh::PetscMatrix<libMesh::Number>& M_A = libMesh::cast_ref<libMesh::PetscMatrix<libMesh::Number>& >(* Sys_BIG.matrix);
	libMesh::PetscMatrix<libMesh::Number>& M_B = libMesh::cast_ref<libMesh::PetscMatrix<libMesh::Number>& >(* Sys_micro.matrix);

	libMesh::PetscMatrix<libMesh::Number>& C_RA = * m_couplingMatrixMap_mediator_BIG[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& C_RB = * m_couplingMatrixMap_mediator_micro[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& C_RR = * m_couplingMatrixMap_mediator_mediator[micro_name];

	// Get the vectors
	libMesh::PetscVector<libMesh::Number>& F_A = libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>& >(* Sys_BIG.rhs);
	libMesh::PetscVector<libMesh::Number>& F_B = libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>& >(* Sys_micro.rhs);

	// Set the solver parameters
	m_LATIN_solver.set_params(k_dA,k_dB,k_cA,k_cB);

	// Set the solver matrices
	m_LATIN_solver.set_matrices(M_A,M_B,C_RA,C_RB,C_RR);

	// Set the solver matrices
	m_LATIN_solver.set_forces(F_A,F_B);

	// Set LATIN parameters (convergence, relaxation ... )
	m_LATIN_solver.set_convergence_limits(eps,convIter);
	m_LATIN_solver.set_relaxation(relax);
};

void carl::coupled_system::set_LATIN_solver(const std::string micro_name, const std::string type_name,
												void fptr_BIG(		libMesh::EquationSystems& es,
																	const std::string& name, weight_parameter_function& weight_mask),
												void fptr_micro(	libMesh::EquationSystems& es,
																	const std::string& name, weight_parameter_function& weight_mask),
												double k_dA, double k_dB, double k_cA, double k_cB,
												double eps, int convIter, double relax)
{
	// Get the systems
	libMesh::EquationSystems& EqSystems_BIG = * m_BIG_EquationSystem.second;
	libMesh::EquationSystems& EqSystems_micro = * m_micro_EquationSystemMap[micro_name];

	libMesh::ImplicitSystem& Sys_BIG =   libMesh::cast_ref<libMesh::ImplicitSystem&>(EqSystems_BIG.get_system(type_name));
	libMesh::ImplicitSystem& Sys_micro = libMesh::cast_ref<libMesh::ImplicitSystem&>(EqSystems_micro.get_system(type_name));

	// Get the weight functions
	weight_parameter_function& alpha_mask = * m_alpha_masks[micro_name];

	// Assemble the systems
	fptr_BIG(EqSystems_BIG,type_name,alpha_mask);
	Sys_BIG.matrix->close();
	Sys_BIG.rhs->close();
	m_bHasAssembled_BIG = true;

	fptr_micro(EqSystems_micro,type_name,alpha_mask);
	Sys_micro.matrix->close();
	Sys_micro.rhs->close();
	m_bHasAssembled_micro[micro_name] = true;

	// Get the matrices
	libMesh::PetscMatrix<libMesh::Number>& M_A = libMesh::cast_ref<libMesh::PetscMatrix<libMesh::Number>& >(* Sys_BIG.matrix);
	libMesh::PetscMatrix<libMesh::Number>& M_B = libMesh::cast_ref<libMesh::PetscMatrix<libMesh::Number>& >(* Sys_micro.matrix);

	libMesh::PetscMatrix<libMesh::Number>& C_RA = * m_couplingMatrixMap_mediator_BIG[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& C_RB = * m_couplingMatrixMap_mediator_micro[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& C_RR = * m_couplingMatrixMap_mediator_mediator[micro_name];

	// Get the vectors
	libMesh::PetscVector<libMesh::Number>& F_A = libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>& >(* Sys_BIG.rhs);
	libMesh::PetscVector<libMesh::Number>& F_B = libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>& >(* Sys_micro.rhs);

	// Set the solver parameters
	m_LATIN_solver.set_params(k_dA,k_dB,k_cA,k_cB);

	// Set the solver matrices
	m_LATIN_solver.set_matrices(M_A,M_B,C_RA,C_RB,C_RR);

	// Set the solver matrices
	m_LATIN_solver.set_forces(F_A,F_B);

	// Set LATIN parameters (convergence, relaxation ... )
	m_LATIN_solver.set_convergence_limits(eps,convIter);
	m_LATIN_solver.set_relaxation(relax);
};

void carl::coupled_system::solve_LATIN(const std::string micro_name, const std::string type_name, const std::string conv_name)
{
	// Solve!
	m_LATIN_solver.solve();

	// Get the solutions
	libMesh::NumericVector<libMesh::Number>& sol_BIG =
			libMesh::cast_ref<libMesh::NumericVector<libMesh::Number>& >(m_LATIN_solver.get_solution_BIG());
	libMesh::NumericVector<libMesh::Number>& sol_micro =
			libMesh::cast_ref<libMesh::NumericVector<libMesh::Number>& >(m_LATIN_solver.get_solution_micro());

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

	print_LATIN_convergence(conv_name);
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

void carl::coupled_system::print_LATIN_convergence(const std::string& filename)
{
	std::ofstream convergenceOut(filename);
	m_LATIN_solver.print_convergence(convergenceOut);
	convergenceOut.close();

}
