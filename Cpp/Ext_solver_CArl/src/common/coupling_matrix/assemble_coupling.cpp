#include "assemble_coupling.h"

namespace carl
{
	// Class members
	void assemble_coupling_matrices::clear()
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
	};

	libMesh::EquationSystems& assemble_coupling_matrices::set_BIG_EquationSystem(const std::string& name, libMesh::MeshBase& BIGMesh)
	{
		libMesh::EquationSystems* EqSystemPtr = new libMesh::EquationSystems(
				BIGMesh);

		m_BIG_EquationSystem.first = name;
		m_BIG_EquationSystem.second = EqSystemPtr;

		return *EqSystemPtr;
	};

	libMesh::EquationSystems& assemble_coupling_matrices::set_Restricted_BIG_EquationSystem(const std::string& name, libMesh::MeshBase& R_BIGMesh)
	{
		m_bHasDefinedMeshRestrictions = true;
		libMesh::EquationSystems* EqSystemPtr = new libMesh::EquationSystems(
				R_BIGMesh);

		m_R_BIG_EquationSystem.first = name;
		m_R_BIG_EquationSystem.second = EqSystemPtr;

		return *EqSystemPtr;
	};

	libMesh::EquationSystems& assemble_coupling_matrices::add_Restricted_micro_EquationSystem(const std::string& name,
			libMesh::MeshBase& R_microMesh)
	{
		m_bHasDefinedMeshRestrictions = true;

		libMesh::EquationSystems* EqSystemPtr = NULL;

		if (!m_R_micro_EquationSystemMap.count(name))
		{
			// Then add a new system
			EqSystemPtr = new libMesh::EquationSystems(R_microMesh);

			m_R_micro_EquationSystemMap.insert(std::make_pair(name, EqSystemPtr));
		}
		else
		{
			// System already exits, return the system pointer
			std::cerr << " *** Warning: restricted micro system " << name
					<< " already exists!" << std::endl;
			EqSystemPtr = m_micro_EquationSystemMap[name];
		}

		return *EqSystemPtr;
	};

	libMesh::EquationSystems& assemble_coupling_matrices::add_inter_EquationSystem(const std::string& name,
			libMesh::MeshBase& interMesh)
	{
		libMesh::EquationSystems* EqSystemPtr = NULL;

		if (!m_inter_EquationSystemMap.count(name))
		{
			// Then add a new system
			EqSystemPtr = new libMesh::EquationSystems(interMesh);
			m_inter_EquationSystemMap.insert(std::make_pair(name, EqSystemPtr));
		}
		else
		{
			// System already exits, return the system pointer
			std::cerr << " *** Warning: inter system " << name
					<< " already exists!" << std::endl;
			EqSystemPtr = m_inter_EquationSystemMap[name];
		}

		return *EqSystemPtr;
	};

	libMesh::EquationSystems& assemble_coupling_matrices::add_mediator_EquationSystem(
			const std::string& name, libMesh::MeshBase& mediatorMesh)
	{
		libMesh::EquationSystems* EqSystemPtr = NULL;

		if (!m_mediator_EquationSystemMap.count(name))
		{
			// Then add a new system
			EqSystemPtr = new libMesh::EquationSystems(mediatorMesh);
			m_mediator_EquationSystemMap.insert(
					std::make_pair(name, EqSystemPtr));
		}
		else
		{
			// System already exits, return the system pointer
			std::cerr << " *** Warning: mediator system " << name
					<< " already exists!" << std::endl;
			EqSystemPtr = m_mediator_EquationSystemMap[name];
		}

		return *EqSystemPtr;
	};

	void assemble_coupling_matrices::set_coupling_parameters(const std::string& name,
			double coupling_rigidity, double coupling_width)
	{
		m_coupling_rigidityMap[name] = coupling_rigidity;
		m_coupling_widthMap[name] = coupling_width;
	};

	void assemble_coupling_matrices::set_corrected_shapes(	const std::vector<std::vector<libMesh::Real> >& 	lambda_weights,
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

	void assemble_coupling_matrices::set_corrected_shape_gradients(	const std::vector<std::vector<libMesh::Real> >& 	lambda_weights,
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

	void assemble_coupling_matrices::get_lambdas(	const unsigned int 							dim,
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

	libMesh::PetscMatrix<libMesh::Number>& assemble_coupling_matrices::get_micro_coupling_matrix(const std::string& name)
	{
		return * m_couplingMatrixMap_mediator_micro[name];
	}

	libMesh::PetscMatrix<libMesh::Number>& assemble_coupling_matrices::get_BIG_coupling_matrix(const std::string& name)
	{
		return * m_couplingMatrixMap_mediator_BIG[name];
	}

	libMesh::PetscMatrix<libMesh::Number>& assemble_coupling_matrices::get_mediator_coupling_matrix(const std::string& name)
	{
		return * m_couplingMatrixMap_mediator_mediator[name];
	}

	void assemble_coupling_matrices::print_matrix_micro_info(const std::string& name)
	{
		libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix =
							* m_couplingMatrixMap_mediator_micro[name];
		std::cout << "| Restrict - Micro matrix -> " << name << std::endl;
		print_matrix(CouplingTestMatrix);
	}

	void assemble_coupling_matrices::print_matrix_BIG_info(const std::string& name)
	{
		libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix =
							* m_couplingMatrixMap_mediator_BIG[name];
		std::cout << "| Restrict - Macro matrix -> " << name << std::endl;
		print_matrix(CouplingTestMatrix);
	}

	void assemble_coupling_matrices::print_matrix_mediator_info(const std::string& name)
	{
		libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix =
							* m_couplingMatrixMap_mediator_mediator[name];
		std::cout << "| Restrict - Restrict matrix -> " << name << std::endl;
		print_matrix(CouplingTestMatrix);
	}

	void assemble_coupling_matrices::print_matrices_matlab(const std::string& name, const std::string& outputRoot)
	{
		libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix_BIG =
								* m_couplingMatrixMap_mediator_BIG[name];
		libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix_micro=
								* m_couplingMatrixMap_mediator_micro[name];
		libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix_mediator=
								* m_couplingMatrixMap_mediator_mediator[name];

		CouplingTestMatrix_BIG.print_matlab(outputRoot + "_macro.m");
		CouplingTestMatrix_micro.print_matlab(outputRoot + "_micro.m");
		CouplingTestMatrix_mediator.print_matlab(outputRoot + "_mediator.m");
	}

	void assemble_coupling_matrices::print_PETSC_matrices(const std::string& name, const std::string& outputRoot)
	{
		libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix_BIG =
								* m_couplingMatrixMap_mediator_BIG[name];
		libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix_micro=
								* m_couplingMatrixMap_mediator_micro[name];
		libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix_mediator=
								* m_couplingMatrixMap_mediator_mediator[name];

		write_PETSC_matrix(CouplingTestMatrix_BIG.mat(), outputRoot + "_macro.petscmat",m_comm->get());
		write_PETSC_matrix(CouplingTestMatrix_micro.mat(), outputRoot + "_micro.petscmat",m_comm->get());
		write_PETSC_matrix(CouplingTestMatrix_mediator.mat(), outputRoot + "_mediator.petscmat",m_comm->get());
	}

	void assemble_coupling_matrices::use_H1_coupling(std::string name)
	{
		m_bUseH1Coupling[name] = true;
	};

	void assemble_coupling_matrices::use_L2_coupling(std::string name)
	{
		m_bUseH1Coupling[name] = false;
	};
};