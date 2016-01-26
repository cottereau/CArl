#include <petscmat.h>
#include <petscsys.h>
#include "coupled_system.h"

// Functions

void carl::lump_matrix(		libMesh::PetscMatrix<libMesh::Number>& matrixInput,
							libMesh::PetscMatrix<libMesh::Number>& matrixOutput)
{
	if(matrixOutput.initialized())
	{
		matrixOutput.clear();
	}

	int M = matrixInput.m();
	int N = matrixInput.n();

	PetscInt local_M, local_N;

	MatGetLocalSize(matrixInput.mat(),&local_M,&local_N);

	// It will be a diagonal matrix, so no need of a heavy preallocation
	matrixOutput.init(M,N,local_M,local_N,1,0);

	libMesh::PetscVector<libMesh::Number> UnityVec(matrixInput.comm(),M,local_M);
	libMesh::PetscVector<libMesh::Number> DummyVector(matrixInput.comm(),M,local_M);

	VecSet(UnityVec.vec(),1);

	UnityVec.close();

	matrixInput.vector_mult(DummyVector,UnityVec);

	MatDiagonalSet(matrixOutput.mat(),DummyVector.vec(),INSERT_VALUES);

	matrixOutput.close();
}

void carl::print_matrix(libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix)
{
	libMesh::Real accumulator = 0;
	std::cout << "| M i,j : " << CouplingTestMatrix.m() << " x " << CouplingTestMatrix.n() << std::endl;
	for(unsigned int iii = 0; iii < CouplingTestMatrix.m(); ++iii)
	{
		for(unsigned int jjj = 0; jjj < CouplingTestMatrix.n(); ++jjj)
		{
			accumulator += CouplingTestMatrix(iii,jjj);
			if(CouplingTestMatrix.m() < 11 && CouplingTestMatrix.n() < 11)
			{
				std::cout << " " << CouplingTestMatrix(iii,jjj);
			}
		}
		if(CouplingTestMatrix.m() < 11 && CouplingTestMatrix.n() < 11)
		{
			std::cout << std::endl;
		}
	}
	std::cout << "|" << std::endl;
	std::cout << "| Sum( C_i,j ) = " << accumulator << std::endl << std::endl;
}


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

	while(!m_inter_EquationSystemMap.empty())
	{
		EqSystem_iterator toClean = m_inter_EquationSystemMap.begin();

		libMesh::EquationSystems *EqSys = toClean->second;
		EqSys->clear();
		delete EqSys;
		EqSys = NULL;

		m_inter_EquationSystemMap.erase(toClean);
	}

	while(!m_restrict_EquationSystemMap.empty())
	{
		EqSystem_iterator toClean = m_restrict_EquationSystemMap.begin();

		libMesh::EquationSystems *EqSys = toClean->second;
		EqSys->clear();
		delete EqSys;
		EqSys = NULL;

		m_restrict_EquationSystemMap.erase(toClean);
	}

	while(!m_couplingMatrixMap_restrict_micro.empty())
	{
		Matrix_iterator toClean = m_couplingMatrixMap_restrict_micro.begin();

		libMesh::PetscMatrix<libMesh::Number> *Mat = toClean->second;
		Mat->clear();
		delete Mat;
		Mat = NULL;

		m_couplingMatrixMap_restrict_micro.erase(toClean);
	}

	while(!m_couplingMatrixMap_restrict_BIG.empty())
	{
		Matrix_iterator toClean = m_couplingMatrixMap_restrict_BIG.begin();

		libMesh::PetscMatrix<libMesh::Number> *Mat = toClean->second;
		Mat->clear();
		delete Mat;
		Mat = NULL;

		m_couplingMatrixMap_restrict_BIG.erase(toClean);
	}

	while(!m_couplingMatrixMap_restrict_restrict.empty())
	{
		Matrix_iterator toClean = m_couplingMatrixMap_restrict_restrict.begin();

		libMesh::PetscMatrix<libMesh::Number> *Mat = toClean->second;
		Mat->clear();
		delete Mat;
		Mat = NULL;

		m_couplingMatrixMap_restrict_restrict.erase(toClean);
	}
};

void carl::coupled_system::assemble_coupling_matrices(	const std::string BIG_name,
														const std::string micro_name,
														const std::string inter_name,
														const std::string restrict_name,
														std::unordered_map<int,int>& equivalence_table_restrict_BIG,
														std::vector<std::pair<int,int> >& intersection_table_restrict_micro,
														std::unordered_multimap<int,int>& intersection_table_inter,
														bool using_same_mesh_restrict_A,
														bool bSameElemsType)
{
	// Addresses to the eq. systems
	libMesh::EquationSystems& restrict_eq_system = * m_restrict_EquationSystemMap[restrict_name];

	libMesh::EquationSystems& BIG_eq_system = * m_BIG_EquationSystem.second;
	libMesh::EquationSystems& micro_eq_system = * m_micro_EquationSystemMap[micro_name];
	libMesh::EquationSystems& inter_eq_system = * m_inter_EquationSystemMap[inter_name];

	// Addresses to the matrices
	libMesh::PetscMatrix<libMesh::Number>& couplingMatrix_restrict_BIG = * m_couplingMatrixMap_restrict_BIG[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& couplingMatrix_restrict_micro = * m_couplingMatrixMap_restrict_micro[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& couplingMatrix_restrict_restrict = * m_couplingMatrixMap_restrict_restrict[micro_name];

	// Addresses to the meshes
	const libMesh::MeshBase& mesh_restrict = restrict_eq_system.get_mesh();

	const libMesh::MeshBase& mesh_BIG = BIG_eq_system.get_mesh();
	const libMesh::MeshBase& mesh_micro = micro_eq_system.get_mesh();
	const libMesh::MeshBase& mesh_inter = inter_eq_system.get_mesh();

	const unsigned int dim_restrict = mesh_restrict.mesh_dimension();

	const unsigned int dim_BIG = mesh_BIG.mesh_dimension();
	const unsigned int dim_micro = mesh_micro.mesh_dimension();
	const unsigned int dim_inter = mesh_inter.mesh_dimension();

	// Systems and vars
	libMesh::LinearImplicitSystem& volume_restrict_system
		= restrict_eq_system.get_system<libMesh::LinearImplicitSystem>("VolTest");

	libMesh::LinearImplicitSystem& volume_BIG_system
		= BIG_eq_system.get_system<libMesh::LinearImplicitSystem>("VolTest");
	libMesh::LinearImplicitSystem& volume_micro_system
		= micro_eq_system.get_system<libMesh::LinearImplicitSystem>("VolTest");
	libMesh::LinearImplicitSystem& volume_inter_system
		= inter_eq_system.get_system<libMesh::LinearImplicitSystem>("VolTest");

	const unsigned int n_components = 1;
	const unsigned int vol_test_restrict = volume_restrict_system.variable_number("SillyVar");

	const unsigned int vol_test_BIG = volume_BIG_system.variable_number ("SillyVar");
	const unsigned int vol_test_micro = volume_micro_system.variable_number ("SillyVar");
	const unsigned int vol_test_inter = volume_inter_system.variable_number ("SillyVar");

	// DoF maps
	const libMesh::DofMap& dof_map_restrict = volume_restrict_system.get_dof_map();

	const libMesh::DofMap& dof_map_BIG = volume_BIG_system.get_dof_map();
	const libMesh::DofMap& dof_map_micro = volume_micro_system.get_dof_map();
	const libMesh::DofMap& dof_map_inter = volume_inter_system.get_dof_map();

	libMesh::FEType fe_type_restrict = dof_map_restrict.variable_type(vol_test_restrict);

	libMesh::FEType fe_type_BIG = dof_map_BIG.variable_type(vol_test_BIG);
	libMesh::FEType fe_type_micro = dof_map_micro.variable_type(vol_test_micro);
	libMesh::FEType fe_type_inter = dof_map_inter.variable_type(vol_test_inter);

	// Set up FE bases and quadratures
	libMesh::UniquePtr<libMesh::FEBase> fe_restrict (libMesh::FEBase::build(dim_restrict, fe_type_restrict));

	libMesh::UniquePtr<libMesh::FEBase> fe_BIG (libMesh::FEBase::build(dim_BIG, fe_type_BIG));
	libMesh::UniquePtr<libMesh::FEBase> fe_micro (libMesh::FEBase::build(dim_micro, fe_type_micro));
	libMesh::UniquePtr<libMesh::FEBase> fe_inter (libMesh::FEBase::build(dim_inter, fe_type_inter));

	libMesh::QGauss qrule_restrict(dim_restrict, fe_type_restrict.default_quadrature_order());
	fe_restrict->attach_quadrature_rule (&qrule_restrict);

	libMesh::QGauss qrule_inter(dim_inter, fe_type_inter.default_quadrature_order());
	fe_inter->attach_quadrature_rule (&qrule_inter);

	// Jacobians
	const std::vector<libMesh::Real>& JxWrestrict = fe_restrict->get_JxW();
	const std::vector<libMesh::Real>& JxW = fe_inter->get_JxW();

	// Shape functions
	const std::vector<std::vector<libMesh::Real> >& phi_restrict = fe_restrict->get_phi();

	const std::vector<std::vector<libMesh::Real> >& phi_BIG = fe_BIG->get_phi();
	const std::vector<std::vector<libMesh::Real> >& phi_micro = fe_micro->get_phi();

	std::vector<libMesh::dof_id_type> dof_indices_restrict;

	std::vector<libMesh::dof_id_type> dof_indices_BIG;
	std::vector<libMesh::dof_id_type> dof_indices_micro;

	// Local matrix
	libMesh::DenseMatrix<libMesh::Number> Me_micro;
	libMesh::DenseMatrix<libMesh::Number> Me_BIG;

	unsigned int n_dofs_restrict;
	unsigned int n_dofs_micro;
	unsigned int n_dofs_BIG;

	if(bSameElemsType)
	{
		const libMesh::Elem* elem_restrict = mesh_restrict.elem(0);
		dof_map_restrict.dof_indices(elem_restrict, dof_indices_restrict);
		n_dofs_restrict   = dof_indices_restrict.size();

		const libMesh::Elem* elem_BIG = mesh_BIG.elem(0);
		dof_map_BIG.dof_indices(elem_BIG, dof_indices_BIG);
		n_dofs_BIG   = dof_indices_BIG.size();

		const libMesh::Elem* elem_micro = mesh_micro.elem(0);
		dof_map_micro.dof_indices(elem_micro, dof_indices_micro);
		n_dofs_micro   = dof_indices_micro.size();

		Me_BIG.resize (n_dofs_restrict, n_dofs_BIG);
		Me_micro.resize (n_dofs_restrict, n_dofs_micro);
	}

	// Vector that will keep the quadrature points
	const std::vector<libMesh::Point>& quad_points_inter = fe_inter->get_xyz();
	const std::vector<libMesh::Point>& quad_points_restrict = fe_restrict->get_xyz();

	// Initialize global matrix
	const unsigned int restrict_M = dof_map_restrict.n_dofs();

	const unsigned int BIG_N = dof_map_BIG.n_dofs();
	const unsigned int micro_N = dof_map_micro.n_dofs();

	const unsigned int restrict_M_local = dof_map_restrict.n_local_dofs();

	const unsigned int BIG_N_local = dof_map_BIG.n_local_dofs();
	const unsigned int micro_N_local = dof_map_micro.n_local_dofs();

	const std::vector<libMesh::dof_id_type>& n_nz = dof_map_restrict.get_n_nz();
	const std::vector<libMesh::dof_id_type>& n_oz = dof_map_restrict.get_n_oz();

	const std::vector<libMesh::dof_id_type>& micro_n_nz = dof_map_micro.get_n_nz();
	const std::vector<libMesh::dof_id_type>& micro_n_oz = dof_map_micro.get_n_oz();

	int max_nnz = *std::max_element(n_nz.begin(),n_nz.end());
	int max_noz = *std::max_element(n_oz.begin(),n_oz.end());

	int micro_max_nnz = *std::max_element(micro_n_nz.begin(),micro_n_nz.end());
	int micro_max_noz = *std::max_element(micro_n_oz.begin(),micro_n_oz.end());

	// TODO : correct memory allocation
	int temp_diag_alloc = max_nnz*micro_max_nnz;
	int temp_off_alloc = max_noz*micro_max_noz + max_nnz*micro_max_noz + max_noz*micro_max_nnz;
	couplingMatrix_restrict_micro.init(		restrict_M, micro_N,
											restrict_M_local, micro_N_local,
											temp_diag_alloc,temp_off_alloc);
	couplingMatrix_restrict_BIG.init(		restrict_M, BIG_N,
											restrict_M_local, BIG_N_local,
											temp_diag_alloc,temp_off_alloc);

	couplingMatrix_restrict_restrict.attach_dof_map(dof_map_restrict);
	couplingMatrix_restrict_restrict.init();

	// Intersection indexes and iterators
	int nb_of_intersections = intersection_table_restrict_micro.size();
	int elem_restrict_idx = -1;

	int elem_BIG_idx = -1;
	int elem_micro_idx = -1;
	int elem_inter_idx = -1;
	std::unordered_multimap<int,int>::iterator it_inter_idx;

	// For each intersection
	for(int iii = 0; iii < nb_of_intersections; ++iii)
	{
		// Get the restrict and micro element pointers
		elem_restrict_idx = intersection_table_restrict_micro[iii].first;
		elem_micro_idx = intersection_table_restrict_micro[iii].second;
		if(using_same_mesh_restrict_A)
		{
			elem_BIG_idx = elem_restrict_idx;
		}
		else
		{
			elem_BIG_idx = equivalence_table_restrict_BIG[elem_restrict_idx];
		}

		const libMesh::Elem* elem_restrict = mesh_restrict.elem(elem_restrict_idx);

		const libMesh::Elem* elem_BIG = mesh_BIG.elem(elem_BIG_idx);
		const libMesh::Elem* elem_micro = mesh_micro.elem(elem_micro_idx);

		// And their dof map indices
		dof_map_restrict.dof_indices(elem_restrict, dof_indices_restrict);
		dof_map_micro.dof_indices(elem_micro, dof_indices_micro);
		if(using_same_mesh_restrict_A)
		{
			dof_indices_BIG = dof_indices_restrict;
		}
		else
		{
			dof_map_BIG.dof_indices(elem_BIG, dof_indices_BIG);
		}

		// Resize dense matrix, if needed
		if(!bSameElemsType)
		{
			// Then must resize the matrix
			n_dofs_restrict   = dof_indices_restrict.size();

			n_dofs_BIG   = dof_indices_BIG.size();
			n_dofs_micro   = dof_indices_micro.size();

			Me_micro.resize (n_dofs_restrict, n_dofs_micro);
			Me_BIG.resize (n_dofs_restrict, n_dofs_BIG);
		}

		Me_micro.zero();
		Me_BIG.zero();

		// Now iterate over the intersections
		auto inter_idx_range = intersection_table_inter.equal_range(iii);

		for(	it_inter_idx = inter_idx_range.first;
				it_inter_idx != inter_idx_range.second;
				++it_inter_idx)
		{
			// Get the intersection mesh pointer
			elem_inter_idx = it_inter_idx->second;

			const libMesh::Elem* elem_inter = mesh_inter.elem(elem_inter_idx);

			// Restart the elements
			fe_inter->reinit(elem_inter);

			fe_restrict->reinit(elem_restrict,&quad_points_inter);
			fe_micro->reinit(elem_micro,&quad_points_inter);

			// For each quadrature point determinate the sub-matrices elements
			for (unsigned int qp=0; qp<qrule_inter.n_points(); qp++)
			{
				// Restrict -> micro coupling
				L2_Coupling(Me_micro,qp,phi_restrict,phi_micro,
							n_dofs_restrict,n_dofs_micro,JxW);

				// Restrict -> restrict coupling
				L2_Coupling(Me_BIG,qp,phi_restrict,phi_restrict,
							n_dofs_restrict,n_dofs_restrict,JxW);
			}
		}

		couplingMatrix_restrict_micro.add_matrix(Me_micro,dof_indices_restrict,dof_indices_micro);
		couplingMatrix_restrict_BIG.add_matrix(Me_BIG,dof_indices_restrict,dof_indices_BIG);
		couplingMatrix_restrict_restrict.add_matrix(Me_BIG,dof_indices_restrict,dof_indices_restrict);
	}

	couplingMatrix_restrict_micro.close();
	couplingMatrix_restrict_BIG.close();
	couplingMatrix_restrict_restrict.close();
};

libMesh::PetscMatrix<libMesh::Number>& carl::coupled_system::get_micro_coupling_matrix(const std::string& name)
{
	return * m_couplingMatrixMap_restrict_micro[name];
}

libMesh::PetscMatrix<libMesh::Number>& carl::coupled_system::get_BIG_coupling_matrix(const std::string& name)
{
	return * m_couplingMatrixMap_restrict_BIG[name];
}

libMesh::PetscMatrix<libMesh::Number>& carl::coupled_system::get_restrict_coupling_matrix(const std::string& name)
{
	return * m_couplingMatrixMap_restrict_restrict[name];
}

void carl::coupled_system::print_matrix_micro_info(const std::string& name)
{
	libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix =
						* m_couplingMatrixMap_restrict_micro[name];
	std::cout << "| Restrict - Micro matrix -> " << name << std::endl;
	print_matrix(CouplingTestMatrix);
}

void carl::coupled_system::print_matrix_BIG_info(const std::string& name)
{
	libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix =
						* m_couplingMatrixMap_restrict_BIG[name];
	std::cout << "| Restrict - Macro matrix -> " << name << std::endl;
	print_matrix(CouplingTestMatrix);
}

void carl::coupled_system::print_matrix_restrict_info(const std::string& name)
{
	libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix =
						* m_couplingMatrixMap_restrict_restrict[name];
	std::cout << "| Restrict - Restrict matrix -> " << name << std::endl;
	print_matrix(CouplingTestMatrix);
}
