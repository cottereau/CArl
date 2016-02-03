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

	while(!m_alpha_masks.empty())
	{
		alpha_mask_iterator toClean = m_alpha_masks.begin();

		weight_parameter_function *alpha = toClean->second;
		alpha->clear();
		delete alpha;
		alpha = NULL;
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

	// TODO : correct memory allocation
	couplingMatrix_restrict_micro.init(		restrict_M, micro_N,
											restrict_M_local, micro_N_local,
											micro_N_local,micro_N - micro_N_local);
	couplingMatrix_restrict_BIG.init(		restrict_M, BIG_N,
											restrict_M_local, BIG_N_local,
											BIG_N_local,BIG_N - BIG_N_local);

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
							n_dofs_restrict,n_dofs_micro,JxW,1);

				// Restrict -> restrict coupling
				L2_Coupling(Me_BIG,qp,phi_restrict,phi_restrict,
							n_dofs_restrict,n_dofs_restrict,JxW,1);
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

void carl::coupled_system::assemble_coupling_elasticity_3D(	const std::string BIG_name,
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

	// First, test if all the systems have an elasticity model and variable sety
	libmesh_assert_msg(micro_eq_system.has_system("Elasticity"),
				" Micro equation systems missing \"Elasticity\" system!");
	libmesh_assert_msg(BIG_eq_system.has_system("Elasticity"),
				" Macro equation systems missing \"Elasticity\" system!");
	libmesh_assert_msg(inter_eq_system.has_system("Elasticity"),
				" Intersection equation systems missing \"Elasticity\" system!");
	libmesh_assert_msg(restrict_eq_system.has_system("Elasticity"),
				" Restricted equation systems missing \"Elasticity\" system!");

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
		= restrict_eq_system.get_system<libMesh::LinearImplicitSystem>("Elasticity");
	libMesh::LinearImplicitSystem& volume_BIG_system
		= BIG_eq_system.get_system<libMesh::LinearImplicitSystem>("Elasticity");
	libMesh::LinearImplicitSystem& volume_micro_system
		= micro_eq_system.get_system<libMesh::LinearImplicitSystem>("Elasticity");
	libMesh::LinearImplicitSystem& volume_inter_system
		= inter_eq_system.get_system<libMesh::LinearImplicitSystem>("Elasticity");

	const unsigned int u_var_restrict = volume_restrict_system.variable_number("u");
	const unsigned int v_var_restrict = volume_restrict_system.variable_number("v");
	const unsigned int w_var_restrict = volume_restrict_system.variable_number("w");

	const unsigned int u_var_BIG = volume_BIG_system.variable_number("u");
	const unsigned int v_var_BIG = volume_BIG_system.variable_number("v");
	const unsigned int w_var_BIG = volume_BIG_system.variable_number("w");

	const unsigned int u_var_micro = volume_micro_system.variable_number("u");
	const unsigned int v_var_micro = volume_micro_system.variable_number("v");
	const unsigned int w_var_micro = volume_micro_system.variable_number("w");

	const unsigned int u_var_inter = volume_inter_system.variable_number ("u");
	const unsigned int v_var_inter = volume_inter_system.variable_number ("v");
	const unsigned int w_var_inter = volume_inter_system.variable_number ("w");

	// DoF maps
	const libMesh::DofMap& dof_map_restrict = volume_restrict_system.get_dof_map();
	const libMesh::DofMap& dof_map_BIG = volume_BIG_system.get_dof_map();
	const libMesh::DofMap& dof_map_micro = volume_micro_system.get_dof_map();
	const libMesh::DofMap& dof_map_inter = volume_inter_system.get_dof_map();

	libMesh::FEType fe_type_restrict = dof_map_restrict.variable_type(u_var_restrict);
	libMesh::FEType fe_type_BIG = dof_map_BIG.variable_type(u_var_BIG);
	libMesh::FEType fe_type_micro = dof_map_micro.variable_type(u_var_micro);
	libMesh::FEType fe_type_inter = dof_map_inter.variable_type(u_var_inter);

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

	// Dof vectors and ranges
	std::vector<libMesh::dof_id_type> dof_indices_restrict;
	std::vector<libMesh::dof_id_type> dof_indices_u_restrict;
	std::vector<libMesh::dof_id_type> dof_indices_v_restrict;
	std::vector<libMesh::dof_id_type> dof_indices_w_restrict;

	std::vector<libMesh::dof_id_type> dof_indices_BIG;
	std::vector<libMesh::dof_id_type> dof_indices_u_BIG;
	std::vector<libMesh::dof_id_type> dof_indices_v_BIG;
	std::vector<libMesh::dof_id_type> dof_indices_w_BIG;

	std::vector<libMesh::dof_id_type> dof_indices_micro;
	std::vector<libMesh::dof_id_type> dof_indices_u_micro;
	std::vector<libMesh::dof_id_type> dof_indices_v_micro;
	std::vector<libMesh::dof_id_type> dof_indices_w_micro;

	unsigned int n_dofs_restrict;
	unsigned int n_dofs_u_restrict;
	unsigned int n_dofs_v_restrict;
	unsigned int n_dofs_w_restrict;

	unsigned int n_dofs_micro;
	unsigned int n_dofs_u_micro;
	unsigned int n_dofs_v_micro;
	unsigned int n_dofs_w_micro;

	unsigned int n_dofs_BIG;
	unsigned int n_dofs_u_BIG;
	unsigned int n_dofs_v_BIG;
	unsigned int n_dofs_w_BIG;

	// Local matrix
	libMesh::DenseMatrix<libMesh::Number> Me_micro;
	libMesh::DenseMatrix<libMesh::Number> Me_BIG;

	libMesh::DenseSubMatrix<libMesh::Number>
	M_micro_uu(Me_micro), M_micro_vv(Me_micro), M_micro_ww(Me_micro);

	libMesh::DenseSubMatrix<libMesh::Number>
	M_micro_uv(Me_micro), M_micro_vw(Me_micro), M_micro_wu(Me_micro),
	M_micro_vu(Me_micro), M_micro_wv(Me_micro), M_micro_uw(Me_micro);

	libMesh::DenseSubMatrix<libMesh::Number>
	M_BIG_uu(Me_BIG), M_BIG_vv(Me_BIG), M_BIG_ww(Me_BIG);

	libMesh::DenseSubMatrix<libMesh::Number>
	M_BIG_uv(Me_BIG), M_BIG_vw(Me_BIG), M_BIG_wu(Me_BIG),
	M_BIG_vu(Me_BIG), M_BIG_wv(Me_BIG), M_BIG_uw(Me_BIG);

	//    If all elements are of the same type, do the index "extraction",
	// the matrices resizes and repositions here
	if(bSameElemsType)
	{
		const libMesh::Elem* elem_restrict = mesh_restrict.elem(0);
		dof_map_restrict.dof_indices(elem_restrict, dof_indices_restrict);
		dof_map_restrict.dof_indices(elem_restrict, dof_indices_u_restrict, u_var_restrict);
		dof_map_restrict.dof_indices(elem_restrict, dof_indices_v_restrict, v_var_restrict);
		dof_map_restrict.dof_indices(elem_restrict, dof_indices_w_restrict, w_var_restrict);
		n_dofs_restrict   = dof_indices_restrict.size();
		n_dofs_u_restrict = dof_indices_u_restrict.size();
		n_dofs_v_restrict = dof_indices_v_restrict.size();
		n_dofs_w_restrict = dof_indices_w_restrict.size();

		const libMesh::Elem* elem_BIG = mesh_BIG.elem(0);
		dof_map_BIG.dof_indices(elem_BIG, dof_indices_BIG);
		dof_map_BIG.dof_indices(elem_BIG, dof_indices_u_BIG, u_var_BIG);
		dof_map_BIG.dof_indices(elem_BIG, dof_indices_v_BIG, v_var_BIG);
		dof_map_BIG.dof_indices(elem_BIG, dof_indices_w_BIG, w_var_BIG);
		n_dofs_BIG   = dof_indices_BIG.size();
		n_dofs_u_BIG = dof_indices_u_BIG.size();
		n_dofs_v_BIG = dof_indices_v_BIG.size();
		n_dofs_w_BIG = dof_indices_w_BIG.size();

		const libMesh::Elem* elem_micro = mesh_micro.elem(0);
		dof_map_micro.dof_indices(elem_micro, dof_indices_micro);
		dof_map_micro.dof_indices(elem_micro, dof_indices_u_micro, u_var_micro);
		dof_map_micro.dof_indices(elem_micro, dof_indices_v_micro, v_var_micro);
		dof_map_micro.dof_indices(elem_micro, dof_indices_w_micro, w_var_micro);
		n_dofs_micro   = dof_indices_micro.size();
		n_dofs_u_micro = dof_indices_u_micro.size();
		n_dofs_v_micro = dof_indices_v_micro.size();
		n_dofs_w_micro = dof_indices_w_micro.size();

		// Resize matrices
		Me_micro.resize (n_dofs_restrict, n_dofs_micro);

		M_micro_uu.reposition (u_var_restrict*n_dofs_u_restrict, u_var_micro*n_dofs_u_micro, n_dofs_u_restrict, n_dofs_u_micro);
//		M_micro_uv.reposition (u_var_restrict*n_dofs_u_restrict, v_var_micro*n_dofs_u_micro, n_dofs_u_restrict, n_dofs_v_micro);
//		M_micro_uw.reposition (u_var_restrict*n_dofs_u_restrict, w_var_micro*n_dofs_u_micro, n_dofs_u_restrict, n_dofs_w_micro);

//		M_micro_vu.reposition (v_var_restrict*n_dofs_u_restrict, u_var_micro*n_dofs_u_micro, n_dofs_v_restrict, n_dofs_u_micro);
		M_micro_vv.reposition (v_var_restrict*n_dofs_u_restrict, v_var_micro*n_dofs_u_micro, n_dofs_v_restrict, n_dofs_v_micro);
//		M_micro_vw.reposition (v_var_restrict*n_dofs_u_restrict, w_var_micro*n_dofs_u_micro, n_dofs_v_restrict, n_dofs_w_micro);

//		M_micro_wu.reposition (w_var_restrict*n_dofs_u_restrict, u_var_micro*n_dofs_u_micro, n_dofs_w_restrict, n_dofs_u_micro);
//		M_micro_wv.reposition (w_var_restrict*n_dofs_u_restrict, v_var_micro*n_dofs_u_micro, n_dofs_w_restrict, n_dofs_v_micro);
		M_micro_ww.reposition (w_var_restrict*n_dofs_u_restrict, w_var_micro*n_dofs_u_micro, n_dofs_w_restrict, n_dofs_w_micro);

		Me_BIG.resize (n_dofs_restrict, n_dofs_BIG);

		M_BIG_uu.reposition (u_var_restrict*n_dofs_u_restrict, u_var_BIG*n_dofs_u_BIG, n_dofs_u_restrict, n_dofs_u_BIG);
//		M_BIG_uv.reposition (u_var_restrict*n_dofs_u_restrict, v_var_BIG*n_dofs_u_BIG, n_dofs_u_restrict, n_dofs_v_BIG);
//		M_BIG_uw.reposition (u_var_restrict*n_dofs_u_restrict, w_var_BIG*n_dofs_u_BIG, n_dofs_u_restrict, n_dofs_w_BIG);
//
//		M_BIG_vu.reposition (v_var_restrict*n_dofs_u_restrict, u_var_BIG*n_dofs_u_BIG, n_dofs_v_restrict, n_dofs_u_BIG);
		M_BIG_vv.reposition (v_var_restrict*n_dofs_u_restrict, v_var_BIG*n_dofs_u_BIG, n_dofs_v_restrict, n_dofs_v_BIG);
//		M_BIG_vw.reposition (v_var_restrict*n_dofs_u_restrict, w_var_BIG*n_dofs_u_BIG, n_dofs_v_restrict, n_dofs_w_BIG);
//
//		M_BIG_wu.reposition (w_var_restrict*n_dofs_u_restrict, u_var_BIG*n_dofs_u_BIG, n_dofs_w_restrict, n_dofs_u_BIG);
//		M_BIG_wv.reposition (w_var_restrict*n_dofs_u_restrict, v_var_BIG*n_dofs_u_BIG, n_dofs_w_restrict, n_dofs_v_BIG);
		M_BIG_ww.reposition (w_var_restrict*n_dofs_u_restrict, w_var_BIG*n_dofs_u_BIG, n_dofs_w_restrict, n_dofs_w_BIG);
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
	couplingMatrix_restrict_micro.init(		restrict_M, micro_N,
											restrict_M_local, micro_N_local,
											micro_N_local,micro_N - micro_N_local);
	couplingMatrix_restrict_BIG.init(		restrict_M, BIG_N,
											restrict_M_local, BIG_N_local,
											BIG_N_local,BIG_N - BIG_N_local);

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
		dof_map_restrict.dof_indices(elem_restrict, dof_indices_u_restrict, u_var_restrict);
		dof_map_restrict.dof_indices(elem_restrict, dof_indices_v_restrict, v_var_restrict);
		dof_map_restrict.dof_indices(elem_restrict, dof_indices_w_restrict, w_var_restrict);

		dof_map_micro.dof_indices(elem_micro, dof_indices_micro);
		dof_map_micro.dof_indices(elem_micro, dof_indices_u_micro, u_var_micro);
		dof_map_micro.dof_indices(elem_micro, dof_indices_v_micro, v_var_micro);
		dof_map_micro.dof_indices(elem_micro, dof_indices_w_micro, w_var_micro);

		if(using_same_mesh_restrict_A)
		{
			dof_indices_BIG = dof_indices_restrict;
			dof_indices_u_BIG = dof_indices_u_restrict;
			dof_indices_v_BIG = dof_indices_v_restrict;
			dof_indices_w_BIG = dof_indices_w_restrict;
		}
		else
		{
			dof_map_BIG.dof_indices(elem_BIG, dof_indices_BIG);
			dof_map_BIG.dof_indices(elem_BIG, dof_indices_u_BIG, u_var_BIG);
			dof_map_BIG.dof_indices(elem_BIG, dof_indices_v_BIG, v_var_BIG);
			dof_map_BIG.dof_indices(elem_BIG, dof_indices_w_BIG, w_var_BIG);
		}

		// Resize dense matrix, if needed
		if(!bSameElemsType)
		{
			n_dofs_restrict   = dof_indices_restrict.size();
			n_dofs_u_restrict = dof_indices_u_restrict.size();
			n_dofs_v_restrict = dof_indices_v_restrict.size();
			n_dofs_w_restrict = dof_indices_w_restrict.size();

			n_dofs_BIG   = dof_indices_BIG.size();
			n_dofs_u_BIG = dof_indices_u_BIG.size();
			n_dofs_v_BIG = dof_indices_v_BIG.size();
			n_dofs_w_BIG = dof_indices_w_BIG.size();

			n_dofs_micro   = dof_indices_micro.size();
			n_dofs_u_micro = dof_indices_u_micro.size();
			n_dofs_v_micro = dof_indices_v_micro.size();
			n_dofs_w_micro = dof_indices_w_micro.size();

			// Resize matrices
			Me_micro.resize (n_dofs_restrict, n_dofs_micro);

			M_micro_uu.reposition (u_var_restrict*n_dofs_u_restrict, u_var_micro*n_dofs_u_micro, n_dofs_u_restrict, n_dofs_u_micro);
//			M_micro_uv.reposition (u_var_restrict*n_dofs_u_restrict, v_var_micro*n_dofs_u_micro, n_dofs_u_restrict, n_dofs_v_micro);
//			M_micro_uw.reposition (u_var_restrict*n_dofs_u_restrict, w_var_micro*n_dofs_u_micro, n_dofs_u_restrict, n_dofs_w_micro);
//
//			M_micro_vu.reposition (v_var_restrict*n_dofs_u_restrict, u_var_micro*n_dofs_u_micro, n_dofs_v_restrict, n_dofs_u_micro);
			M_micro_vv.reposition (v_var_restrict*n_dofs_u_restrict, v_var_micro*n_dofs_u_micro, n_dofs_v_restrict, n_dofs_v_micro);
//			M_micro_vw.reposition (v_var_restrict*n_dofs_u_restrict, w_var_micro*n_dofs_u_micro, n_dofs_v_restrict, n_dofs_w_micro);
//
//			M_micro_wu.reposition (w_var_restrict*n_dofs_u_restrict, u_var_micro*n_dofs_u_micro, n_dofs_w_restrict, n_dofs_u_micro);
//			M_micro_wv.reposition (w_var_restrict*n_dofs_u_restrict, v_var_micro*n_dofs_u_micro, n_dofs_w_restrict, n_dofs_v_micro);
			M_micro_ww.reposition (w_var_restrict*n_dofs_u_restrict, w_var_micro*n_dofs_u_micro, n_dofs_w_restrict, n_dofs_w_micro);

			Me_BIG.resize (n_dofs_restrict, n_dofs_BIG);

			M_BIG_uu.reposition (u_var_restrict*n_dofs_u_restrict, u_var_BIG*n_dofs_u_BIG, n_dofs_u_restrict, n_dofs_u_BIG);
//			M_BIG_uv.reposition (u_var_restrict*n_dofs_u_restrict, v_var_BIG*n_dofs_u_BIG, n_dofs_u_restrict, n_dofs_v_BIG);
//			M_BIG_uw.reposition (u_var_restrict*n_dofs_u_restrict, w_var_BIG*n_dofs_u_BIG, n_dofs_u_restrict, n_dofs_w_BIG);
//
//			M_BIG_vu.reposition (v_var_restrict*n_dofs_u_restrict, u_var_BIG*n_dofs_u_BIG, n_dofs_v_restrict, n_dofs_u_BIG);
			M_BIG_vv.reposition (v_var_restrict*n_dofs_u_restrict, v_var_BIG*n_dofs_u_BIG, n_dofs_v_restrict, n_dofs_v_BIG);
//			M_BIG_vw.reposition (v_var_restrict*n_dofs_u_restrict, w_var_BIG*n_dofs_u_BIG, n_dofs_v_restrict, n_dofs_w_BIG);
//
//			M_BIG_wu.reposition (w_var_restrict*n_dofs_u_restrict, u_var_BIG*n_dofs_u_BIG, n_dofs_w_restrict, n_dofs_u_BIG);
//			M_BIG_wv.reposition (w_var_restrict*n_dofs_u_restrict, v_var_BIG*n_dofs_u_BIG, n_dofs_w_restrict, n_dofs_v_BIG);
			M_BIG_ww.reposition (w_var_restrict*n_dofs_u_restrict, w_var_BIG*n_dofs_u_BIG, n_dofs_w_restrict, n_dofs_w_BIG);
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
				L2_Coupling(M_micro_uu,qp,phi_restrict,phi_micro,
							n_dofs_u_restrict,n_dofs_u_micro,JxW,1);

//				L2_Coupling(M_micro_uv,qp,phi_restrict,phi_micro,
//							n_dofs_u_restrict,n_dofs_v_micro,JxW,1);
//
//				L2_Coupling(M_micro_uw,qp,phi_restrict,phi_micro,
//							n_dofs_u_restrict,n_dofs_w_micro,JxW,1);
//
//				L2_Coupling(M_micro_vu,qp,phi_restrict,phi_micro,
//							n_dofs_v_restrict,n_dofs_u_micro,JxW,1);

				L2_Coupling(M_micro_vv,qp,phi_restrict,phi_micro,
							n_dofs_v_restrict,n_dofs_v_micro,JxW,1);

//				L2_Coupling(M_micro_vw,qp,phi_restrict,phi_micro,
//							n_dofs_v_restrict,n_dofs_w_micro,JxW,1);
//
//				L2_Coupling(M_micro_wu,qp,phi_restrict,phi_micro,
//							n_dofs_w_restrict,n_dofs_u_micro,JxW,1);
//
//				L2_Coupling(M_micro_wv,qp,phi_restrict,phi_micro,
//							n_dofs_w_restrict,n_dofs_v_micro,JxW,1);

				L2_Coupling(M_micro_ww,qp,phi_restrict,phi_micro,
							n_dofs_w_restrict,n_dofs_w_micro,JxW,1);

				// Restrict -> restrict coupling
				L2_Coupling(M_BIG_uu,qp,phi_restrict,phi_restrict,
							n_dofs_u_restrict,n_dofs_u_restrict,JxW,1);

//				L2_Coupling(M_BIG_uv,qp,phi_restrict,phi_restrict,
//							n_dofs_u_restrict,n_dofs_v_restrict,JxW,1);
//
//				L2_Coupling(M_BIG_uw,qp,phi_restrict,phi_restrict,
//							n_dofs_u_restrict,n_dofs_w_restrict,JxW,1);
//
//				L2_Coupling(M_BIG_vu,qp,phi_restrict,phi_restrict,
//							n_dofs_v_restrict,n_dofs_u_restrict,JxW,1);

				L2_Coupling(M_BIG_vv,qp,phi_restrict,phi_restrict,
							n_dofs_v_restrict,n_dofs_v_restrict,JxW,1);

//				L2_Coupling(M_BIG_vw,qp,phi_restrict,phi_restrict,
//							n_dofs_v_restrict,n_dofs_w_restrict,JxW,1);
//
//				L2_Coupling(M_BIG_wu,qp,phi_restrict,phi_restrict,
//							n_dofs_w_restrict,n_dofs_u_restrict,JxW,1);
//
//				L2_Coupling(M_BIG_wv,qp,phi_restrict,phi_restrict,
//							n_dofs_w_restrict,n_dofs_v_restrict,JxW,1);

				L2_Coupling(M_BIG_ww,qp,phi_restrict,phi_restrict,
							n_dofs_w_restrict,n_dofs_w_restrict,JxW,1);
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

void carl::coupled_system::set_LATIN_solver(const std::string micro_name, const std::string type_name)
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

	libMesh::PetscMatrix<libMesh::Number>& C_RA = * m_couplingMatrixMap_restrict_BIG[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& C_RB = * m_couplingMatrixMap_restrict_micro[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& C_RR = * m_couplingMatrixMap_restrict_restrict[micro_name];

	// Get the vectors
	libMesh::PetscVector<libMesh::Number>& F_A = libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>& >(* Sys_BIG.rhs);
	libMesh::PetscVector<libMesh::Number>& F_B = libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>& >(* Sys_micro.rhs);

	// Set the solver parameters
	m_LATIN_solver.set_params(2.5,2.5,2.5,2.5);

	// Set the solver matrices
	m_LATIN_solver.set_matrices(M_A,M_B,C_RA,C_RB,C_RR);

	// Set the solver matrices
	m_LATIN_solver.set_forces(F_A,F_B);
};

void carl::coupled_system::set_LATIN_solver(const std::string micro_name, const std::string type_name,
												void fptr_BIG(		libMesh::EquationSystems& es,
																	const std::string& name, weight_parameter_function& weight_mask),
												void fptr_micro(	libMesh::EquationSystems& es,
																	const std::string& name, weight_parameter_function& weight_mask))
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

	libMesh::PetscMatrix<libMesh::Number>& C_RA = * m_couplingMatrixMap_restrict_BIG[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& C_RB = * m_couplingMatrixMap_restrict_micro[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& C_RR = * m_couplingMatrixMap_restrict_restrict[micro_name];

	// Get the vectors
	libMesh::PetscVector<libMesh::Number>& F_A = libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>& >(* Sys_BIG.rhs);
	libMesh::PetscVector<libMesh::Number>& F_B = libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>& >(* Sys_micro.rhs);

	// Set the solver parameters
	m_LATIN_solver.set_params(2.5,2.5,2.5,2.5);

	// Set the solver matrices
	m_LATIN_solver.set_matrices(M_A,M_B,C_RA,C_RB,C_RR);

	// Set the solver matrices
	m_LATIN_solver.set_forces(F_A,F_B);
};

void carl::coupled_system::solve_LATIN(const std::string micro_name, const std::string type_name)
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
	*(Sys_micro.solution) = sol_micro;
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
