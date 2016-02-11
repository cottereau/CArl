#include "coupled_system.h"

void carl::coupled_system::assemble_coupling_matrices(	const std::string BIG_name,
														const std::string micro_name,
														const std::string inter_name,
														const std::string restrict_name,
														std::unordered_map<int,int>& equivalence_table_restrict_BIG,
														std::vector<std::pair<int,int> >& intersection_table_restrict_micro,
														std::unordered_multimap<int,int>& intersection_table_inter,
														double coupling_const,
														bool using_same_mesh_restrict_A,
														bool bSameElemsType)
{
	// Addresses to the eq. systems
	libMesh::EquationSystems& restrict_eq_system = * m_restrict_EquationSystemMap[restrict_name];
	libMesh::EquationSystems& BIG_eq_system = * m_BIG_EquationSystem.second;
	libMesh::EquationSystems& micro_eq_system = * m_micro_EquationSystemMap[micro_name];
	libMesh::EquationSystems& inter_eq_system = * m_inter_EquationSystemMap[inter_name];

	// First, test if all the systems have an acoustic model and variable set
	homemade_assert_msg(micro_eq_system.has_system("VolTest"),
				" Micro equation systems missing \"VolTest\" system!");
	homemade_assert_msg(BIG_eq_system.has_system("VolTest"),
				" Macro equation systems missing \"VolTest\" system!");
	homemade_assert_msg(inter_eq_system.has_system("VolTest"),
				" Intersection equation systems missing \"VolTest\" system!");
	homemade_assert_msg(restrict_eq_system.has_system("VolTest"),
				" Restricted equation systems missing \"VolTest\" system!");

	// Systems and vars
	libMesh::LinearImplicitSystem& volume_restrict_system
		= restrict_eq_system.get_system<libMesh::LinearImplicitSystem>("VolTest");
	libMesh::LinearImplicitSystem& volume_BIG_system
		= BIG_eq_system.get_system<libMesh::LinearImplicitSystem>("VolTest");
	libMesh::LinearImplicitSystem& volume_micro_system
		= micro_eq_system.get_system<libMesh::LinearImplicitSystem>("VolTest");
	libMesh::LinearImplicitSystem& volume_inter_system
		= inter_eq_system.get_system<libMesh::LinearImplicitSystem>("VolTest");

	// Addresses to the meshes
	const libMesh::MeshBase& mesh_restrict = restrict_eq_system.get_mesh();
	const libMesh::MeshBase& mesh_BIG = BIG_eq_system.get_mesh();
	const libMesh::MeshBase& mesh_micro = micro_eq_system.get_mesh();
	const libMesh::MeshBase& mesh_inter = inter_eq_system.get_mesh();

	// Variable indices and dimensions
	const unsigned int dim_restrict = mesh_restrict.mesh_dimension();
	const unsigned int dim_BIG = mesh_BIG.mesh_dimension();
	const unsigned int dim_micro = mesh_micro.mesh_dimension();
	const unsigned int dim_inter = mesh_inter.mesh_dimension();

	const unsigned int vol_test_restrict = volume_restrict_system.variable_number("SillyVar");
	const unsigned int vol_test_BIG = volume_BIG_system.variable_number ("SillyVar");
	const unsigned int vol_test_micro = volume_micro_system.variable_number ("SillyVar");
	const unsigned int vol_test_inter = volume_inter_system.variable_number ("SillyVar");

	// Addresses to the DoF maps
	const libMesh::DofMap& dof_map_restrict = volume_restrict_system.get_dof_map();
	const libMesh::DofMap& dof_map_BIG = volume_BIG_system.get_dof_map();
	const libMesh::DofMap& dof_map_micro = volume_micro_system.get_dof_map();
	const libMesh::DofMap& dof_map_inter = volume_inter_system.get_dof_map();

	// Finite elements (FE) declaration
	libMesh::FEType fe_type_restrict = dof_map_restrict.variable_type(vol_test_restrict);
	libMesh::FEType fe_type_BIG = dof_map_BIG.variable_type(vol_test_BIG);
	libMesh::FEType fe_type_micro = dof_map_micro.variable_type(vol_test_micro);
	libMesh::FEType fe_type_inter = dof_map_inter.variable_type(vol_test_inter);

	// Test if we are using linear elements
	homemade_assert_msg(fe_type_restrict.order == libMesh::FIRST,
										" Restrict system is not linear!");
	homemade_assert_msg(fe_type_BIG.order == libMesh::FIRST,
										" Macro system is not linear!");
	homemade_assert_msg(fe_type_micro.order == libMesh::FIRST,
										" Micro system is not linear!");
	homemade_assert_msg(fe_type_inter.order == libMesh::FIRST,
										" Intersection system is not linear!");

	// Set up FE bases and quadratures
	libMesh::UniquePtr<libMesh::FEBase> fe_restrict (libMesh::FEBase::build(dim_restrict, fe_type_restrict));
	libMesh::UniquePtr<libMesh::FEBase> fe_BIG (libMesh::FEBase::build(dim_BIG, fe_type_BIG));
	libMesh::UniquePtr<libMesh::FEBase> fe_micro (libMesh::FEBase::build(dim_micro, fe_type_micro));
	libMesh::UniquePtr<libMesh::FEBase> fe_inter (libMesh::FEBase::build(dim_inter, fe_type_inter));

	libMesh::QGauss qrule_restrict(dim_restrict, fe_type_restrict.default_quadrature_order());
	fe_restrict->attach_quadrature_rule (&qrule_restrict);

	libMesh::QGauss qrule_inter(dim_inter, fe_type_inter.default_quadrature_order());
	fe_inter->attach_quadrature_rule (&qrule_inter);

	// Vector that will keep the quadrature points
	const std::vector<libMesh::Point>& quad_points_inter = fe_inter->get_xyz();
	std::vector<libMesh::Point> quad_points_reference;

	// Jacobians
	const std::vector<libMesh::Real>& JxW = fe_inter->get_JxW();

	// Shape functions
	const std::vector<std::vector<libMesh::Real> >& phi_inter = fe_inter->get_phi();
	std::vector<std::vector<libMesh::Real> > corrected_phi_BIG;
	std::vector<std::vector<libMesh::Real> > corrected_phi_micro;
	std::vector<std::vector<libMesh::Real> > corrected_phi_restrict;
	unsigned int n_quadrature_pts = 0;

	// Addresses to the matrices
	libMesh::PetscMatrix<libMesh::Number>& couplingMatrix_restrict_BIG = * m_couplingMatrixMap_restrict_BIG[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& couplingMatrix_restrict_micro = * m_couplingMatrixMap_restrict_micro[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& couplingMatrix_restrict_restrict = * m_couplingMatrixMap_restrict_restrict[micro_name];

	// Dof vectors and ranges
	std::vector<libMesh::dof_id_type> dof_indices_restrict;
	std::vector<libMesh::dof_id_type> dof_indices_BIG;
	std::vector<libMesh::dof_id_type> dof_indices_micro;
	std::vector<libMesh::dof_id_type> dof_indices_inter;

	// Number of DoF's for each system and each variable
	unsigned int n_dofs_restrict;
	unsigned int n_dofs_micro;
	unsigned int n_dofs_BIG;
	unsigned int n_dofs_inter;

	// Local matrix
	libMesh::DenseMatrix<libMesh::Number> Me_micro;
	libMesh::DenseMatrix<libMesh::Number> Me_BIG;
	libMesh::DenseMatrix<libMesh::Number> Me_restrict;

	//    If all elements are of the same type, do the index "extraction",
	// the matrices resizes and repositions here
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

		const libMesh::Elem* elem_inter = mesh_inter.elem(0);
		dof_map_inter.dof_indices(elem_inter, dof_indices_inter);
		n_dofs_inter   = dof_indices_inter.size();

		// Resize matrices
		Me_BIG.resize (n_dofs_restrict, n_dofs_BIG);
		Me_micro.resize (n_dofs_restrict, n_dofs_micro);
		Me_restrict.resize (n_dofs_restrict, n_dofs_restrict);

		// Set up corrected shape vectors
		fe_inter->reinit(elem_inter);
		n_quadrature_pts = fe_inter->n_quadrature_points();

		corrected_phi_BIG.resize(n_dofs_BIG,std::vector<libMesh::Real>(n_quadrature_pts,0));
		corrected_phi_restrict.resize(n_dofs_restrict,std::vector<libMesh::Real>(n_quadrature_pts,0));
		corrected_phi_micro.resize(n_dofs_micro,std::vector<libMesh::Real>(n_quadrature_pts,0));
	}

	// Vectors containing the lambda weights
	std::vector<std::vector<libMesh::Real> > 	lambda_weight_restrict(
													n_dofs_restrict,
													std::vector<libMesh::Real> (n_dofs_inter,0));
	std::vector<std::vector<libMesh::Real> > 	lambda_weight_BIG(
													n_dofs_BIG,
													std::vector<libMesh::Real> (n_dofs_inter,0));
	std::vector<std::vector<libMesh::Real> > 	lambda_weight_micro(
													n_dofs_micro,
													std::vector<libMesh::Real> (n_dofs_inter,0));

	// Initialize global matrix
	const unsigned int restrict_M = dof_map_restrict.n_dofs();

	const unsigned int BIG_N = dof_map_BIG.n_dofs();
	const unsigned int micro_N = dof_map_micro.n_dofs();

	const unsigned int restrict_M_local = dof_map_restrict.n_local_dofs();

	const unsigned int BIG_N_local = dof_map_BIG.n_local_dofs();
	const unsigned int micro_N_local = dof_map_micro.n_local_dofs();

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
			n_dofs_restrict   = dof_indices_restrict.size();
			n_dofs_BIG   = dof_indices_BIG.size();
			n_dofs_micro   = dof_indices_micro.size();

			// Resize matrices
			Me_micro.resize (n_dofs_restrict, n_dofs_micro);
			Me_BIG.resize (n_dofs_restrict, n_dofs_BIG);
			Me_BIG.resize (n_dofs_restrict, n_dofs_restrict);

			// Set up corrected shape vectors
			const libMesh::Elem* elem_inter = mesh_inter.elem(0);
			fe_inter->reinit(elem_inter);
			n_quadrature_pts = fe_inter->n_quadrature_points();

			corrected_phi_BIG.resize(n_dofs_BIG,std::vector<libMesh::Real>(n_quadrature_pts,0));
			corrected_phi_restrict.resize(n_dofs_restrict,std::vector<libMesh::Real>(n_quadrature_pts,0));
			corrected_phi_micro.resize(n_dofs_micro,std::vector<libMesh::Real>(n_quadrature_pts,0));

		}

		Me_micro.zero();
		Me_BIG.zero();
		Me_restrict.zero();

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

			get_lambdas(dim_BIG, fe_type_BIG, elem_BIG,
							quad_points_inter,
							quad_points_reference,
							lambda_weight_BIG);

			get_lambdas(dim_micro, fe_type_micro, elem_micro,
							quad_points_inter,
							quad_points_reference,
							lambda_weight_micro);

			get_lambdas(dim_restrict, fe_type_restrict, elem_restrict,
							quad_points_inter,
							quad_points_reference,
							lambda_weight_restrict);

			set_corrected_shapes(lambda_weight_BIG,phi_inter,corrected_phi_BIG);
			set_corrected_shapes(lambda_weight_micro,phi_inter,corrected_phi_micro);
			set_corrected_shapes(lambda_weight_restrict,phi_inter,corrected_phi_restrict);

			// For each quadrature point determinate the sub-matrices elements
			for (unsigned int qp=0; qp<qrule_inter.n_points(); qp++)
			{
				// Restrict -> micro coupling
				L2_Coupling(Me_micro,qp,corrected_phi_restrict,corrected_phi_micro,
							n_dofs_restrict,n_dofs_micro,JxW,coupling_const);

				// Restrict -> BIG coupling
				L2_Coupling(Me_BIG,qp,corrected_phi_restrict,corrected_phi_BIG,
							n_dofs_restrict,n_dofs_BIG,JxW,coupling_const);

				// Restrict -> restrict coupling
				L2_Coupling(Me_restrict,qp,corrected_phi_restrict,corrected_phi_restrict,
							n_dofs_restrict,n_dofs_restrict,JxW,coupling_const);
			}
		}

		couplingMatrix_restrict_micro.add_matrix(Me_micro,dof_indices_restrict,dof_indices_micro);
		couplingMatrix_restrict_BIG.add_matrix(Me_BIG,dof_indices_restrict,dof_indices_BIG);
		couplingMatrix_restrict_restrict.add_matrix(Me_restrict,dof_indices_restrict,dof_indices_restrict);
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
															double coupling_const,
															bool using_same_mesh_restrict_A,
															bool bSameElemsType)
{
	// Addresses to the eq. systems
	libMesh::EquationSystems& restrict_eq_system = * m_restrict_EquationSystemMap[restrict_name];
	libMesh::EquationSystems& BIG_eq_system = * m_BIG_EquationSystem.second;
	libMesh::EquationSystems& micro_eq_system = * m_micro_EquationSystemMap[micro_name];
	libMesh::EquationSystems& inter_eq_system = * m_inter_EquationSystemMap[inter_name];

	// First, test if all the systems have an elasticity model and variable set
	homemade_assert_msg(micro_eq_system.has_system("Elasticity"),
				" Micro equation systems missing \"Elasticity\" system!");
	homemade_assert_msg(BIG_eq_system.has_system("Elasticity"),
				" Macro equation systems missing \"Elasticity\" system!");
	homemade_assert_msg(inter_eq_system.has_system("Elasticity"),
				" Intersection equation systems missing \"Elasticity\" system!");
	homemade_assert_msg(restrict_eq_system.has_system("Elasticity"),
				" Restricted equation systems missing \"Elasticity\" system!");

	// Systems and vars
	libMesh::LinearImplicitSystem& volume_restrict_system
		= restrict_eq_system.get_system<libMesh::LinearImplicitSystem>("Elasticity");
	libMesh::LinearImplicitSystem& volume_BIG_system
		= BIG_eq_system.get_system<libMesh::LinearImplicitSystem>("Elasticity");
	libMesh::LinearImplicitSystem& volume_micro_system
		= micro_eq_system.get_system<libMesh::LinearImplicitSystem>("Elasticity");
	libMesh::LinearImplicitSystem& volume_inter_system
		= inter_eq_system.get_system<libMesh::LinearImplicitSystem>("Elasticity");

	// Addresses to the meshes
	const libMesh::MeshBase& mesh_restrict = restrict_eq_system.get_mesh();
	const libMesh::MeshBase& mesh_BIG = BIG_eq_system.get_mesh();
	const libMesh::MeshBase& mesh_micro = micro_eq_system.get_mesh();
	const libMesh::MeshBase& mesh_inter = inter_eq_system.get_mesh();

	// Variable indices and dimensions
	const unsigned int dim_restrict = mesh_restrict.mesh_dimension();
	const unsigned int dim_BIG = mesh_BIG.mesh_dimension();
	const unsigned int dim_micro = mesh_micro.mesh_dimension();
	const unsigned int dim_inter = mesh_inter.mesh_dimension();

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

	// Addresses to the  DoF maps
	const libMesh::DofMap& dof_map_restrict = volume_restrict_system.get_dof_map();
	const libMesh::DofMap& dof_map_BIG = volume_BIG_system.get_dof_map();
	const libMesh::DofMap& dof_map_micro = volume_micro_system.get_dof_map();
	const libMesh::DofMap& dof_map_inter = volume_inter_system.get_dof_map();

	// Finite elements (FE) declaration
	libMesh::FEType fe_type_restrict = dof_map_restrict.variable_type(u_var_restrict);
	libMesh::FEType fe_type_BIG = dof_map_BIG.variable_type(u_var_BIG);
	libMesh::FEType fe_type_micro = dof_map_micro.variable_type(u_var_micro);
	libMesh::FEType fe_type_inter = dof_map_inter.variable_type(u_var_inter);

	// Test if we are using linear elements
	homemade_assert_msg(fe_type_restrict.order == libMesh::FIRST,
										" Restrict system is not linear!");
	homemade_assert_msg(fe_type_BIG.order == libMesh::FIRST,
										" Macro system is not linear!");
	homemade_assert_msg(fe_type_micro.order == libMesh::FIRST,
										" Micro system is not linear!");
	homemade_assert_msg(fe_type_inter.order == libMesh::FIRST,
										" Intersection system is not linear!");

	// Set up FE bases and quadratures
	libMesh::UniquePtr<libMesh::FEBase> fe_restrict (libMesh::FEBase::build(dim_restrict, fe_type_restrict));
	libMesh::UniquePtr<libMesh::FEBase> fe_BIG (libMesh::FEBase::build(dim_BIG, fe_type_BIG));
	libMesh::UniquePtr<libMesh::FEBase> fe_micro (libMesh::FEBase::build(dim_micro, fe_type_micro));
	libMesh::UniquePtr<libMesh::FEBase> fe_inter (libMesh::FEBase::build(dim_inter, fe_type_inter));

	libMesh::QGauss qrule_BIG(dim_BIG, fe_type_BIG.default_quadrature_order());
	fe_BIG->attach_quadrature_rule (&qrule_BIG);

	libMesh::QGauss qrule_restrict(dim_restrict, fe_type_restrict.default_quadrature_order());
	fe_restrict->attach_quadrature_rule (&qrule_restrict);

	libMesh::QGauss qrule_micro(dim_micro, fe_type_micro.default_quadrature_order());
	fe_micro->attach_quadrature_rule (&qrule_micro);

	libMesh::QGauss qrule_inter(dim_inter, fe_type_inter.default_quadrature_order());
	fe_inter->attach_quadrature_rule (&qrule_inter);

	// Vector that will keep the quadrature points
	const std::vector<libMesh::Point>& quad_points_inter = fe_inter->get_xyz();
	std::vector<libMesh::Point> quad_points_reference;

	// Jacobians
	const std::vector<libMesh::Real>& JxW = fe_inter->get_JxW();

	// Shape functions
	const std::vector<std::vector<libMesh::Real> >& phi_inter = fe_inter->get_phi();
	std::vector<std::vector<libMesh::Real> > corrected_phi_BIG;
	std::vector<std::vector<libMesh::Real> > corrected_phi_micro;
	std::vector<std::vector<libMesh::Real> > corrected_phi_restrict;
	unsigned int n_quadrature_pts = 0;

	// Addresses to the matrices
	libMesh::PetscMatrix<libMesh::Number>& couplingMatrix_restrict_BIG = * m_couplingMatrixMap_restrict_BIG[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& couplingMatrix_restrict_micro = * m_couplingMatrixMap_restrict_micro[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& couplingMatrix_restrict_restrict = * m_couplingMatrixMap_restrict_restrict[micro_name];

	// DoF vectors and ranges
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

	std::vector<libMesh::dof_id_type> dof_indices_inter;
	std::vector<libMesh::dof_id_type> dof_indices_u_inter;
	std::vector<libMesh::dof_id_type> dof_indices_v_inter;
	std::vector<libMesh::dof_id_type> dof_indices_w_inter;

	// Number of DoF's for each system and each variable
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

	unsigned int n_dofs_inter;
	unsigned int n_dofs_u_inter;
	unsigned int n_dofs_v_inter;
	unsigned int n_dofs_w_inter;

	// Local matrix
	libMesh::DenseMatrix<libMesh::Number> Me_micro;
	libMesh::DenseMatrix<libMesh::Number> Me_BIG;
	libMesh::DenseMatrix<libMesh::Number> Me_restrict;

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

	libMesh::DenseSubMatrix<libMesh::Number>
	M_restrict_uu(Me_restrict), M_restrict_vv(Me_restrict), M_restrict_ww(Me_restrict);

	libMesh::DenseSubMatrix<libMesh::Number>
	M_restrict_uv(Me_restrict), M_restrict_vw(Me_restrict), M_restrict_wu(Me_restrict),
	M_restrict_vu(Me_restrict), M_restrict_wv(Me_restrict), M_restrict_uw(Me_restrict);

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

		const libMesh::Elem* elem_inter = mesh_inter.elem(0);
		dof_map_inter.dof_indices(elem_inter, dof_indices_inter);
		dof_map_inter.dof_indices(elem_inter, dof_indices_u_inter, u_var_inter);
		dof_map_inter.dof_indices(elem_inter, dof_indices_v_inter, v_var_inter);
		dof_map_inter.dof_indices(elem_inter, dof_indices_w_inter, w_var_inter);
		n_dofs_inter   = dof_indices_inter.size();
		n_dofs_u_inter = dof_indices_u_inter.size();
		n_dofs_v_inter = dof_indices_v_inter.size();
		n_dofs_w_inter = dof_indices_w_inter.size();

		// Resize matrices
		Me_micro.resize (n_dofs_restrict, n_dofs_micro);

		M_micro_uu.reposition (u_var_restrict*n_dofs_u_restrict, u_var_micro*n_dofs_u_micro, n_dofs_u_restrict, n_dofs_u_micro);
		M_micro_uv.reposition (u_var_restrict*n_dofs_u_restrict, v_var_micro*n_dofs_u_micro, n_dofs_u_restrict, n_dofs_v_micro);
		M_micro_uw.reposition (u_var_restrict*n_dofs_u_restrict, w_var_micro*n_dofs_u_micro, n_dofs_u_restrict, n_dofs_w_micro);

		M_micro_vu.reposition (v_var_restrict*n_dofs_u_restrict, u_var_micro*n_dofs_u_micro, n_dofs_v_restrict, n_dofs_u_micro);
		M_micro_vv.reposition (v_var_restrict*n_dofs_u_restrict, v_var_micro*n_dofs_u_micro, n_dofs_v_restrict, n_dofs_v_micro);
		M_micro_vw.reposition (v_var_restrict*n_dofs_u_restrict, w_var_micro*n_dofs_u_micro, n_dofs_v_restrict, n_dofs_w_micro);

		M_micro_wu.reposition (w_var_restrict*n_dofs_u_restrict, u_var_micro*n_dofs_u_micro, n_dofs_w_restrict, n_dofs_u_micro);
		M_micro_wv.reposition (w_var_restrict*n_dofs_u_restrict, v_var_micro*n_dofs_u_micro, n_dofs_w_restrict, n_dofs_v_micro);
		M_micro_ww.reposition (w_var_restrict*n_dofs_u_restrict, w_var_micro*n_dofs_u_micro, n_dofs_w_restrict, n_dofs_w_micro);

		Me_BIG.resize (n_dofs_restrict, n_dofs_BIG);

		M_BIG_uu.reposition (u_var_restrict*n_dofs_u_restrict, u_var_BIG*n_dofs_u_BIG, n_dofs_u_restrict, n_dofs_u_BIG);
		M_BIG_uv.reposition (u_var_restrict*n_dofs_u_restrict, v_var_BIG*n_dofs_u_BIG, n_dofs_u_restrict, n_dofs_v_BIG);
		M_BIG_uw.reposition (u_var_restrict*n_dofs_u_restrict, w_var_BIG*n_dofs_u_BIG, n_dofs_u_restrict, n_dofs_w_BIG);

		M_BIG_vu.reposition (v_var_restrict*n_dofs_u_restrict, u_var_BIG*n_dofs_u_BIG, n_dofs_v_restrict, n_dofs_u_BIG);
		M_BIG_vv.reposition (v_var_restrict*n_dofs_u_restrict, v_var_BIG*n_dofs_u_BIG, n_dofs_v_restrict, n_dofs_v_BIG);
		M_BIG_vw.reposition (v_var_restrict*n_dofs_u_restrict, w_var_BIG*n_dofs_u_BIG, n_dofs_v_restrict, n_dofs_w_BIG);

		M_BIG_wu.reposition (w_var_restrict*n_dofs_u_restrict, u_var_BIG*n_dofs_u_BIG, n_dofs_w_restrict, n_dofs_u_BIG);
		M_BIG_wv.reposition (w_var_restrict*n_dofs_u_restrict, v_var_BIG*n_dofs_u_BIG, n_dofs_w_restrict, n_dofs_v_BIG);
		M_BIG_ww.reposition (w_var_restrict*n_dofs_u_restrict, w_var_BIG*n_dofs_u_BIG, n_dofs_w_restrict, n_dofs_w_BIG);

		Me_restrict.resize (n_dofs_restrict, n_dofs_restrict);

		M_restrict_uu.reposition (u_var_restrict*n_dofs_u_restrict, u_var_restrict*n_dofs_u_restrict, n_dofs_u_restrict, n_dofs_u_restrict);
		M_restrict_uv.reposition (u_var_restrict*n_dofs_u_restrict, v_var_restrict*n_dofs_u_restrict, n_dofs_u_restrict, n_dofs_v_restrict);
		M_restrict_uw.reposition (u_var_restrict*n_dofs_u_restrict, w_var_restrict*n_dofs_u_restrict, n_dofs_u_restrict, n_dofs_w_restrict);

		M_restrict_vu.reposition (v_var_restrict*n_dofs_u_restrict, u_var_restrict*n_dofs_u_restrict, n_dofs_v_restrict, n_dofs_u_restrict);
		M_restrict_vv.reposition (v_var_restrict*n_dofs_u_restrict, v_var_restrict*n_dofs_u_restrict, n_dofs_v_restrict, n_dofs_v_restrict);
		M_restrict_vw.reposition (v_var_restrict*n_dofs_u_restrict, w_var_restrict*n_dofs_u_restrict, n_dofs_v_restrict, n_dofs_w_restrict);

		M_restrict_wu.reposition (w_var_restrict*n_dofs_u_restrict, u_var_restrict*n_dofs_u_restrict, n_dofs_w_restrict, n_dofs_u_restrict);
		M_restrict_wv.reposition (w_var_restrict*n_dofs_u_restrict, v_var_restrict*n_dofs_u_restrict, n_dofs_w_restrict, n_dofs_v_restrict);
		M_restrict_ww.reposition (w_var_restrict*n_dofs_u_restrict, w_var_restrict*n_dofs_u_restrict, n_dofs_w_restrict, n_dofs_w_restrict);

		// Restart the element
		fe_inter->reinit(elem_inter);
		n_quadrature_pts = fe_inter->n_quadrature_points();

		corrected_phi_BIG.resize(n_dofs_u_BIG,std::vector<libMesh::Real>(n_quadrature_pts,0));
		corrected_phi_restrict.resize(n_dofs_u_restrict,std::vector<libMesh::Real>(n_quadrature_pts,0));
		corrected_phi_micro.resize(n_dofs_u_micro,std::vector<libMesh::Real>(n_quadrature_pts,0));
	}

	// Vectors containing the lambda weights
	std::vector<std::vector<libMesh::Real> > 	lambda_weight_restrict(
													n_dofs_u_restrict,
													std::vector<libMesh::Real> (n_dofs_u_inter,0));
	std::vector<std::vector<libMesh::Real> > 	lambda_weight_BIG(
													n_dofs_u_BIG,
													std::vector<libMesh::Real> (n_dofs_u_inter,0));
	std::vector<std::vector<libMesh::Real> > 	lambda_weight_micro(
													n_dofs_u_micro,
													std::vector<libMesh::Real> (n_dofs_u_inter,0));

	// Initialize global matrix
	const unsigned int restrict_M = dof_map_restrict.n_dofs();

	const unsigned int BIG_N = dof_map_BIG.n_dofs();
	const unsigned int micro_N = dof_map_micro.n_dofs();

	const unsigned int restrict_M_local = dof_map_restrict.n_local_dofs();

	const unsigned int BIG_N_local = dof_map_BIG.n_local_dofs();
	const unsigned int micro_N_local = dof_map_micro.n_local_dofs();

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
			M_micro_uv.reposition (u_var_restrict*n_dofs_u_restrict, v_var_micro*n_dofs_u_micro, n_dofs_u_restrict, n_dofs_v_micro);
			M_micro_uw.reposition (u_var_restrict*n_dofs_u_restrict, w_var_micro*n_dofs_u_micro, n_dofs_u_restrict, n_dofs_w_micro);

			M_micro_vu.reposition (v_var_restrict*n_dofs_u_restrict, u_var_micro*n_dofs_u_micro, n_dofs_v_restrict, n_dofs_u_micro);
			M_micro_vv.reposition (v_var_restrict*n_dofs_u_restrict, v_var_micro*n_dofs_u_micro, n_dofs_v_restrict, n_dofs_v_micro);
			M_micro_vw.reposition (v_var_restrict*n_dofs_u_restrict, w_var_micro*n_dofs_u_micro, n_dofs_v_restrict, n_dofs_w_micro);

			M_micro_wu.reposition (w_var_restrict*n_dofs_u_restrict, u_var_micro*n_dofs_u_micro, n_dofs_w_restrict, n_dofs_u_micro);
			M_micro_wv.reposition (w_var_restrict*n_dofs_u_restrict, v_var_micro*n_dofs_u_micro, n_dofs_w_restrict, n_dofs_v_micro);
			M_micro_ww.reposition (w_var_restrict*n_dofs_u_restrict, w_var_micro*n_dofs_u_micro, n_dofs_w_restrict, n_dofs_w_micro);

			Me_BIG.resize (n_dofs_restrict, n_dofs_BIG);

			M_BIG_uu.reposition (u_var_restrict*n_dofs_u_restrict, u_var_BIG*n_dofs_u_BIG, n_dofs_u_restrict, n_dofs_u_BIG);
			M_BIG_uv.reposition (u_var_restrict*n_dofs_u_restrict, v_var_BIG*n_dofs_u_BIG, n_dofs_u_restrict, n_dofs_v_BIG);
			M_BIG_uw.reposition (u_var_restrict*n_dofs_u_restrict, w_var_BIG*n_dofs_u_BIG, n_dofs_u_restrict, n_dofs_w_BIG);

			M_BIG_vu.reposition (v_var_restrict*n_dofs_u_restrict, u_var_BIG*n_dofs_u_BIG, n_dofs_v_restrict, n_dofs_u_BIG);
			M_BIG_vv.reposition (v_var_restrict*n_dofs_u_restrict, v_var_BIG*n_dofs_u_BIG, n_dofs_v_restrict, n_dofs_v_BIG);
			M_BIG_vw.reposition (v_var_restrict*n_dofs_u_restrict, w_var_BIG*n_dofs_u_BIG, n_dofs_v_restrict, n_dofs_w_BIG);

			M_BIG_wu.reposition (w_var_restrict*n_dofs_u_restrict, u_var_BIG*n_dofs_u_BIG, n_dofs_w_restrict, n_dofs_u_BIG);
			M_BIG_wv.reposition (w_var_restrict*n_dofs_u_restrict, v_var_BIG*n_dofs_u_BIG, n_dofs_w_restrict, n_dofs_v_BIG);
			M_BIG_ww.reposition (w_var_restrict*n_dofs_u_restrict, w_var_BIG*n_dofs_u_BIG, n_dofs_w_restrict, n_dofs_w_BIG);

			Me_restrict.resize (n_dofs_restrict, n_dofs_restrict);

			M_restrict_uu.reposition (u_var_restrict*n_dofs_u_restrict, u_var_restrict*n_dofs_u_restrict, n_dofs_u_restrict, n_dofs_u_restrict);
			M_restrict_uv.reposition (u_var_restrict*n_dofs_u_restrict, v_var_restrict*n_dofs_u_restrict, n_dofs_u_restrict, n_dofs_v_restrict);
			M_restrict_uw.reposition (u_var_restrict*n_dofs_u_restrict, w_var_restrict*n_dofs_u_restrict, n_dofs_u_restrict, n_dofs_w_restrict);

			M_restrict_vu.reposition (v_var_restrict*n_dofs_u_restrict, u_var_restrict*n_dofs_u_restrict, n_dofs_v_restrict, n_dofs_u_restrict);
			M_restrict_vv.reposition (v_var_restrict*n_dofs_u_restrict, v_var_restrict*n_dofs_u_restrict, n_dofs_v_restrict, n_dofs_v_restrict);
			M_restrict_vw.reposition (v_var_restrict*n_dofs_u_restrict, w_var_restrict*n_dofs_u_restrict, n_dofs_v_restrict, n_dofs_w_restrict);

			M_restrict_wu.reposition (w_var_restrict*n_dofs_u_restrict, u_var_restrict*n_dofs_u_restrict, n_dofs_w_restrict, n_dofs_u_restrict);
			M_restrict_wv.reposition (w_var_restrict*n_dofs_u_restrict, v_var_restrict*n_dofs_u_restrict, n_dofs_w_restrict, n_dofs_v_restrict);
			M_restrict_ww.reposition (w_var_restrict*n_dofs_u_restrict, w_var_restrict*n_dofs_u_restrict, n_dofs_w_restrict, n_dofs_w_restrict);

			// Set up corrected shape vectors
			const libMesh::Elem* elem_inter = mesh_inter.elem(0);
			fe_inter->reinit(elem_inter);
			n_quadrature_pts = fe_inter->n_quadrature_points();

			corrected_phi_BIG.resize(n_dofs_u_BIG,std::vector<libMesh::Real>(n_quadrature_pts,0));
			corrected_phi_restrict.resize(n_dofs_u_restrict,std::vector<libMesh::Real>(n_quadrature_pts,0));
			corrected_phi_micro.resize(n_dofs_u_micro,std::vector<libMesh::Real>(n_quadrature_pts,0));
		}

		Me_restrict.zero();
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

			// Restart the element
			fe_inter->reinit(elem_inter);

			get_lambdas(dim_BIG, fe_type_BIG, elem_BIG,
							quad_points_inter,
							quad_points_reference,
							lambda_weight_BIG);

			get_lambdas(dim_micro, fe_type_micro, elem_micro,
							quad_points_inter,
							quad_points_reference,
							lambda_weight_micro);

			get_lambdas(dim_restrict, fe_type_restrict, elem_restrict,
							quad_points_inter,
							quad_points_reference,
							lambda_weight_restrict);

			set_corrected_shapes(lambda_weight_BIG,phi_inter,corrected_phi_BIG);
			set_corrected_shapes(lambda_weight_micro,phi_inter,corrected_phi_micro);
			set_corrected_shapes(lambda_weight_restrict,phi_inter,corrected_phi_restrict);

			// For each quadrature point determinate the sub-matrices elements
			for (unsigned int qp=0; qp<qrule_inter.n_points(); qp++)
			{
				// Restrict -> micro coupling
				L2_Coupling(M_micro_uu,qp,corrected_phi_BIG,corrected_phi_micro,
							n_dofs_u_restrict,n_dofs_u_micro,JxW,coupling_const);

				L2_Coupling(M_micro_vv,qp,corrected_phi_BIG,corrected_phi_micro,
							n_dofs_v_restrict,n_dofs_v_micro,JxW,coupling_const);

				L2_Coupling(M_micro_ww,qp,corrected_phi_BIG,corrected_phi_micro,
							n_dofs_w_restrict,n_dofs_w_micro,JxW,coupling_const);

				// Restrict -> BIG coupling
				L2_Coupling(M_BIG_uu,qp,corrected_phi_restrict,corrected_phi_BIG,
							n_dofs_u_restrict,n_dofs_u_BIG,JxW,coupling_const);

				L2_Coupling(M_BIG_vv,qp,corrected_phi_restrict,corrected_phi_BIG,
							n_dofs_v_restrict,n_dofs_v_BIG,JxW,coupling_const);

				L2_Coupling(M_BIG_ww,qp,corrected_phi_restrict,corrected_phi_BIG,
							n_dofs_w_restrict,n_dofs_w_BIG,JxW,coupling_const);

				// Restrict -> Restrict coupling
				L2_Coupling(M_restrict_uu,qp,corrected_phi_restrict,corrected_phi_restrict,
							n_dofs_u_restrict,n_dofs_u_restrict,JxW,coupling_const);

				L2_Coupling(M_restrict_vv,qp,corrected_phi_restrict,corrected_phi_restrict,
							n_dofs_v_restrict,n_dofs_v_restrict,JxW,coupling_const);

				L2_Coupling(M_restrict_ww,qp,corrected_phi_restrict,corrected_phi_restrict,
							n_dofs_w_restrict,n_dofs_w_restrict,JxW,coupling_const);
			}
		}

		couplingMatrix_restrict_micro.add_matrix(Me_micro,dof_indices_restrict,dof_indices_micro);
		couplingMatrix_restrict_BIG.add_matrix(Me_BIG,dof_indices_restrict,dof_indices_BIG);
		couplingMatrix_restrict_restrict.add_matrix(Me_restrict,dof_indices_restrict,dof_indices_restrict);
	}

	couplingMatrix_restrict_micro.close();
	couplingMatrix_restrict_BIG.close();
	couplingMatrix_restrict_restrict.close();
};
