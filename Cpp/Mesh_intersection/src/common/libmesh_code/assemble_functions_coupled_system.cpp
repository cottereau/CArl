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
	libMesh::LinearImplicitSystem& volume_restrict_system = restrict_eq_system.get_system<libMesh::LinearImplicitSystem>("Elasticity");
	libMesh::LinearImplicitSystem& volume_BIG_system = BIG_eq_system.get_system<libMesh::LinearImplicitSystem>("Elasticity");
	libMesh::LinearImplicitSystem& volume_micro_system = micro_eq_system.get_system<libMesh::LinearImplicitSystem>("Elasticity");
	libMesh::LinearImplicitSystem& volume_inter_system = inter_eq_system.get_system<libMesh::LinearImplicitSystem>("Elasticity");

	libMesh_fe_addresses_3 restrict_addresses(volume_restrict_system);
	libMesh_fe_addresses_3 BIG_addresses(volume_BIG_system);
	libMesh_fe_addresses_3 micro_addresses(volume_micro_system);
	libMesh_fe_addresses_3 inter_addresses(volume_inter_system);

	// Vector that will keep the quadrature points
	const std::vector<libMesh::Point>& quad_points_inter = inter_addresses.fe_unique_ptr->get_xyz();
	std::vector<libMesh::Point> quad_points_reference;

	// Jacobians
	const std::vector<libMesh::Real>& JxW = inter_addresses.fe_unique_ptr->get_JxW();

	// Shape functions
	const std::vector<std::vector<libMesh::Real> >& phi_inter = inter_addresses.fe_unique_ptr->get_phi();
	std::vector<std::vector<libMesh::Real> > corrected_phi_BIG;
	std::vector<std::vector<libMesh::Real> > corrected_phi_micro;
	std::vector<std::vector<libMesh::Real> > corrected_phi_restrict;

	unsigned int n_quadrature_pts = 0;

	// Addresses to the matrices
	libMesh::PetscMatrix<libMesh::Number>& couplingMatrix_restrict_BIG = * m_couplingMatrixMap_restrict_BIG[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& couplingMatrix_restrict_micro = * m_couplingMatrixMap_restrict_micro[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& couplingMatrix_restrict_restrict = * m_couplingMatrixMap_restrict_restrict[micro_name];

	// DoF vectors and ranges
	restrict_addresses.set_DoFs();
	BIG_addresses.set_DoFs();
	micro_addresses.set_DoFs();
	inter_addresses.set_DoFs();

	// Local matrix
	coupling_matrices_3 Me_restrict_micro;
	coupling_matrices_3 Me_restrict_BIG;
	coupling_matrices_3 Me_restrict_restrict;

	//    If all elements are of the same type, do the index "extraction",
	// the matrices resizes and repositions here
	if(bSameElemsType)
	{
		Me_restrict_micro.set_matrices(restrict_addresses,micro_addresses);
		Me_restrict_BIG.set_matrices(restrict_addresses,BIG_addresses);
		Me_restrict_restrict.set_matrices(restrict_addresses,restrict_addresses);

		// Restart the element
		const libMesh::Elem* elem_inter = micro_addresses.mesh.elem(0);
		inter_addresses.fe_unique_ptr->reinit(elem_inter);
		n_quadrature_pts = inter_addresses.fe_unique_ptr->n_quadrature_points();

		corrected_phi_BIG.resize(BIG_addresses.n_dofs_u,std::vector<libMesh::Real>(n_quadrature_pts,0));
		corrected_phi_restrict.resize(restrict_addresses.n_dofs_u,std::vector<libMesh::Real>(n_quadrature_pts,0));
		corrected_phi_micro.resize(micro_addresses.n_dofs_u,std::vector<libMesh::Real>(n_quadrature_pts,0));
	}

	// Vectors containing the lambda weights
	std::vector<std::vector<libMesh::Real> > 	lambda_weight_restrict(
													restrict_addresses.n_dofs_u,
													std::vector<libMesh::Real> (inter_addresses.n_dofs_u,0));
	std::vector<std::vector<libMesh::Real> > 	lambda_weight_BIG(
													BIG_addresses.n_dofs_u,
													std::vector<libMesh::Real> (inter_addresses.n_dofs_u,0));
	std::vector<std::vector<libMesh::Real> > 	lambda_weight_micro(
													micro_addresses.n_dofs_u,
													std::vector<libMesh::Real> (inter_addresses.n_dofs_u,0));

	// Initialize global matrix
	const unsigned int restrict_M = restrict_addresses.dof_map.n_dofs();

	const unsigned int BIG_N = BIG_addresses.dof_map.n_dofs();
	const unsigned int micro_N = micro_addresses.dof_map.n_dofs();

	const unsigned int restrict_M_local = restrict_addresses.dof_map.n_local_dofs();

	const unsigned int BIG_N_local = BIG_addresses.dof_map.n_local_dofs();
	const unsigned int micro_N_local = micro_addresses.dof_map.n_local_dofs();

	couplingMatrix_restrict_micro.init(		restrict_M, micro_N,
											restrict_M_local, micro_N_local,
											micro_N_local,micro_N - micro_N_local);
	couplingMatrix_restrict_BIG.init(		restrict_M, BIG_N,
											restrict_M_local, BIG_N_local,
											BIG_N_local,BIG_N - BIG_N_local);

	couplingMatrix_restrict_restrict.attach_dof_map(restrict_addresses.dof_map);
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

		const libMesh::Elem* elem_restrict = restrict_addresses.mesh.elem(elem_restrict_idx);
		const libMesh::Elem* elem_BIG = BIG_addresses.mesh.elem(elem_BIG_idx);
		const libMesh::Elem* elem_micro = micro_addresses.mesh.elem(elem_micro_idx);

		restrict_addresses.set_DoFs(elem_restrict_idx);
		BIG_addresses.set_DoFs(elem_BIG_idx);
		micro_addresses.set_DoFs(elem_micro_idx);

		// Resize dense matrix, if needed
		if(!bSameElemsType)
		{
			Me_restrict_micro.set_matrices(restrict_addresses,micro_addresses);
			Me_restrict_BIG.set_matrices(restrict_addresses,BIG_addresses);
			Me_restrict_restrict.set_matrices(restrict_addresses,restrict_addresses);

			// Set up corrected shape vectors
			const libMesh::Elem* elem_inter = micro_addresses.mesh.elem(0);
			inter_addresses.fe_unique_ptr->reinit(elem_inter);
			n_quadrature_pts = inter_addresses.fe_unique_ptr->n_quadrature_points();

			corrected_phi_BIG.resize(BIG_addresses.n_dofs_u,std::vector<libMesh::Real>(n_quadrature_pts,0));
			corrected_phi_restrict.resize(restrict_addresses.n_dofs_u,std::vector<libMesh::Real>(n_quadrature_pts,0));
			corrected_phi_micro.resize(micro_addresses.n_dofs_u,std::vector<libMesh::Real>(n_quadrature_pts,0));
		}

		Me_restrict_micro.zero();
		Me_restrict_BIG.zero();
		Me_restrict_restrict.zero();

		// Now iterate over the intersections
		auto inter_idx_range = intersection_table_inter.equal_range(iii);

		for(	it_inter_idx = inter_idx_range.first;
				it_inter_idx != inter_idx_range.second;
				++it_inter_idx)
		{
			// Get the intersection mesh pointer
			elem_inter_idx = it_inter_idx->second;
			const libMesh::Elem* elem_inter = inter_addresses.mesh.elem(elem_inter_idx);

			// Restart the element
			inter_addresses.fe_unique_ptr->reinit(elem_inter);

			get_lambdas(BIG_addresses.dim, BIG_addresses.fe_type, elem_BIG,
							quad_points_inter,
							quad_points_reference,
							lambda_weight_BIG);

			get_lambdas(micro_addresses.dim, micro_addresses.fe_type, elem_micro,
							quad_points_inter,
							quad_points_reference,
							lambda_weight_micro);

			get_lambdas(restrict_addresses.dim, restrict_addresses.fe_type, elem_restrict,
							quad_points_inter,
							quad_points_reference,
							lambda_weight_restrict);

			set_corrected_shapes(lambda_weight_BIG,phi_inter,corrected_phi_BIG);
			set_corrected_shapes(lambda_weight_micro,phi_inter,corrected_phi_micro);
			set_corrected_shapes(lambda_weight_restrict,phi_inter,corrected_phi_restrict);

			// For each quadrature point determinate the sub-matrices elements
			for (unsigned int qp=0; qp < inter_addresses.qrule.n_points(); qp++)
			{
				// Restrict -> micro coupling
				L2_Coupling(Me_restrict_micro.Me_uu,qp,corrected_phi_restrict,corrected_phi_micro,
							restrict_addresses.n_dofs_u,micro_addresses.n_dofs_u,JxW,coupling_const);

				L2_Coupling(Me_restrict_micro.Me_vv,qp,corrected_phi_restrict,corrected_phi_micro,
							restrict_addresses.n_dofs_v,micro_addresses.n_dofs_v,JxW,coupling_const);

				L2_Coupling(Me_restrict_micro.Me_ww,qp,corrected_phi_restrict,corrected_phi_micro,
							restrict_addresses.n_dofs_w,micro_addresses.n_dofs_w,JxW,coupling_const);

				// Restrict -> BIG coupling
				L2_Coupling(Me_restrict_BIG.Me_uu,qp,corrected_phi_restrict,corrected_phi_BIG,
							restrict_addresses.n_dofs_u,BIG_addresses.n_dofs_u,JxW,coupling_const);

				L2_Coupling(Me_restrict_BIG.Me_vv,qp,corrected_phi_restrict,corrected_phi_BIG,
							restrict_addresses.n_dofs_v,BIG_addresses.n_dofs_v,JxW,coupling_const);

				L2_Coupling(Me_restrict_BIG.Me_ww,qp,corrected_phi_restrict,corrected_phi_BIG,
							restrict_addresses.n_dofs_w,BIG_addresses.n_dofs_w,JxW,coupling_const);

				// Restrict -> Restrict coupling
				L2_Coupling(Me_restrict_restrict.Me_uu,qp,corrected_phi_restrict,corrected_phi_restrict,
							restrict_addresses.n_dofs_u,restrict_addresses.n_dofs_u,JxW,coupling_const);

				L2_Coupling(Me_restrict_restrict.Me_vv,qp,corrected_phi_restrict,corrected_phi_restrict,
							restrict_addresses.n_dofs_v,restrict_addresses.n_dofs_v,JxW,coupling_const);

				L2_Coupling(Me_restrict_restrict.Me_ww,qp,corrected_phi_restrict,corrected_phi_restrict,
							restrict_addresses.n_dofs_w,restrict_addresses.n_dofs_w,JxW,coupling_const);
			}
		}

		couplingMatrix_restrict_micro.add_matrix(Me_restrict_micro.Me,restrict_addresses.dof_indices,micro_addresses.dof_indices);
		couplingMatrix_restrict_BIG.add_matrix(Me_restrict_BIG.Me,restrict_addresses.dof_indices,BIG_addresses.dof_indices);
		couplingMatrix_restrict_restrict.add_matrix(Me_restrict_restrict.Me,restrict_addresses.dof_indices,restrict_addresses.dof_indices);
	}

	couplingMatrix_restrict_micro.close();
	couplingMatrix_restrict_BIG.close();
	couplingMatrix_restrict_restrict.close();
};
