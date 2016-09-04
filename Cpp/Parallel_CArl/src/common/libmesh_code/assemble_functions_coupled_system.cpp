#include "coupled_system.h"

void carl::coupled_system::assemble_coupling_matrices(
		const std::string BIG_name,
		const std::string micro_name,
		const std::string inter_name,
		const std::string mediator_name,

		std::unordered_map<int, int>& equivalence_table_mediator_BIG,
		std::vector<std::pair<int, int> >& intersection_table_mediator_micro,
		std::unordered_multimap<int, int>& intersection_table_inter,

		double coupling_const,

		bool using_same_mesh_mediator_A,
		bool bSameElemsType)
{
	// Addresses to the eq. systems
	libMesh::EquationSystems& mediator_eq_system =
			*m_mediator_EquationSystemMap[mediator_name];
	libMesh::EquationSystems& BIG_eq_system = *m_BIG_EquationSystem.second;
	libMesh::EquationSystems& micro_eq_system =
			*m_micro_EquationSystemMap[micro_name];
	libMesh::EquationSystems& inter_eq_system =
			*m_inter_EquationSystemMap[inter_name];

	// First, test if all the systems have an acoustic model and variable set
	homemade_assert_msg(micro_eq_system.has_system("VolTest"),
			" Micro equation systems missing \"VolTest\" system!");
	homemade_assert_msg(BIG_eq_system.has_system("VolTest"),
			" Macro equation systems missing \"VolTest\" system!");
	homemade_assert_msg(inter_eq_system.has_system("VolTest"),
			" Intersection equation systems missing \"VolTest\" system!");
	homemade_assert_msg(mediator_eq_system.has_system("VolTest"),
			" Mediatored equation systems missing \"VolTest\" system!");

	// Systems and vars
	libMesh::LinearImplicitSystem& volume_mediator_system =
			mediator_eq_system.get_system<libMesh::LinearImplicitSystem>(
					"VolTest");
	libMesh::LinearImplicitSystem& volume_BIG_system = BIG_eq_system.get_system<
			libMesh::LinearImplicitSystem>("VolTest");
	libMesh::LinearImplicitSystem& volume_micro_system =
			micro_eq_system.get_system<libMesh::LinearImplicitSystem>(
					"VolTest");
	libMesh::LinearImplicitSystem& volume_inter_system =
			inter_eq_system.get_system<libMesh::LinearImplicitSystem>(
					"VolTest");

	// Addresses to the meshes
	const libMesh::MeshBase& mesh_mediator = mediator_eq_system.get_mesh();
	const libMesh::MeshBase& mesh_BIG = BIG_eq_system.get_mesh();
	const libMesh::MeshBase& mesh_micro = micro_eq_system.get_mesh();
	const libMesh::MeshBase& mesh_inter = inter_eq_system.get_mesh();

	// Variable indices and dimensions
	const unsigned int dim_mediator = mesh_mediator.mesh_dimension();
	const unsigned int dim_BIG = mesh_BIG.mesh_dimension();
	const unsigned int dim_micro = mesh_micro.mesh_dimension();
	const unsigned int dim_inter = mesh_inter.mesh_dimension();

	const unsigned int vol_test_mediator =
			volume_mediator_system.variable_number("SillyVar");
	const unsigned int vol_test_BIG = volume_BIG_system.variable_number(
			"SillyVar");
	const unsigned int vol_test_micro = volume_micro_system.variable_number(
			"SillyVar");
	const unsigned int vol_test_inter = volume_inter_system.variable_number(
			"SillyVar");

	// Addresses to the DoF maps
	const libMesh::DofMap& dof_map_mediator =
			volume_mediator_system.get_dof_map();
	const libMesh::DofMap& dof_map_BIG = volume_BIG_system.get_dof_map();
	const libMesh::DofMap& dof_map_micro = volume_micro_system.get_dof_map();
	const libMesh::DofMap& dof_map_inter = volume_inter_system.get_dof_map();

	// Finite elements (FE) declaration
	libMesh::FEType fe_type_mediator = dof_map_mediator.variable_type(
			vol_test_mediator);
	libMesh::FEType fe_type_BIG = dof_map_BIG.variable_type(vol_test_BIG);
	libMesh::FEType fe_type_micro = dof_map_micro.variable_type(vol_test_micro);
	libMesh::FEType fe_type_inter = dof_map_inter.variable_type(vol_test_inter);

	// Test if we are using linear elements
	homemade_assert_msg(fe_type_mediator.order == libMesh::FIRST,
			" Mediator system is not linear!");
	homemade_assert_msg(fe_type_BIG.order == libMesh::FIRST,
			" Macro system is not linear!");
	homemade_assert_msg(fe_type_micro.order == libMesh::FIRST,
			" Micro system is not linear!");
	homemade_assert_msg(fe_type_inter.order == libMesh::FIRST,
			" Intersection system is not linear!");

	// Set up FE bases and quadratures
	libMesh::UniquePtr<libMesh::FEBase> fe_mediator(
			libMesh::FEBase::build(dim_mediator, fe_type_mediator));
	libMesh::UniquePtr<libMesh::FEBase> fe_BIG(
			libMesh::FEBase::build(dim_BIG, fe_type_BIG));
	libMesh::UniquePtr<libMesh::FEBase> fe_micro(
			libMesh::FEBase::build(dim_micro, fe_type_micro));
	libMesh::UniquePtr<libMesh::FEBase> fe_inter(
			libMesh::FEBase::build(dim_inter, fe_type_inter));

	libMesh::QGauss qrule_mediator(dim_mediator,
			fe_type_mediator.default_quadrature_order());
	fe_mediator->attach_quadrature_rule(&qrule_mediator);

	libMesh::QGauss qrule_inter(dim_inter,
			fe_type_inter.default_quadrature_order());
	fe_inter->attach_quadrature_rule(&qrule_inter);

	// Vector that will keep the quadrature points
	const std::vector<libMesh::Point>& quad_points_inter = fe_inter->get_xyz();
	std::vector<libMesh::Point> quad_points_reference;

	// Jacobians
	const std::vector<libMesh::Real>& JxW = fe_inter->get_JxW();

	// Shape functions
	const std::vector<std::vector<libMesh::Real> >& phi_inter =
			fe_inter->get_phi();
	std::vector<std::vector<libMesh::Real> > corrected_phi_BIG;
	std::vector<std::vector<libMesh::Real> > corrected_phi_micro;
	std::vector<std::vector<libMesh::Real> > corrected_phi_mediator;
	unsigned int n_quadrature_pts = 0;

	// Addresses to the matrices
	libMesh::PetscMatrix<libMesh::Number>& couplingMatrix_mediator_BIG =
			*m_couplingMatrixMap_mediator_BIG[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& couplingMatrix_mediator_micro =
			*m_couplingMatrixMap_mediator_micro[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& couplingMatrix_mediator_mediator =
			*m_couplingMatrixMap_mediator_mediator[micro_name];

	// Dof vectors and ranges
	std::vector<libMesh::dof_id_type> dof_indices_mediator;
	std::vector<libMesh::dof_id_type> dof_indices_BIG;
	std::vector<libMesh::dof_id_type> dof_indices_micro;
	std::vector<libMesh::dof_id_type> dof_indices_inter;

	// Number of DoF's for each system and each variable
	unsigned int n_dofs_mediator;
	unsigned int n_dofs_micro;
	unsigned int n_dofs_BIG;
	unsigned int n_dofs_inter;

	// Local matrix
	libMesh::DenseMatrix<libMesh::Number> Me_micro;
	libMesh::DenseMatrix<libMesh::Number> Me_BIG;
	libMesh::DenseMatrix<libMesh::Number> Me_mediator;

	//    If all elements are of the same type, do the index "extraction",
	// the matrices resizes and repositions here
	const libMesh::Elem* dummy_elem_mediator = mesh_mediator.elem(0);
	dof_map_mediator.dof_indices(dummy_elem_mediator, dof_indices_mediator);
	n_dofs_mediator = dof_indices_mediator.size();

	const libMesh::Elem* dummy_elem_BIG = mesh_BIG.elem(0);
	dof_map_BIG.dof_indices(dummy_elem_BIG, dof_indices_BIG);
	n_dofs_BIG = dof_indices_BIG.size();

	const libMesh::Elem* dummy_elem_micro = mesh_micro.elem(0);
	dof_map_micro.dof_indices(dummy_elem_micro, dof_indices_micro);
	n_dofs_micro = dof_indices_micro.size();

	const libMesh::Elem* dummy_elem_inter = mesh_inter.elem(0);
	dof_map_inter.dof_indices(dummy_elem_inter, dof_indices_inter);
	n_dofs_inter = dof_indices_inter.size();

	if (bSameElemsType)
	{
		// Resize matrices
		Me_BIG.resize(n_dofs_mediator, n_dofs_BIG);
		Me_micro.resize(n_dofs_mediator, n_dofs_micro);
		Me_mediator.resize(n_dofs_mediator, n_dofs_mediator);

		// Set up corrected shape vectors
		fe_inter->reinit(dummy_elem_inter);
		n_quadrature_pts = fe_inter->n_quadrature_points();

		corrected_phi_BIG.resize(n_dofs_BIG,
				std::vector<libMesh::Real>(n_quadrature_pts, 0));
		corrected_phi_mediator.resize(n_dofs_mediator,
				std::vector<libMesh::Real>(n_quadrature_pts, 0));
		corrected_phi_micro.resize(n_dofs_micro,
				std::vector<libMesh::Real>(n_quadrature_pts, 0));
	}

	// Vectors containing the lambda weights
	std::vector<std::vector<libMesh::Real> > lambda_weight_mediator(
			n_dofs_mediator, std::vector<libMesh::Real>(n_dofs_inter, 0));
	std::vector<std::vector<libMesh::Real> > lambda_weight_BIG(n_dofs_BIG,
			std::vector<libMesh::Real>(n_dofs_inter, 0));
	std::vector<std::vector<libMesh::Real> > lambda_weight_micro(n_dofs_micro,
			std::vector<libMesh::Real>(n_dofs_inter, 0));

	// Initialize global matrix
	const unsigned int mediator_M = dof_map_mediator.n_dofs();

	const unsigned int BIG_N = dof_map_BIG.n_dofs();
	const unsigned int micro_N = dof_map_micro.n_dofs();

	const unsigned int mediator_M_local = dof_map_mediator.n_local_dofs();

	const unsigned int BIG_N_local = dof_map_BIG.n_local_dofs();
	const unsigned int micro_N_local = dof_map_micro.n_local_dofs();

	couplingMatrix_mediator_micro.init(mediator_M, micro_N, mediator_M_local,
			micro_N_local, micro_N_local, micro_N - micro_N_local);
	couplingMatrix_mediator_BIG.init(mediator_M, BIG_N, mediator_M_local,
			BIG_N_local, BIG_N_local, BIG_N - BIG_N_local);

	couplingMatrix_mediator_mediator.attach_dof_map(dof_map_mediator);
	couplingMatrix_mediator_mediator.init();

	// Intersection indexes and iterators
	int nb_of_intersections = intersection_table_mediator_micro.size();
	int elem_mediator_idx = -1;

	int elem_BIG_idx = -1;
	int elem_micro_idx = -1;
	int elem_inter_idx = -1;
	std::unordered_multimap<int, int>::iterator it_inter_idx;

	// For each intersection
	for (int iii = 0; iii < nb_of_intersections; ++iii)
	{
		// Get the mediator and micro element pointers
		elem_mediator_idx = intersection_table_mediator_micro[iii].first;
		elem_micro_idx = intersection_table_mediator_micro[iii].second;
		if (using_same_mesh_mediator_A)
		{
			elem_BIG_idx = elem_mediator_idx;
		}
		else
		{
			elem_BIG_idx = equivalence_table_mediator_BIG[elem_mediator_idx];
		}

		const libMesh::Elem* elem_mediator = mesh_mediator.elem(
				elem_mediator_idx);
		const libMesh::Elem* elem_BIG = mesh_BIG.elem(elem_BIG_idx);
		const libMesh::Elem* elem_micro = mesh_micro.elem(elem_micro_idx);

		// And their dof map indices
		dof_map_mediator.dof_indices(elem_mediator, dof_indices_mediator);

		dof_map_micro.dof_indices(elem_micro, dof_indices_micro);

		if (using_same_mesh_mediator_A)
		{
			dof_indices_BIG = dof_indices_mediator;
		}
		else
		{
			dof_map_BIG.dof_indices(elem_BIG, dof_indices_BIG);
		}

		// Resize dense matrix, if needed
		if (!bSameElemsType)
		{
			n_dofs_mediator = dof_indices_mediator.size();
			n_dofs_BIG = dof_indices_BIG.size();
			n_dofs_micro = dof_indices_micro.size();

			// Resize matrices
			Me_micro.resize(n_dofs_mediator, n_dofs_micro);
			Me_BIG.resize(n_dofs_mediator, n_dofs_BIG);
			Me_BIG.resize(n_dofs_mediator, n_dofs_mediator);

			// Set up corrected shape vectors
			const libMesh::Elem* elem_inter = mesh_inter.elem(0);
			fe_inter->reinit(elem_inter);
			n_quadrature_pts = fe_inter->n_quadrature_points();

			corrected_phi_BIG.resize(n_dofs_BIG,
					std::vector<libMesh::Real>(n_quadrature_pts, 0));
			corrected_phi_mediator.resize(n_dofs_mediator,
					std::vector<libMesh::Real>(n_quadrature_pts, 0));
			corrected_phi_micro.resize(n_dofs_micro,
					std::vector<libMesh::Real>(n_quadrature_pts, 0));

		}

		Me_micro.zero();
		Me_BIG.zero();
		Me_mediator.zero();

		// Now iterate over the intersections
		auto inter_idx_range = intersection_table_inter.equal_range(iii);

		for (it_inter_idx = inter_idx_range.first;
				it_inter_idx != inter_idx_range.second; ++it_inter_idx)
		{
			// Get the intersection mesh pointer
			elem_inter_idx = it_inter_idx->second;

			const libMesh::Elem* elem_inter = mesh_inter.elem(elem_inter_idx);

			// Restart the elements
			fe_inter->reinit(elem_inter);

			get_lambdas(dim_BIG, fe_type_BIG, elem_BIG, quad_points_inter,
					quad_points_reference, lambda_weight_BIG);

			get_lambdas(dim_micro, fe_type_micro, elem_micro, quad_points_inter,
					quad_points_reference, lambda_weight_micro);

			get_lambdas(dim_mediator, fe_type_mediator, elem_mediator,
					quad_points_inter, quad_points_reference,
					lambda_weight_mediator);

			set_corrected_shapes(lambda_weight_BIG, phi_inter,
					corrected_phi_BIG);
			set_corrected_shapes(lambda_weight_micro, phi_inter,
					corrected_phi_micro);
			set_corrected_shapes(lambda_weight_mediator, phi_inter,
					corrected_phi_mediator);

			// For each quadrature point determinate the sub-matrices elements
			for (unsigned int qp = 0; qp < qrule_inter.n_points(); qp++)
			{
				// Mediator -> micro coupling
				L2_Coupling(Me_micro, qp, corrected_phi_mediator,
						corrected_phi_micro, n_dofs_mediator, n_dofs_micro, JxW,
						coupling_const);

				// Mediator -> BIG coupling
				L2_Coupling(Me_BIG, qp, corrected_phi_mediator,
						corrected_phi_BIG, n_dofs_mediator, n_dofs_BIG, JxW,
						coupling_const);

				// Mediator -> mediator coupling
				L2_Coupling(Me_mediator, qp, corrected_phi_mediator,
						corrected_phi_mediator, n_dofs_mediator,
						n_dofs_mediator, JxW, coupling_const);
			}
		}

		couplingMatrix_mediator_micro.add_matrix(Me_micro, dof_indices_mediator,
				dof_indices_micro);
		couplingMatrix_mediator_BIG.add_matrix(Me_BIG, dof_indices_mediator,
				dof_indices_BIG);
		couplingMatrix_mediator_mediator.add_matrix(Me_mediator,
				dof_indices_mediator, dof_indices_mediator);
	}

	couplingMatrix_mediator_micro.close();
	couplingMatrix_mediator_BIG.close();
	couplingMatrix_mediator_mediator.close();
};

void carl::coupled_system::assemble_coupling_elasticity_3D(
		const std::string BIG_name,
		const std::string micro_name,
		const std::string inter_name,
		const std::string mediator_name,

		std::unordered_map<int, int>& equivalence_table_mediator_BIG,
		std::vector<std::pair<int, int> >& intersection_table_mediator_micro,
		std::unordered_multimap<int, int>& intersection_table_inter,

		bool using_same_mesh_mediator_A,
		bool bSameElemsType)
{
	// Addresses to the eq. systems
	libMesh::EquationSystems& mediator_eq_system =
			*m_mediator_EquationSystemMap[mediator_name];
	libMesh::EquationSystems& BIG_eq_system = *m_BIG_EquationSystem.second;
	libMesh::EquationSystems& micro_eq_system =
			*m_micro_EquationSystemMap[micro_name];
	libMesh::EquationSystems& inter_eq_system =
			*m_inter_EquationSystemMap[inter_name];

	// First, test if all the systems have an elasticity model and variable set
	homemade_assert_msg(micro_eq_system.has_system("Elasticity") || micro_eq_system.has_system("NonlinearElasticity"),
			" Micro equation systems missing \"Elasticity\" system!");
	homemade_assert_msg(BIG_eq_system.has_system("Elasticity"),
			" Macro equation systems missing \"Elasticity\" system!");
	homemade_assert_msg(inter_eq_system.has_system("Elasticity"),
			" Intersection equation systems missing \"Elasticity\" system!");
	homemade_assert_msg(mediator_eq_system.has_system("Elasticity"),
			" Mediatored equation systems missing \"Elasticity\" system!");

	// Systems and vars
	libMesh::LinearImplicitSystem& volume_mediator_system =
			mediator_eq_system.get_system<libMesh::LinearImplicitSystem>(
					"Elasticity");
	libMesh::LinearImplicitSystem& volume_BIG_system = BIG_eq_system.get_system<
			libMesh::LinearImplicitSystem>("Elasticity");
	libMesh::LinearImplicitSystem& volume_micro_system =
			micro_eq_system.get_system<libMesh::LinearImplicitSystem>(
					"Elasticity");
	libMesh::LinearImplicitSystem& volume_inter_system =
			inter_eq_system.get_system<libMesh::LinearImplicitSystem>(
					"Elasticity");

	libMesh_fe_addresses_3 mediator_addresses(volume_mediator_system);
	libMesh_fe_addresses_3 BIG_addresses(volume_BIG_system);
	libMesh_fe_addresses_3 micro_addresses(volume_micro_system);
	libMesh_fe_addresses_3 inter_addresses(volume_inter_system);

	// Vector that will keep the quadrature points
	const std::vector<libMesh::Point>& quad_points_inter =
			inter_addresses.fe_unique_ptr->get_xyz();
	std::vector<libMesh::Point> quad_points_reference;

	// Jacobians
	const std::vector<libMesh::Real>& JxW =
			inter_addresses.fe_unique_ptr->get_JxW();

	// Shape functions
	const std::vector<std::vector<libMesh::Real> >& phi_inter =
			inter_addresses.fe_unique_ptr->get_phi();
	std::vector<std::vector<libMesh::Real> > corrected_phi_BIG;
	std::vector<std::vector<libMesh::Real> > corrected_phi_micro;
	std::vector<std::vector<libMesh::Real> > corrected_phi_mediator;

	// Shape functions gradients
	const std::vector<std::vector<libMesh::RealGradient> >& dphi_inter =
			inter_addresses.fe_unique_ptr->get_dphi();
	std::vector<std::vector<libMesh::RealGradient> > corrected_dphi_BIG;
	std::vector<std::vector<libMesh::RealGradient> > corrected_dphi_micro;
	std::vector<std::vector<libMesh::RealGradient> > corrected_dphi_mediator;

	unsigned int n_quadrature_pts = 0;

	// Addresses to the matrices
	libMesh::PetscMatrix<libMesh::Number>& couplingMatrix_mediator_BIG =
			*m_couplingMatrixMap_mediator_BIG[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& couplingMatrix_mediator_micro =
			*m_couplingMatrixMap_mediator_micro[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& couplingMatrix_mediator_mediator =
			*m_couplingMatrixMap_mediator_mediator[micro_name];

	// Values of the coupling constants
	double L2_coupling_const = m_coupling_constantMap[micro_name]
			/ (m_coupling_lengthMap[micro_name]
					* m_coupling_lengthMap[micro_name]);
	double H1_coupling_const = m_coupling_constantMap[micro_name];

	// DoF vectors and ranges
	mediator_addresses.set_DoFs();
	BIG_addresses.set_DoFs();
	micro_addresses.set_DoFs();
	inter_addresses.set_DoFs();

	// Local matrix
	coupling_matrices_3 Me_mediator_micro;
	coupling_matrices_3 Me_mediator_BIG;
	coupling_matrices_3 Me_mediator_mediator;

	//    If all elements are of the same type, do the index "extraction",
	// the matrices resizes and repositions here
	if (bSameElemsType)
	{
		Me_mediator_micro.set_matrices(mediator_addresses, micro_addresses);
		Me_mediator_BIG.set_matrices(mediator_addresses, BIG_addresses);
		Me_mediator_mediator.set_matrices(mediator_addresses,
				mediator_addresses);

		// Restart the element
		const libMesh::Elem* elem_inter = micro_addresses.mesh.elem(0);
		inter_addresses.fe_unique_ptr->reinit(elem_inter);
		n_quadrature_pts = inter_addresses.fe_unique_ptr->n_quadrature_points();

		corrected_phi_BIG.resize(BIG_addresses.n_dofs_u,
				std::vector<libMesh::Real>(n_quadrature_pts, 0));
		corrected_phi_mediator.resize(mediator_addresses.n_dofs_u,
				std::vector<libMesh::Real>(n_quadrature_pts, 0));
		corrected_phi_micro.resize(micro_addresses.n_dofs_u,
				std::vector<libMesh::Real>(n_quadrature_pts, 0));

		if (m_bUseH1Coupling[micro_name])
		{
			corrected_dphi_BIG.resize(BIG_addresses.n_dofs_u,
					std::vector<libMesh::RealGradient>(n_quadrature_pts));
			corrected_dphi_mediator.resize(mediator_addresses.n_dofs_u,
					std::vector<libMesh::RealGradient>(n_quadrature_pts));
			corrected_dphi_micro.resize(micro_addresses.n_dofs_u,
					std::vector<libMesh::RealGradient>(n_quadrature_pts));
		}
	}

	// Vectors containing the lambda weights
	std::vector<std::vector<libMesh::Real> > lambda_weight_mediator(
			mediator_addresses.n_dofs_u,
			std::vector<libMesh::Real>(inter_addresses.n_dofs_u, 0));
	std::vector<std::vector<libMesh::Real> > lambda_weight_BIG(
			BIG_addresses.n_dofs_u,
			std::vector<libMesh::Real>(inter_addresses.n_dofs_u, 0));
	std::vector<std::vector<libMesh::Real> > lambda_weight_micro(
			micro_addresses.n_dofs_u,
			std::vector<libMesh::Real>(inter_addresses.n_dofs_u, 0));

	// Initialize global matrix
	const unsigned int mediator_M = mediator_addresses.dof_map.n_dofs();

	const unsigned int BIG_N = BIG_addresses.dof_map.n_dofs();
	const unsigned int micro_N = micro_addresses.dof_map.n_dofs();

	const unsigned int mediator_M_local =
			mediator_addresses.dof_map.n_local_dofs();

	const unsigned int BIG_N_local = BIG_addresses.dof_map.n_local_dofs();
	const unsigned int micro_N_local = micro_addresses.dof_map.n_local_dofs();

	couplingMatrix_mediator_micro.init(mediator_M, micro_N, mediator_M_local,
			micro_N_local, micro_N_local, micro_N - micro_N_local);
	couplingMatrix_mediator_BIG.init(mediator_M, BIG_N, mediator_M_local,
			BIG_N_local, BIG_N_local, BIG_N - BIG_N_local);

	couplingMatrix_mediator_mediator.attach_dof_map(mediator_addresses.dof_map);
	couplingMatrix_mediator_mediator.init();

	// Intersection indexes and iterators
	int nb_of_intersections = intersection_table_mediator_micro.size();
	int elem_mediator_idx = -1;

	int elem_BIG_idx = -1;
	int elem_micro_idx = -1;
	int elem_inter_idx = -1;
	std::unordered_multimap<int, int>::iterator it_inter_idx;

	// DEBUG CODE !!!
	int DEBUG_max_inter_idx = -1;
	int DEBUG_nb_of_inter_elems = 0;
	double DEBUG_vol = 0;

	// For each intersection
	for (int iii = 0; iii < nb_of_intersections; ++iii)
	{
		// Get the mediator and micro element pointers
		elem_mediator_idx = intersection_table_mediator_micro[iii].first;
		elem_micro_idx = intersection_table_mediator_micro[iii].second;

		if (using_same_mesh_mediator_A)
		{
			elem_BIG_idx = elem_mediator_idx;
		}
		else
		{
			elem_BIG_idx = equivalence_table_mediator_BIG[elem_mediator_idx];
		}

		const libMesh::Elem* elem_mediator = mediator_addresses.mesh.elem(
				elem_mediator_idx);
		const libMesh::Elem* elem_BIG = BIG_addresses.mesh.elem(elem_BIG_idx);
		const libMesh::Elem* elem_micro = micro_addresses.mesh.elem(
				elem_micro_idx);

		mediator_addresses.set_DoFs(elem_mediator_idx);
		BIG_addresses.set_DoFs(elem_BIG_idx);
		micro_addresses.set_DoFs(elem_micro_idx);

		// Resize dense matrix, if needed
		if (!bSameElemsType)
		{
			Me_mediator_micro.set_matrices(mediator_addresses, micro_addresses);
			Me_mediator_BIG.set_matrices(mediator_addresses, BIG_addresses);
			Me_mediator_mediator.set_matrices(mediator_addresses,
					mediator_addresses);

			// Set up corrected shape vectors
			const libMesh::Elem* elem_inter = micro_addresses.mesh.elem(0);
			inter_addresses.fe_unique_ptr->reinit(elem_inter);
			n_quadrature_pts =
					inter_addresses.fe_unique_ptr->n_quadrature_points();

			corrected_phi_BIG.resize(BIG_addresses.n_dofs_u,
					std::vector<libMesh::Real>(n_quadrature_pts, 0));
			corrected_phi_mediator.resize(mediator_addresses.n_dofs_u,
					std::vector<libMesh::Real>(n_quadrature_pts, 0));
			corrected_phi_micro.resize(micro_addresses.n_dofs_u,
					std::vector<libMesh::Real>(n_quadrature_pts, 0));

			if (m_bUseH1Coupling[micro_name])
			{
				corrected_dphi_BIG.resize(BIG_addresses.n_dofs_u,
						std::vector<libMesh::RealGradient>(n_quadrature_pts));
				corrected_dphi_mediator.resize(mediator_addresses.n_dofs_u,
						std::vector<libMesh::RealGradient>(n_quadrature_pts));
				corrected_dphi_micro.resize(micro_addresses.n_dofs_u,
						std::vector<libMesh::RealGradient>(n_quadrature_pts));
			}
		}

		Me_mediator_micro.zero();
		Me_mediator_BIG.zero();
		Me_mediator_mediator.zero();

		// Now iterate over the intersections
		auto inter_idx_range = intersection_table_inter.equal_range(iii);

		for (it_inter_idx = inter_idx_range.first;
				it_inter_idx != inter_idx_range.second; ++it_inter_idx)
		{
			// Get the intersection mesh pointer
			elem_inter_idx = it_inter_idx->second;

			// DEBUG!!!
			if (MASTER_debug_coupling_assemble)
			{
				if (elem_inter_idx > DEBUG_max_inter_idx)
				{
					DEBUG_max_inter_idx = elem_inter_idx;
				}
			}

			const libMesh::Elem* elem_inter = inter_addresses.mesh.elem(
					elem_inter_idx);

			// DEBUG !!!
			if (MASTER_debug_coupling_assemble)
			{
				DEBUG_vol += elem_inter->volume();
				++DEBUG_nb_of_inter_elems;
			}

			// Restart the element
			inter_addresses.fe_unique_ptr->reinit(elem_inter);

			get_lambdas(BIG_addresses.dim, BIG_addresses.fe_type, elem_BIG,
					quad_points_inter, quad_points_reference,
					lambda_weight_BIG);

			get_lambdas(micro_addresses.dim, micro_addresses.fe_type,
					elem_micro, quad_points_inter, quad_points_reference,
					lambda_weight_micro);

			get_lambdas(mediator_addresses.dim, mediator_addresses.fe_type,
					elem_mediator, quad_points_inter, quad_points_reference,
					lambda_weight_mediator);

			set_corrected_shapes(lambda_weight_BIG, phi_inter,
					corrected_phi_BIG);
			set_corrected_shapes(lambda_weight_micro, phi_inter,
					corrected_phi_micro);
			set_corrected_shapes(lambda_weight_mediator, phi_inter,
					corrected_phi_mediator);

			if (m_bUseH1Coupling[micro_name])
			{
				set_corrected_shape_gradients(lambda_weight_BIG, dphi_inter,
						corrected_dphi_BIG);
				set_corrected_shape_gradients(lambda_weight_micro, dphi_inter,
						corrected_dphi_micro);
				set_corrected_shape_gradients(lambda_weight_mediator,
						dphi_inter, corrected_dphi_mediator);
			}

			// For each quadrature point determinate the sub-matrices elements
			for (unsigned int qp = 0; qp < inter_addresses.qrule.n_points();
					qp++)
			{
				// Mediator -> micro coupling
				Me_mediator_micro.build_L2_coupling_matrix(mediator_addresses,
						micro_addresses, qp, corrected_phi_mediator,
						corrected_phi_micro, JxW, L2_coupling_const);

				// Mediator -> BIG coupling
				Me_mediator_BIG.build_L2_coupling_matrix(mediator_addresses,
						BIG_addresses, qp, corrected_phi_mediator,
						corrected_phi_BIG, JxW, L2_coupling_const);

				// Mediator -> Mediator coupling
				Me_mediator_mediator.build_L2_coupling_matrix(
						mediator_addresses, mediator_addresses, qp,
						corrected_phi_mediator, corrected_phi_mediator, JxW,
						L2_coupling_const);

				if (m_bUseH1Coupling[micro_name])
				{
					// Then we must also build the strain terms

					// Mediator -> micro coupling
					Me_mediator_micro.add_H1_coupling_matrix(mediator_addresses,
							micro_addresses, qp, corrected_dphi_mediator,
							corrected_dphi_micro, JxW, H1_coupling_const);

					Me_mediator_BIG.add_H1_coupling_matrix(mediator_addresses,
							BIG_addresses, qp, corrected_dphi_mediator,
							corrected_dphi_BIG, JxW, H1_coupling_const);

					Me_mediator_mediator.add_H1_coupling_matrix(
							mediator_addresses, mediator_addresses, qp,
							corrected_dphi_mediator, corrected_dphi_mediator,
							JxW, H1_coupling_const);
				}
			}
		}

		couplingMatrix_mediator_micro.add_matrix(Me_mediator_micro.Me,
				mediator_addresses.dof_indices, micro_addresses.dof_indices);
		couplingMatrix_mediator_BIG.add_matrix(Me_mediator_BIG.Me,
				mediator_addresses.dof_indices, BIG_addresses.dof_indices);
		couplingMatrix_mediator_mediator.add_matrix(Me_mediator_mediator.Me,
				mediator_addresses.dof_indices, mediator_addresses.dof_indices);
	}

	couplingMatrix_mediator_micro.close();
	couplingMatrix_mediator_BIG.close();
	couplingMatrix_mediator_mediator.close();

	if (MASTER_debug_coupling_assemble)
	{
		std::cout << "> COUPLING DEBUG !!! " << std::endl;
		std::cout << ">" << std::endl;
		std::cout << ">    nb_of_intersections     = " << nb_of_intersections
				<< std::endl;
		std::cout << ">    DEBUG_max_inter_idx     = " << DEBUG_max_inter_idx
				<< std::endl;
		std::cout << ">    DEBUG_nb_of_inter_elems = "
				<< DEBUG_nb_of_inter_elems << std::endl;
		std::cout << ">    DEBUG_vol               = " << DEBUG_vol << std::endl
				<< std::endl;
	}
};

void carl::coupled_system::assemble_coupling_elasticity_3D_parallel(
		const std::string BIG_name,
		const std::string micro_name,
		const std::string inter_name,
		const std::string mediator_name,

		const libMesh::MeshBase& mesh_R_BIG,
		const libMesh::MeshBase& mesh_R_micro,

		const std::unordered_map<int,std::pair<int,int> >&
										full_intersection_pairs_map,
		const std::unordered_map<int,std::pair<int,int> >&
										full_intersection_restricted_pairs_map,
		const std::unordered_map<int,int>&
										local_intersection_meshI_to_inter_map,
		const std::string BIG_type,
		const std::string micro_type,
		bool bSameElemsType)
{
	// TODO : make it possible to invert the algorithm's systems!
	// Addresses to the eq. systems
	libMesh::EquationSystems& mediator_eq_system =
			*m_mediator_EquationSystemMap[mediator_name];
	libMesh::EquationSystems& BIG_eq_system =
			*m_BIG_EquationSystem.second;
	libMesh::EquationSystems& micro_eq_system =
			*m_micro_EquationSystemMap[micro_name];
	libMesh::EquationSystems& inter_eq_system =
			*m_inter_EquationSystemMap[inter_name];

	libMesh::EquationSystems& R_BIG_eq_system =
			*m_R_BIG_EquationSystem.second;
	libMesh::EquationSystems& R_micro_eq_system =
			*m_R_micro_EquationSystemMap[micro_name];

	// First, test if all the systems have an elasticity model and variable set
	homemade_assert_msg(micro_eq_system.has_system(micro_type),
			" Micro equation systems is missing a system type!");
	homemade_assert_msg(BIG_eq_system.has_system(BIG_type),
			" Macro equation systems is missing a system type!");
	homemade_assert_msg(R_micro_eq_system.has_system("Elasticity"),
			" Restricted micro equation systems missing \"Elasticity\" system!");
	homemade_assert_msg(R_BIG_eq_system.has_system("Elasticity"),
			" Restricted macro equation systems missing \"Elasticity\" system!");
	homemade_assert_msg(inter_eq_system.has_system("Elasticity"),
			" Intersection equation systems missing \"Elasticity\" system!");
	homemade_assert_msg(mediator_eq_system.has_system("Elasticity"),
			" Mediatored equation systems missing \"Elasticity\" system!");

	// Systems and vars
	libMesh::ImplicitSystem& volume_mediator_system =
			libMesh::cast_ref<libMesh::ImplicitSystem&>(mediator_eq_system.get_system("Elasticity"));

	libMesh::ImplicitSystem& volume_BIG_system =
			libMesh::cast_ref<libMesh::ImplicitSystem&>(BIG_eq_system.get_system<libMesh::LinearImplicitSystem>(BIG_type));

	libMesh::ImplicitSystem& volume_micro_system =
			libMesh::cast_ref<libMesh::ImplicitSystem&>(micro_eq_system.get_system<libMesh::LinearImplicitSystem>(micro_type));

	libMesh::LinearImplicitSystem& volume_inter_system =
			inter_eq_system.get_system<libMesh::LinearImplicitSystem>("Elasticity");

	libMesh::LinearImplicitSystem& volume_R_BIG_system =
			R_BIG_eq_system.get_system<libMesh::LinearImplicitSystem>("Elasticity");

	libMesh::LinearImplicitSystem& volume_R_micro_system =
			R_micro_eq_system.get_system<libMesh::LinearImplicitSystem>("Elasticity");

	libMesh_fe_addresses_3 mediator_addresses(volume_mediator_system);
	libMesh_fe_addresses_3 BIG_addresses(volume_BIG_system);
	libMesh_fe_addresses_3 micro_addresses(volume_micro_system);
	libMesh_fe_addresses_3 inter_addresses(volume_inter_system);

	libMesh_fe_addresses_3 R_BIG_addresses(volume_R_BIG_system);
	libMesh_fe_addresses_3 R_micro_addresses(volume_R_micro_system);

	// Vector that will keep the quadrature points
	const std::vector<libMesh::Point>& quad_points_inter =
			inter_addresses.fe_unique_ptr->get_xyz();
	std::vector<libMesh::Point> quad_points_reference;

	// Jacobians
	const std::vector<libMesh::Real>& JxW =
			inter_addresses.fe_unique_ptr->get_JxW();

	// Shape functions
	const std::vector<std::vector<libMesh::Real> >& phi_inter =
			inter_addresses.fe_unique_ptr->get_phi();
	std::vector<std::vector<libMesh::Real> > corrected_phi_R_BIG;
	std::vector<std::vector<libMesh::Real> > corrected_phi_R_micro;
	std::vector<std::vector<libMesh::Real> > corrected_phi_mediator;

	// Shape functions gradients
	const std::vector<std::vector<libMesh::RealGradient> >& dphi_inter =
			inter_addresses.fe_unique_ptr->get_dphi();
	std::vector<std::vector<libMesh::RealGradient> > corrected_dphi_R_BIG;
	std::vector<std::vector<libMesh::RealGradient> > corrected_dphi_R_micro;
	std::vector<std::vector<libMesh::RealGradient> > corrected_dphi_mediator;

	unsigned int n_quadrature_pts = 0;

	// Addresses to the matrices
	libMesh::PetscMatrix<libMesh::Number>& couplingMatrix_mediator_BIG =
			*m_couplingMatrixMap_mediator_BIG[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& couplingMatrix_mediator_micro =
			*m_couplingMatrixMap_mediator_micro[micro_name];
	libMesh::PetscMatrix<libMesh::Number>& couplingMatrix_mediator_mediator =
			*m_couplingMatrixMap_mediator_mediator[micro_name];

	// Values of the coupling constants
	double L2_coupling_const = m_coupling_constantMap[micro_name]
			/ (m_coupling_lengthMap[micro_name]
					* m_coupling_lengthMap[micro_name]);
	double H1_coupling_const = m_coupling_constantMap[micro_name];

	// DoF vectors and ranges
	mediator_addresses.set_DoFs();
	BIG_addresses.set_DoFs();
	micro_addresses.set_DoFs();
	inter_addresses.set_DoFs();

	R_BIG_addresses.set_DoFs();
	R_micro_addresses.set_DoFs();

	// Local matrix
	coupling_matrices_3 Me_mediator_R_micro;
	coupling_matrices_3 Me_mediator_R_BIG;
	coupling_matrices_3 Me_mediator_mediator;

	//    If all elements are of the same type, do the index "extraction",
	// the matrices resizes and repositions here
	if (bSameElemsType)
	{
		Me_mediator_R_micro.set_matrices(mediator_addresses, R_micro_addresses);
		Me_mediator_R_BIG.set_matrices(mediator_addresses, R_BIG_addresses);
		Me_mediator_mediator.set_matrices(mediator_addresses,
				mediator_addresses);

		// Restart the element

		const libMesh::Elem* elem_inter =
				*inter_addresses.mesh.active_local_elements_begin();
		inter_addresses.fe_unique_ptr->reinit(elem_inter);
		n_quadrature_pts = inter_addresses.fe_unique_ptr->n_quadrature_points();

		corrected_phi_R_BIG.resize(R_BIG_addresses.n_dofs_u,
				std::vector<libMesh::Real>(n_quadrature_pts, 0));
		corrected_phi_R_micro.resize(R_micro_addresses.n_dofs_u,
				std::vector<libMesh::Real>(n_quadrature_pts, 0));
		corrected_phi_mediator.resize(mediator_addresses.n_dofs_u,
				std::vector<libMesh::Real>(n_quadrature_pts, 0));

		if (m_bUseH1Coupling[micro_name])
		{
			corrected_dphi_R_BIG.resize(R_BIG_addresses.n_dofs_u,
					std::vector<libMesh::RealGradient>(n_quadrature_pts));
			corrected_dphi_mediator.resize(mediator_addresses.n_dofs_u,
					std::vector<libMesh::RealGradient>(n_quadrature_pts));
			corrected_dphi_R_micro.resize(R_micro_addresses.n_dofs_u,
					std::vector<libMesh::RealGradient>(n_quadrature_pts));
		}
	}

	// Vectors containing the lambda weights
	std::vector<std::vector<libMesh::Real> > lambda_weight_mediator(
			mediator_addresses.n_dofs_u,
			std::vector<libMesh::Real>(inter_addresses.n_dofs_u, 0));
	std::vector<std::vector<libMesh::Real> > lambda_weight_R_BIG(
			R_BIG_addresses.n_dofs_u,
			std::vector<libMesh::Real>(inter_addresses.n_dofs_u, 0));
	std::vector<std::vector<libMesh::Real> > lambda_weight_R_micro(
			R_micro_addresses.n_dofs_u,
			std::vector<libMesh::Real>(inter_addresses.n_dofs_u, 0));

	// Initialize global matrix
	const unsigned int mediator_M = mediator_addresses.dof_map.n_dofs();

	const unsigned int BIG_N = BIG_addresses.dof_map.n_dofs();
	const unsigned int micro_N = micro_addresses.dof_map.n_dofs();

	const unsigned int mediator_M_local =
			mediator_addresses.dof_map.n_local_dofs();

	const unsigned int BIG_N_local = BIG_addresses.dof_map.n_local_dofs();
	const unsigned int micro_N_local = micro_addresses.dof_map.n_local_dofs();

	couplingMatrix_mediator_micro.init(mediator_M, micro_N, mediator_M_local,
			micro_N_local, micro_N_local, micro_N - micro_N_local);
	couplingMatrix_mediator_BIG.init(mediator_M, BIG_N, mediator_M_local,
			BIG_N_local, BIG_N_local, BIG_N - BIG_N_local);

	couplingMatrix_mediator_mediator.attach_dof_map(mediator_addresses.dof_map);
	couplingMatrix_mediator_mediator.init();

	// Intersection indexes and iterators
	int inter_idx = -1;
	int elem_restrict_BIG_idx = -1;
	int elem_restrict_micro_idx = -1;

	int elem_BIG_idx = -1;
	int elem_micro_idx = -1;
	int elem_mediator_idx = -1;
	int elem_inter_idx = -1;

	std::pair<int, int> restrict_idx_pair;
	std::pair<int, int> idx_pair;

	// DEBUG CODE !!!
	int DEBUG_max_inter_idx = -1;
	int DEBUG_nb_of_inter_elems = 0;
	double DEBUG_vol = 0;

	// For each intersection
	libMesh::MeshBase::const_element_iterator			inter_elIt =
			inter_addresses.mesh.active_local_elements_begin();
	const libMesh::MeshBase::const_element_iterator		end_inter_elIt =
			inter_addresses.mesh.active_local_elements_end();

	for( ; inter_elIt != end_inter_elIt; ++inter_elIt )
	{
		const libMesh::Elem* elem_inter = *inter_elIt;

		// Get the intersection element idx
		elem_inter_idx = elem_inter->id();
		inter_idx = local_intersection_meshI_to_inter_map.at(elem_inter_idx);

		restrict_idx_pair = full_intersection_restricted_pairs_map.at(inter_idx);
		elem_restrict_BIG_idx = restrict_idx_pair.first;
		elem_restrict_micro_idx = restrict_idx_pair.second;

		idx_pair = full_intersection_pairs_map.at(inter_idx);
		elem_BIG_idx = idx_pair.first;
		elem_micro_idx = idx_pair.second;

		// TODO For now, we suppose that BIG contains the mediator. So ...
		elem_mediator_idx = elem_restrict_BIG_idx;

		const libMesh::Elem* elem_mediator =
				mediator_addresses.mesh.elem(elem_mediator_idx);
		const libMesh::Elem* elem_R_BIG =
				R_BIG_addresses.mesh.elem(elem_restrict_BIG_idx);
		const libMesh::Elem* elem_R_micro =
				R_micro_addresses.mesh.elem(elem_restrict_micro_idx);

		mediator_addresses.set_DoFs(elem_mediator_idx);
		R_BIG_addresses.set_DoFs(elem_restrict_BIG_idx);
		R_micro_addresses.set_DoFs(elem_restrict_micro_idx);

		BIG_addresses.set_DoFs(elem_BIG_idx);
		micro_addresses.set_DoFs(elem_micro_idx);

		// Resize dense matrix, if needed
		if (!bSameElemsType)
		{
			Me_mediator_R_micro.set_matrices(mediator_addresses,
					R_micro_addresses);
			Me_mediator_R_BIG.set_matrices(mediator_addresses,
					R_BIG_addresses);
			Me_mediator_mediator.set_matrices(mediator_addresses,
					mediator_addresses);

			// Set up corrected shape vectors
			inter_addresses.fe_unique_ptr->reinit(elem_inter);
			n_quadrature_pts =
					inter_addresses.fe_unique_ptr->n_quadrature_points();

			corrected_phi_R_BIG.resize(R_BIG_addresses.n_dofs_u,
					std::vector<libMesh::Real>(n_quadrature_pts, 0));
			corrected_phi_mediator.resize(mediator_addresses.n_dofs_u,
					std::vector<libMesh::Real>(n_quadrature_pts, 0));
			corrected_phi_R_micro.resize(R_micro_addresses.n_dofs_u,
					std::vector<libMesh::Real>(n_quadrature_pts, 0));

			if (m_bUseH1Coupling[micro_name])
			{
				corrected_dphi_R_BIG.resize(R_BIG_addresses.n_dofs_u,
						std::vector<libMesh::RealGradient>(n_quadrature_pts));
				corrected_dphi_mediator.resize(mediator_addresses.n_dofs_u,
						std::vector<libMesh::RealGradient>(n_quadrature_pts));
				corrected_dphi_R_micro.resize(R_micro_addresses.n_dofs_u,
						std::vector<libMesh::RealGradient>(n_quadrature_pts));
			}
		}

		Me_mediator_R_micro.zero();
		Me_mediator_R_BIG.zero();
		Me_mediator_mediator.zero();

		// DEBUG!!!
		if (MASTER_debug_coupling_assemble)
		{
			if (elem_inter_idx > DEBUG_max_inter_idx)
			{
				DEBUG_max_inter_idx = elem_inter_idx;
			}
		}

		// DEBUG !!!
		if (MASTER_debug_coupling_assemble)
		{
			DEBUG_vol += elem_inter->volume();
			++DEBUG_nb_of_inter_elems;
		}

		// Restart the element
		inter_addresses.fe_unique_ptr->reinit(elem_inter);

		get_lambdas(R_BIG_addresses.dim, R_BIG_addresses.fe_type, elem_R_BIG,
				quad_points_inter, quad_points_reference,
				lambda_weight_R_BIG);

		get_lambdas(R_micro_addresses.dim, R_micro_addresses.fe_type,
				elem_R_micro, quad_points_inter, quad_points_reference,
				lambda_weight_R_micro);

		get_lambdas(mediator_addresses.dim, mediator_addresses.fe_type,
				elem_mediator, quad_points_inter, quad_points_reference,
				lambda_weight_mediator);

		set_corrected_shapes(lambda_weight_R_BIG, phi_inter,
				corrected_phi_R_BIG);
		set_corrected_shapes(lambda_weight_R_micro, phi_inter,
				corrected_phi_R_micro);
		set_corrected_shapes(lambda_weight_mediator, phi_inter,
				corrected_phi_mediator);

		if (m_bUseH1Coupling[micro_name])
		{
			set_corrected_shape_gradients(lambda_weight_R_BIG, dphi_inter,
					corrected_dphi_R_BIG);
			set_corrected_shape_gradients(lambda_weight_R_micro, dphi_inter,
					corrected_dphi_R_micro);
			set_corrected_shape_gradients(lambda_weight_mediator,
					dphi_inter, corrected_dphi_mediator);
		}

		// For each quadrature point determinate the sub-matrices elements
		for (unsigned int qp = 0; qp < inter_addresses.qrule.n_points();
				qp++)
		{
			// Mediator -> micro coupling
			Me_mediator_R_micro.build_L2_coupling_matrix(mediator_addresses,
					R_micro_addresses, qp, corrected_phi_mediator,
					corrected_phi_R_micro, JxW, L2_coupling_const);

			// Mediator -> BIG coupling
			Me_mediator_R_BIG.build_L2_coupling_matrix(mediator_addresses,
					R_BIG_addresses, qp, corrected_phi_mediator,
					corrected_phi_R_BIG, JxW, L2_coupling_const);

			// Mediator -> Mediator coupling
			Me_mediator_mediator.build_L2_coupling_matrix(
					mediator_addresses, mediator_addresses, qp,
					corrected_phi_mediator, corrected_phi_mediator, JxW,
					L2_coupling_const);

			if (m_bUseH1Coupling[micro_name])
			{
				// Then we must also build the strain terms

				// Mediator -> micro coupling
				Me_mediator_R_micro.add_H1_coupling_matrix(mediator_addresses,
						R_micro_addresses, qp, corrected_dphi_mediator,
						corrected_dphi_R_micro, JxW, H1_coupling_const);

				Me_mediator_R_BIG.add_H1_coupling_matrix(mediator_addresses,
						R_BIG_addresses, qp, corrected_dphi_mediator,
						corrected_dphi_R_BIG, JxW, H1_coupling_const);

				Me_mediator_mediator.add_H1_coupling_matrix(
						mediator_addresses, mediator_addresses, qp,
						corrected_dphi_mediator, corrected_dphi_mediator,
						JxW, H1_coupling_const);
			}
		}

		couplingMatrix_mediator_micro.add_matrix(Me_mediator_R_micro.Me,
				mediator_addresses.dof_indices, micro_addresses.dof_indices);
		couplingMatrix_mediator_BIG.add_matrix(Me_mediator_R_BIG.Me,
				mediator_addresses.dof_indices, BIG_addresses.dof_indices);
		couplingMatrix_mediator_mediator.add_matrix(Me_mediator_mediator.Me,
				mediator_addresses.dof_indices, mediator_addresses.dof_indices);
	}

	const libMesh::Parallel::Communicator& WorldComm = couplingMatrix_mediator_micro.comm();
	WorldComm.barrier();

	couplingMatrix_mediator_micro.close();
	couplingMatrix_mediator_BIG.close();
	couplingMatrix_mediator_mediator.close();

	if (MASTER_debug_coupling_assemble)
	{
		std::cout << "> COUPLING DEBUG !!! " << std::endl;
		std::cout << ">" << std::endl;
		std::cout << ">    DEBUG_max_inter_idx     = " << DEBUG_max_inter_idx
				<< std::endl;
		std::cout << ">    DEBUG_nb_of_inter_elems = "
				<< DEBUG_nb_of_inter_elems << std::endl;
		std::cout << ">    DEBUG_vol               = " << DEBUG_vol << std::endl
				<< std::endl;
	}
};
