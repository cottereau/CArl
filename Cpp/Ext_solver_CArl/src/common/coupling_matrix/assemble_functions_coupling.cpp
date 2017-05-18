#include "assemble_coupling.h"

namespace carl
{

	void assemble_coupling_matrices::prepare_coupling_preallocation(
			libMesh::PetscMatrix<libMesh::Number>& coupling_matrix,
			libMesh_fe_addresses_3& row_addresses,
			libMesh_fe_addresses_3& col_addresses,
			const std::unordered_multimap<int,int>&  inter_table
			)
{
	const libMesh::Parallel::Communicator& WorldComm = row_addresses.mesh.comm();
	int rank = WorldComm.rank();
	int nodes = WorldComm.size();

	// Local and global dimensions setup
	const unsigned int row_M = row_addresses.dof_map.n_dofs();
	const int col_N = col_addresses.dof_map.n_dofs();
	const unsigned int row_M_local = row_addresses.dof_map.n_local_dofs();
	const int col_N_local = col_addresses.dof_map.n_local_dofs();

	// Vectors that will be used to calculate the preallocation. There are two
	// "versions": a local one, starting at zero and which will contain this
	// processor's preallocation after syncing, and a global one, which will be
	// send to other processors
	std::vector<unsigned int> local_n_nz(row_M_local,0);
	std::vector<unsigned int> local_n_oz(row_M_local,0);

	// Local limits of the DoF's
	std::vector<int> begin_row_DoF(nodes);
	std::vector<int> end_row_DoF(nodes);
	for(int iii = 0; iii < nodes; ++iii)
	{
		begin_row_DoF[iii] = row_addresses.dof_map.first_dof(iii);
		end_row_DoF[iii] = row_addresses.dof_map.end_dof(iii);
	}

	std::vector<int> begin_col_DoF(nodes);
	std::vector<int> end_col_DoF(nodes);
	for(int iii = 0; iii < nodes; ++iii)
	{
		begin_col_DoF[iii] = col_addresses.dof_map.first_dof(iii);
		end_col_DoF[iii] = col_addresses.dof_map.end_dof(iii);
	}

	int local_begin_row_DoF = begin_row_DoF[rank];
	int local_end_row_DoF = end_row_DoF[rank];

	// Iterate over the row's mesh cells

	int row_mesh_elem_idx = 0;
	int dummy_col_dof;

	std::vector<std::vector<int> > to_add_per_proc(row_M, std::vector<int>(nodes,0));

	std::vector<int> row_dof_indices_procs;
	libMesh::MeshBase::const_element_iterator it_row_elem =
			row_addresses.mesh.active_local_elements_begin();
	const libMesh::MeshBase::const_element_iterator		it_row_elem_end =
			row_addresses.mesh.active_local_elements_end();

	// What the code has to do:
	/*
	 * 		- For each element associated to the rows, get the intersecting
	 * 		  elements and find out which processors own each of its DoFs. Use
	 * 		  this to populate the "to_add_per_proc" table.
	 *
	 * 		- Once it is filled, sync the "to_add_per_proc" between the
	 * 		  processors.
	 *
	 * 		- On each processor, iterate over the local row DoFs and add the
	 * 		  info from "to_add_per_proc" to either "local_n_nz" or "local_n_oz",
	 * 		  depending on the processor.
	 *
	 * 		- ... pray that it works?
	 */

	std::unordered_set<int> col_DoFs_to_allocate;
	col_DoFs_to_allocate.reserve(col_addresses.mesh.n_nodes());

	for( ; it_row_elem != it_row_elem_end; ++it_row_elem )
	{
		// For each local mediator mesh elements ...
		const libMesh::Elem* elem_mediator = *it_row_elem;

		// ... extract the element's index ...
		row_mesh_elem_idx = elem_mediator->id();

		// ... set its DoFs ...
		row_addresses.set_DoFs(row_mesh_elem_idx);

		// ... clear the col_DoFs set
		col_DoFs_to_allocate.clear();

		// ... get the DoFs from the col mesh ...
		auto range_elem_col_idx = inter_table.equal_range(row_mesh_elem_idx);

		for(auto it = range_elem_col_idx.first; it != range_elem_col_idx.second; ++it)
		{
			col_addresses.set_DoFs(it->second);

			for(unsigned int iii = 0; iii < col_addresses.n_dofs; ++iii)
			{
				col_DoFs_to_allocate.insert(col_addresses.dof_indices[iii]);
			}
		}

		// ... and finally increase the preallocation
		for(auto it = col_DoFs_to_allocate.begin(); it != col_DoFs_to_allocate.end(); ++it)
		{
			dummy_col_dof = *it;

			for(int nnn = 0; nnn < nodes; ++nnn)
			{
				if(	dummy_col_dof > begin_col_DoF[nnn] &&
					dummy_col_dof < end_col_DoF[nnn])
				{
					for(unsigned int jjj = 0; jjj < row_addresses.dof_indices.size(); ++jjj)
					{
						++to_add_per_proc[row_addresses.dof_indices[jjj]][nnn];
					}
					break;
				}
			}
		}
	}

	// Guarantee that everyone is in the same page
	WorldComm.barrier();

	// Do a (lot) of syncs
	for(int nnn = 0; nnn < nodes; ++nnn)
	{
		for(int iii = begin_row_DoF[nnn]; iii < end_row_DoF[nnn]; ++iii)
		{
			MPI_reduce_vector(to_add_per_proc[iii],nnn,WorldComm);
		}
	}

	int dummy_total_nz = 0;
	int dummy_total_oz = 0;

	// Now, calculate the preallocation
	for(int iii = local_begin_row_DoF; iii < local_end_row_DoF; ++iii)
	{
		for(int nnn = 0; nnn < rank; ++nnn)
		{
			local_n_oz[iii - local_begin_row_DoF] += to_add_per_proc[iii][nnn];
			dummy_total_oz += to_add_per_proc[iii][nnn];
		}

		local_n_nz[iii - local_begin_row_DoF] += to_add_per_proc[iii][rank];
		dummy_total_nz += to_add_per_proc[iii][rank];

		for(int nnn = rank + 1; nnn < nodes; ++nnn)
		{
			local_n_oz[iii - local_begin_row_DoF] += to_add_per_proc[iii][nnn];
			dummy_total_oz += to_add_per_proc[iii][nnn];
		}
	}

	dummy_total_nz = 0;
	dummy_total_oz = 0;

	for(int iii = local_begin_row_DoF; iii < local_end_row_DoF; ++iii)
	{
		local_n_nz[iii - local_begin_row_DoF] = 1.1*local_n_nz[iii - local_begin_row_DoF];
		local_n_nz[iii - local_begin_row_DoF] = std::max(static_cast<unsigned int>(30),local_n_nz[iii - local_begin_row_DoF]);
		local_n_nz[iii - local_begin_row_DoF] = std::min(static_cast<unsigned int>(col_N_local),local_n_nz[iii - local_begin_row_DoF]);

		local_n_oz[iii - local_begin_row_DoF] = 1.1*local_n_oz[iii - local_begin_row_DoF];
		local_n_oz[iii - local_begin_row_DoF] = std::max(static_cast<unsigned int>(10),local_n_oz[iii - local_begin_row_DoF]);
		local_n_oz[iii - local_begin_row_DoF] = std::min(static_cast<unsigned int>(col_N - col_N_local),local_n_oz[iii - local_begin_row_DoF]);

		dummy_total_nz += local_n_nz[iii - local_begin_row_DoF];
		dummy_total_oz += local_n_oz[iii - local_begin_row_DoF];
	}

	coupling_matrix.init(row_M, col_N, row_M_local, col_N_local, local_n_nz, local_n_oz);
}

	void assemble_coupling_matrices::assemble_coupling_elasticity_3D_parallel(
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
			const std::unordered_multimap<int,int>& inter_table_mediator_BIG,
			const std::unordered_multimap<int,int>& inter_table_mediator_micro,

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
		libMesh::System& volume_mediator_system =
				libMesh::cast_ref<libMesh::System&>(mediator_eq_system.get_system("Elasticity"));

		libMesh::System& volume_BIG_system =
				libMesh::cast_ref<libMesh::System&>(BIG_eq_system.get_system<libMesh::ExplicitSystem>(BIG_type));

		libMesh::System& volume_micro_system =
				libMesh::cast_ref<libMesh::System&>(micro_eq_system.get_system<libMesh::ExplicitSystem>(micro_type));

		libMesh::System& volume_inter_system =
				inter_eq_system.get_system<libMesh::ExplicitSystem>("Elasticity");

		libMesh::System& volume_R_BIG_system =
				R_BIG_eq_system.get_system<libMesh::ExplicitSystem>("Elasticity");

		libMesh::System& volume_R_micro_system =
				R_micro_eq_system.get_system<libMesh::ExplicitSystem>("Elasticity");

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

		const libMesh::Parallel::Communicator& WorldComm = couplingMatrix_mediator_micro.comm();

		// Values of the coupling constants
		double L2_coupling_const = m_coupling_rigidityMap[micro_name]
				/ (m_coupling_widthMap[micro_name]
						* m_coupling_widthMap[micro_name]);
		double H1_coupling_const = m_coupling_rigidityMap[micro_name];

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
		prepare_coupling_preallocation(
				couplingMatrix_mediator_BIG,
				mediator_addresses,
				BIG_addresses,
				inter_table_mediator_BIG
				);

		prepare_coupling_preallocation(
				couplingMatrix_mediator_micro,
				mediator_addresses,
				micro_addresses,
				inter_table_mediator_micro
				);
		
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

		WorldComm.barrier();

		couplingMatrix_mediator_micro.close();
		couplingMatrix_mediator_BIG.close();
		couplingMatrix_mediator_mediator.close();

		print_matrix_dim(couplingMatrix_mediator_micro);
		print_matrix_dim(couplingMatrix_mediator_BIG);
		print_matrix_dim(couplingMatrix_mediator_mediator);

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

	void assemble_coupling_matrices::check_coupling_construction_3D_parallel(
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
			const std::unordered_multimap<int,int>& inter_table_mediator_BIG,
			const std::unordered_multimap<int,int>& inter_table_mediator_micro,

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
		libMesh::System& volume_mediator_system =
				libMesh::cast_ref<libMesh::System&>(mediator_eq_system.get_system("Elasticity"));

		libMesh::System& volume_BIG_system =
				libMesh::cast_ref<libMesh::System&>(BIG_eq_system.get_system<libMesh::ExplicitSystem>(BIG_type));

		libMesh::System& volume_micro_system =
				libMesh::cast_ref<libMesh::System&>(micro_eq_system.get_system<libMesh::ExplicitSystem>(micro_type));

		libMesh::System& volume_inter_system =
				inter_eq_system.get_system<libMesh::ExplicitSystem>("Elasticity");

		libMesh::System& volume_R_BIG_system =
				R_BIG_eq_system.get_system<libMesh::ExplicitSystem>("Elasticity");

		libMesh::System& volume_R_micro_system =
				R_micro_eq_system.get_system<libMesh::ExplicitSystem>("Elasticity");

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
		libMesh::PetscMatrix<libMesh::Number> couplingMatrix_mediator_BIG(mesh_R_BIG.comm());
		libMesh::PetscMatrix<libMesh::Number> couplingMatrix_mediator_micro(mesh_R_BIG.comm());
		libMesh::PetscMatrix<libMesh::Number> couplingMatrix_mediator_mediator(mesh_R_BIG.comm());

		const libMesh::Parallel::Communicator& WorldComm = couplingMatrix_mediator_micro.comm();

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


				corrected_dphi_R_BIG.resize(R_BIG_addresses.n_dofs_u,
						std::vector<libMesh::RealGradient>(n_quadrature_pts));
				corrected_dphi_mediator.resize(mediator_addresses.n_dofs_u,
						std::vector<libMesh::RealGradient>(n_quadrature_pts));
				corrected_dphi_R_micro.resize(R_micro_addresses.n_dofs_u,
						std::vector<libMesh::RealGradient>(n_quadrature_pts));

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
		prepare_coupling_preallocation(
				couplingMatrix_mediator_BIG,
				mediator_addresses,
				BIG_addresses,
				inter_table_mediator_BIG
				);

		prepare_coupling_preallocation(
				couplingMatrix_mediator_micro,
				mediator_addresses,
				micro_addresses,
				inter_table_mediator_micro
				);

		couplingMatrix_mediator_mediator.attach_dof_map(mediator_addresses.dof_map);
		couplingMatrix_mediator_mediator.init();

		std::cout << " ----------------- " << std::endl;
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


					corrected_dphi_R_BIG.resize(R_BIG_addresses.n_dofs_u,
							std::vector<libMesh::RealGradient>(n_quadrature_pts));
					corrected_dphi_mediator.resize(mediator_addresses.n_dofs_u,
							std::vector<libMesh::RealGradient>(n_quadrature_pts));
					corrected_dphi_R_micro.resize(R_micro_addresses.n_dofs_u,
							std::vector<libMesh::RealGradient>(n_quadrature_pts));
			}

			Me_mediator_R_micro.zero();
			Me_mediator_R_BIG.zero();
			Me_mediator_mediator.zero();

			// DEBUG!!!

				if (elem_inter_idx > DEBUG_max_inter_idx)
				{
					DEBUG_max_inter_idx = elem_inter_idx;
				}


			// DEBUG !!!

				DEBUG_vol += std::abs(elem_inter->volume());
				++DEBUG_nb_of_inter_elems;


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

				set_corrected_shape_gradients(lambda_weight_R_BIG, dphi_inter,
						corrected_dphi_R_BIG);
				set_corrected_shape_gradients(lambda_weight_R_micro, dphi_inter,
						corrected_dphi_R_micro);
				set_corrected_shape_gradients(lambda_weight_mediator,
						dphi_inter, corrected_dphi_mediator);


			// For each quadrature point determinate the sub-matrices elements
			for (unsigned int qp = 0; qp < inter_addresses.qrule.n_points();
					qp++)
			{
				// Mediator -> micro coupling
				Me_mediator_R_micro.build_L2_coupling_matrix(mediator_addresses,
						R_micro_addresses, qp, corrected_phi_mediator,
						corrected_phi_R_micro, JxW, 1);

				// Mediator -> BIG coupling
				Me_mediator_R_BIG.build_L2_coupling_matrix(mediator_addresses,
						R_BIG_addresses, qp, corrected_phi_mediator,
						corrected_phi_R_BIG, JxW, 1);

				// Mediator -> Mediator coupling
				Me_mediator_mediator.build_L2_coupling_matrix(
						mediator_addresses, mediator_addresses, qp,
						corrected_phi_mediator, corrected_phi_mediator, JxW,
						1);

					// Then we must also build the strain terms

					// Mediator -> micro coupling
					Me_mediator_R_micro.add_H1_coupling_matrix(mediator_addresses,
							R_micro_addresses, qp, corrected_dphi_mediator,
							corrected_dphi_R_micro, JxW, 0);

					Me_mediator_R_BIG.add_H1_coupling_matrix(mediator_addresses,
							R_BIG_addresses, qp, corrected_dphi_mediator,
							corrected_dphi_R_BIG, JxW, 0);

					Me_mediator_mediator.add_H1_coupling_matrix(
							mediator_addresses, mediator_addresses, qp,
							corrected_dphi_mediator, corrected_dphi_mediator,
							JxW, 0);
			}

			couplingMatrix_mediator_micro.add_matrix(Me_mediator_R_micro.Me,
					mediator_addresses.dof_indices, micro_addresses.dof_indices);
			couplingMatrix_mediator_BIG.add_matrix(Me_mediator_R_BIG.Me,
					mediator_addresses.dof_indices, BIG_addresses.dof_indices);
			couplingMatrix_mediator_mediator.add_matrix(Me_mediator_mediator.Me,
					mediator_addresses.dof_indices, mediator_addresses.dof_indices);
		}

		WorldComm.barrier();

		couplingMatrix_mediator_micro.close();
		couplingMatrix_mediator_BIG.close();
		couplingMatrix_mediator_mediator.close();

		std::cout << " ----------------- " << std::endl;

		Vec row_sum_micro, row_sum_BIG, row_sum_mediator;
		MatCreateVecs(couplingMatrix_mediator_micro.mat(),&row_sum_micro,NULL);
		MatCreateVecs(couplingMatrix_mediator_BIG.mat(),&row_sum_BIG,NULL);
		MatCreateVecs(couplingMatrix_mediator_mediator.mat(),&row_sum_mediator,NULL);

		MatGetRowSum(couplingMatrix_mediator_micro.mat(),row_sum_micro);
		MatGetRowSum(couplingMatrix_mediator_BIG.mat(),row_sum_BIG);
		MatGetRowSum(couplingMatrix_mediator_mediator.mat(),row_sum_mediator);

		PetscScalar coupl_vol_micro = 0, coupl_vol_BIG = 0, coupl_vol_mediator = 0;

		VecSum(row_sum_micro,&coupl_vol_micro);
		VecSum(row_sum_BIG,&coupl_vol_BIG);
		VecSum(row_sum_mediator,&coupl_vol_mediator);

		mesh_R_BIG.comm().sum(DEBUG_vol);

			std::cout << "> COUPLING DEBUG !!! " << std::endl;
			std::cout << ">" << std::endl;
			std::cout << ">    DEBUG_max_inter_idx     = " << DEBUG_max_inter_idx
					<< std::endl;
			std::cout << ">    DEBUG_nb_of_inter_elems = "
					<< DEBUG_nb_of_inter_elems << std::endl;
			std::cout << ">    DEBUG_mesh_vol              = " << DEBUG_vol
					<< std::endl;
			std::cout << ">    DEBUG_micro_vol              = " << coupl_vol_micro << std::endl;
			std::cout << ">    DEBUG_big_vol              = " << coupl_vol_BIG << std::endl;
			std::cout << ">    DEBUG_mediator_vol              = " << coupl_vol_mediator << std::endl
					<< std::endl;

			VecDestroy(&row_sum_micro);
			VecDestroy(&row_sum_BIG);
			VecDestroy(&row_sum_mediator);
	};
};