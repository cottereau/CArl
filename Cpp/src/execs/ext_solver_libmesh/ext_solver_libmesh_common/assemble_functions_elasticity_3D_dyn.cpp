#include "assemble_functions_elasticity_3D.h"
#include "assemble_functions_elasticity_3D_dyn.h"

void set_homogeneous_physical_properties_dyn(libMesh::EquationSystems& es, std::string& physicalParamsFile)
{
	const libMesh::Parallel::Communicator& SysComm = es.comm();
	int rank = SysComm.rank();
	int nodes = SysComm.size();

	// Read the random data info
	double inputE;
	double inputMu;
    double inputRho;

	if(rank == 0)
	{
		std::ifstream physicalParamsIFS(physicalParamsFile);
		physicalParamsIFS >> inputE >> inputMu >> inputRho;
        std::cout << " >> DYNAMIC PROPERTIES: " << inputE<<inputMu<<inputRho; 
	}

	if(nodes > 1)
	{
		SysComm.broadcast(inputE);
		SysComm.broadcast(inputMu);
		SysComm.broadcast(inputRho);
	}

	// Mesh pointer
	const libMesh::MeshBase& mesh = es.get_mesh();

	// Physical system and its "variables"
	libMesh::ExplicitSystem& physical_param_system = es.get_system<libMesh::ExplicitSystem>("PhysicalConstants");
	const libMesh::DofMap& physical_dof_map = physical_param_system.get_dof_map();

	unsigned int physical_consts[3];
	physical_consts[0] = physical_param_system.variable_number ("E");
	physical_consts[1] = physical_param_system.variable_number ("mu");
	physical_consts[2] = physical_param_system.variable_number ("Rho");

	std::vector<libMesh::dof_id_type> physical_dof_indices_var;
	libMesh::MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
	const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	for ( ; el != end_el; ++el)
	{
		const libMesh::Elem* elem = *el;

		// Young modulus, E
		physical_dof_map.dof_indices(elem, physical_dof_indices_var, physical_consts[0]);
		libMesh::dof_id_type dof_index = physical_dof_indices_var[0];

		if( (physical_param_system.solution->first_local_index() <= dof_index) &&
		(dof_index < physical_param_system.solution->last_local_index()) )
		{
			physical_param_system.solution->set(dof_index, inputE);
		}

		// Mu
		physical_dof_map.dof_indices (elem, physical_dof_indices_var, physical_consts[1]);

		dof_index = physical_dof_indices_var[0];

		if( (physical_param_system.solution->first_local_index() <= dof_index) &&
		(dof_index < physical_param_system.solution->last_local_index()) )
		{
			physical_param_system.solution->set(dof_index, inputMu);
		}
        
		// Rho 
		physical_dof_map.dof_indices (elem, physical_dof_indices_var, physical_consts[2]);

		dof_index = physical_dof_indices_var[1];

		if( (physical_param_system.solution->first_local_index() <= dof_index) &&
		(dof_index < physical_param_system.solution->last_local_index()) )
		{
			physical_param_system.solution->set(dof_index, inputRho);
		}
	}

	physical_param_system.solution->close();
	physical_param_system.update();
}

void set_heterogeneous_physical_properties_dyn(libMesh::EquationSystems& es, std::string& physicalParamsFile)
{
	const libMesh::Parallel::Communicator& SysComm = es.comm();
	int rank = SysComm.rank();
	int nodes = SysComm.size();

	// Read the random data info
	std::vector<double> inputE;
	std::vector<double> inputMu;
	std::vector<double> inputRho;
	std::vector<int> 	inputIdx;

	int NbOfSubdomains = -1;

	double meanE = 0;
	double meanMu = 0;
	double meanRho = 0;

	if(rank == 0)
	{
		std::ifstream physicalParamsIFS(physicalParamsFile);
		physicalParamsIFS >> NbOfSubdomains;
		inputE.resize(NbOfSubdomains);
		inputMu.resize(NbOfSubdomains);
		inputRho.resize(NbOfSubdomains);
		inputIdx.resize(NbOfSubdomains);

		for(int iii = 0; iii < NbOfSubdomains; ++iii)
		{
			physicalParamsIFS >> inputE[iii];
			physicalParamsIFS >> inputMu[iii];
			physicalParamsIFS >> inputRho[iii];
			physicalParamsIFS >> inputIdx[iii];

			meanE += inputE[iii];
			meanMu += inputMu[iii];
			meanRho += inputRho[iii];
		}
		meanE /= NbOfSubdomains;
		meanMu /= NbOfSubdomains;
		meanRho /= NbOfSubdomains;
		physicalParamsIFS.close();
	}

	if(nodes > 1)
	{
		SysComm.broadcast(NbOfSubdomains);
		SysComm.broadcast(meanE);
		SysComm.broadcast(meanMu);
		SysComm.broadcast(meanRho);

		if(rank != 0)
		{
			inputE.resize(NbOfSubdomains);
			inputMu.resize(NbOfSubdomains);
			inputRho.resize(NbOfSubdomains);
			inputIdx.resize(NbOfSubdomains);
		}
		SysComm.broadcast(inputE);
		SysComm.broadcast(inputMu);
		SysComm.broadcast(inputRho);
		SysComm.broadcast(inputIdx);
	}

	// Mesh pointer
	const libMesh::MeshBase& mesh = es.get_mesh();

	// Physical system and its "variables"
	libMesh::ExplicitSystem& physical_param_system = es.get_system<libMesh::ExplicitSystem>("PhysicalConstants");
	const libMesh::DofMap& physical_dof_map = physical_param_system.get_dof_map();

	unsigned int physical_consts[3];
	physical_consts[0] = physical_param_system.variable_number ("E");
	physical_consts[1] = physical_param_system.variable_number ("mu");
	physical_consts[2] = physical_param_system.variable_number ("rho");

	std::vector<libMesh::dof_id_type> physical_dof_indices_var;
	libMesh::MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
	const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	int currentSubdomain = -1;

	for ( ; el != end_el; ++el)
	{
		const libMesh::Elem* elem = *el;

		currentSubdomain = elem->subdomain_id()-1;

		// Young modulus, E
		physical_dof_map.dof_indices(elem, physical_dof_indices_var, physical_consts[0]);
		libMesh::dof_id_type dof_index = physical_dof_indices_var[0];

		if( (physical_param_system.solution->first_local_index() <= dof_index) &&
		(dof_index < physical_param_system.solution->last_local_index()) )
		{
			physical_param_system.solution->set(dof_index, inputE[currentSubdomain]);
		}

		// Mu
		physical_dof_map.dof_indices (elem, physical_dof_indices_var, physical_consts[1]);

		dof_index = physical_dof_indices_var[0];

		if( (physical_param_system.solution->first_local_index() <= dof_index) &&
		(dof_index < physical_param_system.solution->last_local_index()) )
		{
			physical_param_system.solution->set(dof_index, inputMu[currentSubdomain]);
		}
		// Rho
		physical_dof_map.dof_indices (elem, physical_dof_indices_var, physical_consts[2]);

		dof_index = physical_dof_indices_var[1];

		if( (physical_param_system.solution->first_local_index() <= dof_index) &&
		(dof_index < physical_param_system.solution->last_local_index()) )
		{
			physical_param_system.solution->set(dof_index, inputRho[currentSubdomain]);
		}
	}

	physical_param_system.solution->close();
	physical_param_system.update();
}

// Update Mass Matrix
void Update_SubM_isotropic(	libMesh::DenseSubMatrix<libMesh::Number>& SubM,
        unsigned int qp,
        const std::vector<std::vector<libMesh::Real> >& phi,
        const unsigned int n_u_dofs,
        const std::vector<libMesh::Real>& JxW,
        libMesh::Number rho,
        double cte
        )

{
    for (unsigned int iii=0; iii<n_u_dofs; iii++)
    {
        for (unsigned int jjj=0; jjj<n_u_dofs; jjj++)
        {
            SubM(iii,jjj) += cte * JxW[qp]*rho*phi[iii][qp]*phi[jjj][qp];
        }
    }
};

// Update Damping Matrix
void Update_SubC_isotropic(	libMesh::DenseSubMatrix<libMesh::Number>& SubC,
        libMesh::DenseSubMatrix<libMesh::Number>& SubM, 
        libMesh::DenseSubMatrix<libMesh::Number>& SubK, 
        unsigned int qp,
        const unsigned int n_u_dofs,
        double cte,
        double alpha, 
        double beta
        )
{
    for (unsigned int iii=0; iii<n_u_dofs; iii++)
    {
        for (unsigned int jjj=0; jjj<n_u_dofs; jjj++)
        {
            SubC(iii,jjj) = cte*(alpha*SubM(iii,jjj)+beta*SubK(iii,jjj));
        }
    }
};


void assemble_elasticity_with_weight_dyn(	libMesh::EquationSystems& es,
							const std::string& system_name,
							weight_parameter_function& weight_mask,
							WeightFunctionSystemType system_type)
{
	libmesh_assert_equal_to(system_name, "Transient");

	libMesh::PerfLog perf_log ("Matrix Assembly (Homogeneous elasticity)",MASTER_bPerfLog_assemble_fem);

	perf_log.push("Preamble");

	const libMesh::MeshBase& mesh = es.get_mesh();

	const unsigned int dim = mesh.mesh_dimension();

	// - Set up physical properties system ------------------------------------
	libMesh::ExplicitSystem& physical_param_system = es.get_system<libMesh::ExplicitSystem>("PhysicalConstants");
	const unsigned int young_var = physical_param_system.variable_number ("E");
	const unsigned int mu_var = physical_param_system.variable_number ("mu");
	const unsigned int rho_var = physical_param_system.variable_number ("rho");

	const libMesh::DofMap& physical_dof_map = physical_param_system.get_dof_map();
	std::vector<libMesh::dof_id_type> physical_dof_indices_var;

	// The DoF and values of the physical system
	const libMesh::Elem*  phys_elem   = *(mesh.active_local_elements_begin());
	physical_dof_map.dof_indices(phys_elem, physical_dof_indices_var, young_var);
	libMesh::Number localE = physical_param_system.current_solution(physical_dof_indices_var[0]);

	physical_dof_map.dof_indices(phys_elem, physical_dof_indices_var, mu_var);
	libMesh::Number localMu = physical_param_system.current_solution(physical_dof_indices_var[0]);

	physical_dof_map.dof_indices(phys_elem, physical_dof_indices_var, rho_var);
	libMesh::Number localRho = physical_param_system.current_solution(physical_dof_indices_var[0]);
	// - Set up elasticity system ---------------------------------------------
	//libMesh::LinearImplicitSystem& system = es.get_system<libMesh::LinearImplicitSystem>("Elasticity");
    //libMesh::LinearImplicitSystem& system = es.get_system<libMesh::LinearImplicitSystem>("Elasticity");
    libMesh::NewmarkSystem& system = es.get_system<libMesh::NewmarkSystem>(system_name);	

    const unsigned int n_components = 3;
	const unsigned int u_var = system.variable_number ("u");
	const unsigned int v_var = system.variable_number ("v");
	const unsigned int w_var = system.variable_number ("w");

	const libMesh::DofMap& dof_map = system.get_dof_map();
	libMesh::FEType fe_type = dof_map.variable_type(u_var);

    // 
    libMesh::SparseMatrix<libMesh::Number>&  stiffness = system.get_matrix("stiffness");
    libMesh::SparseMatrix<libMesh::Number>&  damping   = system.get_matrix("damping");
    libMesh::SparseMatrix<libMesh::Number>&  mass      = system.get_matrix("mass");
    libMesh::NumericVector<libMesh::Number>& force     = system.get_vector("force");
    libMesh::SparseMatrix<libMesh::Number>&  matrix    = *system.matrix;
    libMesh::DenseMatrix<libMesh::Number>    zero_matrix;
    //
	// Set up pointers to FEBase's of dimension dim and FE type fe_type
	// -> 3D elements
	libMesh::UniquePtr<libMesh::FEBase> fe (libMesh::FEBase::build(dim, fe_type));
	libMesh::QGauss qrule (dim, fe_type.default_quadrature_order());
	fe->attach_quadrature_rule (&qrule);

	// -> Faces
	libMesh::UniquePtr<libMesh::FEBase> fe_face (libMesh::FEBase::build(dim, fe_type));
	libMesh::QGauss qface(dim-1, fe_type.default_quadrature_order());
	fe_face->attach_quadrature_rule (&qface);

	// Jacobian
	const std::vector<libMesh::Real>& JxW = fe->get_JxW();

	// Shape functions
	const std::vector<std::vector<libMesh::Real> >& phi = fe->get_phi();

	// Shape functions derivatives
	const std::vector<std::vector<libMesh::RealGradient> >& dphi = fe->get_dphi();

	// Quadrature points
	const std::vector<libMesh::Point>& qp_points = fe->get_xyz();

	// Weights for the Arlequin method
	double weight_alpha = 1;

	libMesh::DenseMatrix<libMesh::Number> Me;
	libMesh::DenseMatrix<libMesh::Number> Ce;
	libMesh::DenseMatrix<libMesh::Number> Ke;
	libMesh::DenseVector<libMesh::Number> Fe;

	libMesh::DenseSubMatrix<libMesh::Number>
	Muu(Me), Mvv(Me), Mww(Me);

	libMesh::DenseSubMatrix<libMesh::Number>
	Cuu(Ce), Cuv(Ce), Cuw(Ce),
	Cvu(Ce), Cvv(Ce), Cvw(Ce),
	Cwu(Ce), Cwv(Ce), Cww(Ce);

	libMesh::DenseSubMatrix<libMesh::Number>
	Kuu(Ke), Kuv(Ke), Kuw(Ke),
	Kvu(Ke), Kvv(Ke), Kvw(Ke),
	Kwu(Ke), Kwv(Ke), Kww(Ke);

	libMesh::DenseSubVector<libMesh::Number>
	Fu(Fe),
	Fv(Fe),
	Fw(Fe);

	std::vector<libMesh::dof_id_type> dof_indices;
	std::vector<libMesh::dof_id_type> dof_indices_u;
	std::vector<libMesh::dof_id_type> dof_indices_v;
	std::vector<libMesh::dof_id_type> dof_indices_w;

	libMesh::MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
	const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	perf_log.pop("Preamble");

	///// For each element
	///for ( ; el != end_el; ++el)
	///{
	///	// Get its pointer
	///	const libMesh::Elem* elem = *el;

	///	perf_log.push("Define DoF");
	///	// The total DoF indices, and those associated to each variable
	///	dof_map.dof_indices (elem, dof_indices);
	///	dof_map.dof_indices (elem, dof_indices_u, u_var);
	///	dof_map.dof_indices (elem, dof_indices_v, v_var);
	///	dof_map.dof_indices (elem, dof_indices_w, w_var);

	///	const unsigned int n_dofs   = dof_indices.size();
	///	const unsigned int n_u_dofs = dof_indices_u.size();
	///	const unsigned int n_v_dofs = dof_indices_v.size();
	///	const unsigned int n_w_dofs = dof_indices_w.size();

	///	perf_log.pop("Define DoF");

	///	// Restart the FE to the "geometry" of the element
	///	// -> Determines quadrature points, shape functions ...
	///	// !!! User can change the points to be used (can use other mesh's points
	///	//		instead of the quadrature points)
	///	fe->reinit (elem);

	///	perf_log.push("Matrix manipulations");
	///	Me.resize (n_dofs, n_dofs);
	///	Ce.resize (n_dofs, n_dofs);
	///	Ke.resize (n_dofs, n_dofs);
    ///    zero_matrix.resize(n_dofs,n_dofs);
	///	Fe.resize (n_dofs);

	///	// Set the positions of the sub-matrices
	///	Muu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
	///	Mvv.reposition (v_var*n_u_dofs, v_var*n_u_dofs, n_v_dofs, n_v_dofs);
	///	Mww.reposition (w_var*n_u_dofs, w_var*n_u_dofs, n_w_dofs, n_w_dofs);

	///	Cuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
	///	Cuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
	///	Cuw.reposition (u_var*n_u_dofs, w_var*n_u_dofs, n_u_dofs, n_w_dofs);

	///	Cvu.reposition (v_var*n_u_dofs, u_var*n_u_dofs, n_v_dofs, n_u_dofs);
	///	Cvv.reposition (v_var*n_u_dofs, v_var*n_u_dofs, n_v_dofs, n_v_dofs);
	///	Cvw.reposition (v_var*n_u_dofs, w_var*n_u_dofs, n_v_dofs, n_w_dofs);

	///	Cwu.reposition (w_var*n_u_dofs, u_var*n_u_dofs, n_w_dofs, n_u_dofs);
	///	Cwv.reposition (w_var*n_u_dofs, v_var*n_u_dofs, n_w_dofs, n_v_dofs);
	///	Cww.reposition (w_var*n_u_dofs, w_var*n_u_dofs, n_w_dofs, n_w_dofs);

	///	Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
	///	Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
	///	Kuw.reposition (u_var*n_u_dofs, w_var*n_u_dofs, n_u_dofs, n_w_dofs);

	///	Kvu.reposition (v_var*n_u_dofs, u_var*n_u_dofs, n_v_dofs, n_u_dofs);
	///	Kvv.reposition (v_var*n_u_dofs, v_var*n_u_dofs, n_v_dofs, n_v_dofs);
	///	Kvw.reposition (v_var*n_u_dofs, w_var*n_u_dofs, n_v_dofs, n_w_dofs);

	///	Kwu.reposition (w_var*n_u_dofs, u_var*n_u_dofs, n_w_dofs, n_u_dofs);
	///	Kwv.reposition (w_var*n_u_dofs, v_var*n_u_dofs, n_w_dofs, n_v_dofs);
	///	Kww.reposition (w_var*n_u_dofs, w_var*n_u_dofs, n_w_dofs, n_w_dofs);

	///	Fu.reposition (u_var*n_u_dofs, n_u_dofs);
	///	Fv.reposition (v_var*n_u_dofs, n_v_dofs);
	///	Fw.reposition (w_var*n_u_dofs, n_w_dofs);
	///	perf_log.push("Matrix manipulations");

    ///    // For each quadrature point determinate the sub-matrices elements
	///	for (unsigned int qp=0; qp<qrule.n_points(); qp++)
	///	{

	///		perf_log.push("Mass","Matrix element calculations");
	///		weight_alpha = weight_mask.get_alpha(qp_points[qp],system_type);
    ///        
	///		// Mass Matrix
	///		Update_SubM_isotropic(Muu, qp, phi, n_u_dofs, JxW, localRho, weight_alpha);
	///		Update_SubM_isotropic(Mvv, qp, phi, n_u_dofs, JxW, localRho, weight_alpha);
	///		Update_SubM_isotropic(Mww, qp, phi, n_u_dofs, JxW, localRho, weight_alpha);

	///		perf_log.pop("Mass","Matrix element calculations");
    ///        
	///		// Internal tension
	///		perf_log.push("Rigidity","Matrix element calculations");
	///		Update_SubK_isotropic(Kuu, qp, 0, 0, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
	///		Update_SubK_isotropic(Kuv, qp, 0, 1, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
	///		Update_SubK_isotropic(Kuw, qp, 0, 2, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);

	///		Update_SubK_isotropic(Kvu, qp, 1, 0, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
	///		Update_SubK_isotropic(Kvv, qp, 1, 1, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
	///		Update_SubK_isotropic(Kvw, qp, 1, 2, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);

	///		Update_SubK_isotropic(Kwu, qp, 2, 0, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
	///		Update_SubK_isotropic(Kwv, qp, 2, 1, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
	///		Update_SubK_isotropic(Kww, qp, 2, 2, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);

    ///        // damping ??
	///		Update_SubC_isotropic(Cuu, Muu, Kuu, qp, n_u_dofs, weight_alpha,0.0,0.0);
	///		Update_SubC_isotropic(Cuv, Muu, Kuv, qp, n_u_dofs, weight_alpha,0.0,0.0);
	///		Update_SubC_isotropic(Cuw, Muu, Kuw, qp, n_u_dofs, weight_alpha,0.0,0.0);

	///		Update_SubC_isotropic(Cvu, Mvv, Kvu, qp, n_u_dofs, weight_alpha,0.0,0.0);
	///		Update_SubC_isotropic(Cvv, Mvv, Kvv, qp, n_u_dofs, weight_alpha,0.0,0.0);
	///		Update_SubC_isotropic(Cvw, Mvv, Kvw, qp, n_u_dofs, weight_alpha,0.0,0.0);

	///		Update_SubC_isotropic(Cwu, Mww, Kwu, qp, n_u_dofs, weight_alpha,0.0,0.0);
	///		Update_SubC_isotropic(Cwv, Mww, Kwv, qp, n_u_dofs, weight_alpha,0.0,0.0);
	///		Update_SubC_isotropic(Cww, Mww, Kww, qp, n_u_dofs, weight_alpha,0.0,0.0);
    ///        
	///		perf_log.pop("Rigidity","Matrix element calculations");
	///	}

	///	// Apply constraints
	///	perf_log.push("Constraints","Matrix element calculations");
	///	dof_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
	///	perf_log.pop("Constraints","Matrix element calculations");

	///	perf_log.push("Adding elements");
	///	system.mass->add_matrix (Me, dof_indices);
	///	system.stiffness->add_matrix (Ke, dof_indices);
	///	system.damping->add_matrix (Ce, dof_indices);

	///	system.rhs->add_vector    (Fe, dof_indices);
	///	perf_log.pop("Adding elements");
	///}

	///system.stiffness->close();
    ///system.mass->close();
    ///system.damping->close();
	///system.rhs->close();
}
////
//void assemble_elasticity_with_weight_and_traction_dyn(libMesh::EquationSystems& es,
//							const std::string& system_name, 
//							weight_parameter_function& weight_mask,
//							WeightFunctionSystemType system_type,
//							int traction_boundary_id,
//							std::vector<double> traction_density)
//
////{
//	libmesh_assert_equal_to(system_name, "Elasticity");
//
//	libMesh::PerfLog perf_log ("Matrix Assembly (Homogeneous elasticity)",MASTER_bPerfLog_assemble_fem);
//
//	perf_log.push("Preamble");
//
//	const libMesh::MeshBase& mesh = es.get_mesh();
//
//	const unsigned int dim = mesh.mesh_dimension();
//
//	// - Set up physical properties system ------------------------------------
//	libMesh::ExplicitSystem& physical_param_system = es.get_system<libMesh::ExplicitSystem>("PhysicalConstants");
//	const unsigned int young_var = physical_param_system.variable_number ("E");
//	const unsigned int mu_var = physical_param_system.variable_number ("mu");
//	const unsigned int rho_var = physical_param_system.variable_number ("rho");
//
//	const libMesh::DofMap& physical_dof_map = physical_param_system.get_dof_map();
//	std::vector<libMesh::dof_id_type> physical_dof_indices_var;
//
//	// The DoF and values of the physical system
//	const libMesh::Elem*  phys_elem   = *(mesh.active_local_elements_begin());
//	physical_dof_map.dof_indices(phys_elem, physical_dof_indices_var, young_var);
//	libMesh::Number localE = physical_param_system.current_solution(physical_dof_indices_var[0]);
//
//	physical_dof_map.dof_indices(phys_elem, physical_dof_indices_var, mu_var);
//	libMesh::Number localMu = physical_param_system.current_solution(physical_dof_indices_var[0]);
//
//	physical_dof_map.dof_indices(phys_elem, physical_dof_indices_var, rho_var);
//	libMesh::Number localRho = physical_param_system.current_solution(physical_dof_indices_var[0]);
//	// - Set up elasticity system ---------------------------------------------
//	libMesh::LinearImplicitSystem& system = es.get_system<libMesh::LinearImplicitSystem>("Elasticity");
//
//	const unsigned int n_components = 3;
//	const unsigned int u_var = system.variable_number ("u");
//	const unsigned int v_var = system.variable_number ("v");
//	const unsigned int w_var = system.variable_number ("w");
//
//	const libMesh::DofMap& dof_map = system.get_dof_map();
//	libMesh::FEType fe_type = dof_map.variable_type(u_var);
//
//	// Set up pointers to FEBase's of dimension dim and FE type fe_type
//	// -> 3D elements
//	libMesh::UniquePtr<libMesh::FEBase> fe (libMesh::FEBase::build(dim, fe_type));
//	libMesh::QGauss qrule (dim, fe_type.default_quadrature_order());
//	fe->attach_quadrature_rule (&qrule);
//
//	// -> Faces
//	libMesh::UniquePtr<libMesh::FEBase> fe_face (libMesh::FEBase::build(dim, fe_type));
//	libMesh::QGauss qface(dim-1, fe_type.default_quadrature_order());
//	fe_face->attach_quadrature_rule (&qface);
//
//	// -> Face traction
//	libMesh::DenseVector<libMesh::Number> g_vec(3);
//	g_vec(0) = traction_density[0];
//	g_vec(1) = traction_density[1];
//	g_vec(2) = traction_density[2];
//
//	// Jacobian
//	const std::vector<libMesh::Real>& JxW = fe->get_JxW();
//
//	// Shape functions
//	const std::vector<std::vector<libMesh::Real> >& phi = fe->get_phi();
//   
//	// Shape functions derivatives
//	const std::vector<std::vector<libMesh::RealGradient> >& dphi = fe->get_dphi();
//
//	// Quadrature points
//	const std::vector<libMesh::Point>& qp_points = fe->get_xyz();
//
//	// Face Jacobian
//	const std::vector<libMesh::Real> & JxW_face = fe_face->get_JxW();
//
//	// Face shape functions
//	const std::vector<std::vector<libMesh::Real> > & phi_face = fe_face->get_phi();
//
//	// Face quadrature points
//	const std::vector<libMesh::Point>& qp_face_points = fe_face->get_xyz();
//
//
//	// Weights for the Arlequin method
//	double weight_alpha = 1;
//
//	libMesh::DenseMatrix<libMesh::Number> Me;
//	libMesh::DenseMatrix<libMesh::Number> Ke;
//	libMesh::DenseVector<libMesh::Number> Fe;
//
//	libMesh::DenseSubMatrix<libMesh::Number>
//	Muu(Me), Mvv(Me), Mww(Me);
//
//	libMesh::DenseSubMatrix<libMesh::Number>
//	Kuu(Ke), Kuv(Ke), Kuw(Ke),
//	Kvu(Ke), Kvv(Ke), Kvw(Ke),
//	Kwu(Ke), Kwv(Ke), Kww(Ke);
//
//	libMesh::DenseSubVector<libMesh::Number>
//	Fu(Fe),
//	Fv(Fe),
//	Fw(Fe);
//
//	std::vector<libMesh::dof_id_type> dof_indices;
//	std::vector<libMesh::dof_id_type> dof_indices_u;
//	std::vector<libMesh::dof_id_type> dof_indices_v;
//	std::vector<libMesh::dof_id_type> dof_indices_w;
//
//	libMesh::MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
//	const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
//
//	perf_log.pop("Preamble");
//
//	// For each element
//	for ( ; el != end_el; ++el)
//	{
//		// Get its pointer
//		const libMesh::Elem* elem = *el;
//
//		perf_log.push("Define DoF");
//		// The total DoF indices, and those associated to each variable
//		dof_map.dof_indices (elem, dof_indices);
//		dof_map.dof_indices (elem, dof_indices_u, u_var);
//		dof_map.dof_indices (elem, dof_indices_v, v_var);
//		dof_map.dof_indices (elem, dof_indices_w, w_var);
//
//		const unsigned int n_dofs   = dof_indices.size();
//		const unsigned int n_u_dofs = dof_indices_u.size();
//		const unsigned int n_v_dofs = dof_indices_v.size();
//		const unsigned int n_w_dofs = dof_indices_w.size();
//
//		perf_log.pop("Define DoF");
//
//		// Restart the FE to the "geometry" of the element
//		// -> Determines quadrature points, shape functions ...
//		// !!! User can change the points to be used (can use other mesh's points
//		//		instead of the quadrature points)
//		fe->reinit (elem);
//
//		perf_log.push("Matrix manipulations");
//		Me.resize (n_dofs, n_dofs);
//		Ke.resize (n_dofs, n_dofs);
//		Fe.resize (n_dofs);
//
//		// Set the positions of the sub-matrices
//		Muu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
//		Mvv.reposition (v_var*n_u_dofs, v_var*n_u_dofs, n_v_dofs, n_v_dofs);
//		Mww.reposition (w_var*n_u_dofs, w_var*n_u_dofs, n_w_dofs, n_w_dofs);
//        
//		// Set the positions of the sub-matrices
//		Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
//		Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
//		Kuw.reposition (u_var*n_u_dofs, w_var*n_u_dofs, n_u_dofs, n_w_dofs);
//
//		Kvu.reposition (v_var*n_u_dofs, u_var*n_u_dofs, n_v_dofs, n_u_dofs);
//		Kvv.reposition (v_var*n_u_dofs, v_var*n_u_dofs, n_v_dofs, n_v_dofs);
//		Kvw.reposition (v_var*n_u_dofs, w_var*n_u_dofs, n_v_dofs, n_w_dofs);
//
//		Kwu.reposition (w_var*n_u_dofs, u_var*n_u_dofs, n_w_dofs, n_u_dofs);
//		Kwv.reposition (w_var*n_u_dofs, v_var*n_u_dofs, n_w_dofs, n_v_dofs);
//		Kww.reposition (w_var*n_u_dofs, w_var*n_u_dofs, n_w_dofs, n_w_dofs);
//
//		Fu.reposition (u_var*n_u_dofs, n_u_dofs);
//		Fv.reposition (v_var*n_u_dofs, n_v_dofs);
//		Fw.reposition (w_var*n_u_dofs, n_w_dofs);
//		perf_log.push("Matrix manipulations");
//
//		// For each quadrature point determinate the sub-matrices elements
//		for (unsigned int qp=0; qp<qrule.n_points(); qp++)
//		{
//
//			perf_log.push("Rigidity","Matrix element calculations");
//			weight_alpha = weight_mask.get_alpha(qp_points[qp],system_type);
//
//			// Mass Matrix
//			Update_SubM_isotropic(Muu, qp, phi, n_u_dofs, JxW, localRho, weight_alpha);
//			Update_SubM_isotropic(Mvv, qp, phi, n_u_dofs, JxW, localRho, weight_alpha);
//			Update_SubM_isotropic(Mww, qp, phi, n_u_dofs, JxW, localRho, weight_alpha);
//
//			perf_log.pop("Mass","Matrix element calculations");
//            
//			// Internal tension
//			perf_log.push("Rigidity","Matrix element calculations");
//			Update_SubK_isotropic(Kuu, qp, 0, 0, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
//			Update_SubK_isotropic(Kuv, qp, 0, 1, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
//			Update_SubK_isotropic(Kuw, qp, 0, 2, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
//
//			Update_SubK_isotropic(Kvu, qp, 1, 0, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
//			Update_SubK_isotropic(Kvv, qp, 1, 1, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
//			Update_SubK_isotropic(Kvw, qp, 1, 2, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
//
//			Update_SubK_isotropic(Kwu, qp, 2, 0, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
//			Update_SubK_isotropic(Kwv, qp, 2, 1, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
//			Update_SubK_isotropic(Kww, qp, 2, 2, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
//
//			perf_log.pop("Rigidity","Matrix element calculations");
//		}
//
//		// Assemble the traction
//		for (unsigned int side=0; side<elem->n_sides(); side++)
//		{
//			if (elem->neighbor_ptr(side) == libmesh_nullptr)
//			{
//				fe_face->reinit(elem, side);
//
//				// Apply a traction
//				for (unsigned int qp=0; qp<qface.n_points(); qp++)
//				{
//					weight_alpha = weight_mask.get_alpha(qp_face_points[qp],system_type);
//
//					if (mesh.get_boundary_info().has_boundary_id(elem, side, traction_boundary_id))
//					{	
//						for (unsigned int dof_iii=0; dof_iii<n_u_dofs; dof_iii++)
//						{
//							Fu(dof_iii) += weight_alpha * JxW_face[qp] * (g_vec(0) * phi_face[dof_iii][qp]);
//							Fv(dof_iii) += weight_alpha * JxW_face[qp] * (g_vec(1) * phi_face[dof_iii][qp]);
//							Fw(dof_iii) += weight_alpha * JxW_face[qp] * (g_vec(2) * phi_face[dof_iii][qp]);
//						}
//					}
//				}
//			}
//		}
//
//		// Apply constraints
//		perf_log.push("Constraints","Matrix element calculations");
//		dof_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
//		dof_map.constrain_element_matrix(Me, dof_indices);
//		perf_log.pop("Constraints","Matrix element calculations");
//
//		perf_log.push("Adding elements");
//		system.matrix->add_matrix (Me, dof_indices);
//		system.matrix->add_matrix (Ke, dof_indices);
//		system.rhs->add_vector    (Fe, dof_indices);
//		perf_log.pop("Adding elements");
//	}
//
//	system.matrix->close();
//	system.rhs->close();
//};
////
//void assemble_elasticity_heterogeneous_with_weight_dyn(	libMesh::EquationSystems& es,
//							const std::string& system_name,
//							weight_parameter_function& weight_mask,
//							WeightFunctionSystemType system_type)
//{
//	libmesh_assert_equal_to (system_name, "Elasticity");
//
//	libMesh::PerfLog perf_log ("Matrix Assembly (Heterogeneous elasticity)",MASTER_bPerfLog_assemble_fem);
//
//	perf_log.push("Preamble");
//	// Set up mesh
//	const libMesh::MeshBase& mesh = es.get_mesh();
//
//	const unsigned int dim = mesh.mesh_dimension();
//
//	// - Set up physical properties system ------------------------------------
//	libMesh::Number localE = -1;
//	libMesh::Number localMu = -1;
//	libMesh::Number localRho = -1;
//
//	libMesh::ExplicitSystem& physical_param_system = es.get_system<libMesh::ExplicitSystem>("PhysicalConstants");
//	const unsigned int young_var = physical_param_system.variable_number ("E");
//	const unsigned int mu_var = physical_param_system.variable_number ("mu");
//	const unsigned int rho_var = physical_param_system.variable_number ("rho");
//
//	const libMesh::DofMap& physical_dof_map = physical_param_system.get_dof_map();
//	std::vector<libMesh::dof_id_type> physical_dof_indices_var;
//
//	// - Set up elasticity system ---------------------------------------------
//	libMesh::NewmarkSystem& system = es.get_system<libMesh::NewmarkSystem>("Transient");
//
//	const unsigned int n_components = 3;
//	const unsigned int u_var = system.variable_number ("u");
//	const unsigned int v_var = system.variable_number ("v");
//	const unsigned int w_var = system.variable_number ("w");
//
//	const libMesh::DofMap& dof_map = system.get_dof_map();
//	libMesh::FEType fe_type = dof_map.variable_type(u_var);
//
//	// Set up pointers to FEBase's of dimension dim and FE type fe_type
//	// -> 3D elements
//	libMesh::UniquePtr<libMesh::FEBase> fe (libMesh::FEBase::build(dim, fe_type));
//	libMesh::QGauss qrule (dim, fe_type.default_quadrature_order());
//	fe->attach_quadrature_rule (&qrule);
//
//	// -> Faces
//	libMesh::UniquePtr<libMesh::FEBase> fe_face (libMesh::FEBase::build(dim, fe_type));
//	libMesh::QGauss qface(dim-1, fe_type.default_quadrature_order());
//	fe_face->attach_quadrature_rule (&qface);
//
//	// Jacobian
//	const std::vector<libMesh::Real>& JxW = fe->get_JxW();
//
//	// Shape functions
//	const std::vector<std::vector<libMesh::Real> >& phi = fe->get_phi();
//    
//	// Shape functions derivatives
//	const std::vector<std::vector<libMesh::RealGradient> >& dphi = fe->get_dphi();
//
//	// Quadrature points
//	const std::vector<libMesh::Point>& qp_points = fe->get_xyz();
//
//	// Weights for the Arlequin method
//	double weight_alpha = 1;
//
//	libMesh::DenseMatrix<libMesh::Number> Me;
//	libMesh::DenseMatrix<libMesh::Number> Ke;
//	libMesh::DenseVector<libMesh::Number> Fe;
//
//	libMesh::DenseSubMatrix<libMesh::Number>
//	Muu(Me), Mvv(Me), Mww(Me);
//
//	libMesh::DenseSubMatrix<libMesh::Number>
//	Kuu(Ke), Kuv(Ke), Kuw(Ke),
//	Kvu(Ke), Kvv(Ke), Kvw(Ke),
//	Kwu(Ke), Kwv(Ke), Kww(Ke);
//
//	libMesh::DenseSubVector<libMesh::Number>
//	Fu(Fe),
//	Fv(Fe),
//	Fw(Fe);
//
//	std::vector<libMesh::dof_id_type> dof_indices;
//	std::vector<libMesh::dof_id_type> dof_indices_u;
//	std::vector<libMesh::dof_id_type> dof_indices_v;
//	std::vector<libMesh::dof_id_type> dof_indices_w;
//
//	libMesh::MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
//	const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
//
//	perf_log.pop("Preamble");
//
//	// For each element
//	for ( ; el != end_el; ++el)
//	{
//		// Get its pointer
//		const libMesh::Elem* elem = *el;
//
//		perf_log.push("Define DoF");
//		// The total DoF indices, and those associated to each variable
//		dof_map.dof_indices (elem, dof_indices);
//		dof_map.dof_indices (elem, dof_indices_u, u_var);
//		dof_map.dof_indices (elem, dof_indices_v, v_var);
//		dof_map.dof_indices (elem, dof_indices_w, w_var);
//
//		const unsigned int n_dofs   = dof_indices.size();
//		const unsigned int n_u_dofs = dof_indices_u.size();
//		const unsigned int n_v_dofs = dof_indices_v.size();
//		const unsigned int n_w_dofs = dof_indices_w.size();
//
//		perf_log.pop("Define DoF");
//
//		perf_log.push("Define physical params");
//		// The DoF and values of the physical system
//		physical_dof_map.dof_indices(elem, physical_dof_indices_var, young_var);
//		localE = physical_param_system.current_solution(physical_dof_indices_var[0]);
//
//		physical_dof_map.dof_indices(elem, physical_dof_indices_var, mu_var);
//		localMu = physical_param_system.current_solution(physical_dof_indices_var[0]);
//
//		physical_dof_map.dof_indices(elem, physical_dof_indices_var, rho_var);
//		localRho = physical_param_system.current_solution(physical_dof_indices_var[0]);
//		perf_log.pop("Define physical params");
//
//		// Restart the FE to the "geometry" of the element
//		// -> Determines quadrature points, shape functions ...
//		// !!! User can change the points to be used (can use other mesh's points
//		//		instead of the quadrature points)
//		fe->reinit (elem);
//
//		perf_log.push("Matrix manipulations");
//		Me.resize (n_dofs, n_dofs);
//		Ke.resize (n_dofs, n_dofs);
//		Fe.resize (n_dofs);
//
//		Muu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
//		Mvv.reposition (v_var*n_u_dofs, v_var*n_u_dofs, n_v_dofs, n_v_dofs);
//		Mww.reposition (w_var*n_u_dofs, w_var*n_u_dofs, n_w_dofs, n_w_dofs);
//       
//		// Set the positions of the sub-matrices
//		Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
//		Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
//		Kuw.reposition (u_var*n_u_dofs, w_var*n_u_dofs, n_u_dofs, n_w_dofs);
//
//		Kvu.reposition (v_var*n_u_dofs, u_var*n_u_dofs, n_v_dofs, n_u_dofs);
//		Kvv.reposition (v_var*n_u_dofs, v_var*n_u_dofs, n_v_dofs, n_v_dofs);
//		Kvw.reposition (v_var*n_u_dofs, w_var*n_u_dofs, n_v_dofs, n_w_dofs);
//
//		Kwu.reposition (w_var*n_u_dofs, u_var*n_u_dofs, n_w_dofs, n_u_dofs);
//		Kwv.reposition (w_var*n_u_dofs, v_var*n_u_dofs, n_w_dofs, n_v_dofs);
//		Kww.reposition (w_var*n_u_dofs, w_var*n_u_dofs, n_w_dofs, n_w_dofs);
//
//		Fu.reposition (u_var*n_u_dofs, n_u_dofs);
//		Fv.reposition (v_var*n_u_dofs, n_v_dofs);
//		Fw.reposition (w_var*n_u_dofs, n_w_dofs);
//		perf_log.push("Matrix manipulations");
//
//		perf_log.push("Matrix element calculations");
//		// For each quadrature point determinate the sub-matrices elements
//		for (unsigned int qp=0; qp<qrule.n_points(); qp++)
//		{
//
//			perf_log.push("Rigidity","Matrix element calculations");
//			weight_alpha = weight_mask.get_alpha(qp_points[qp],system_type);
//
//            // Mass matrix
//			Update_SubM_isotropic(Muu, qp, phi, n_u_dofs, JxW, localRho, weight_alpha);
//			Update_SubM_isotropic(Mvv, qp, phi, n_u_dofs, JxW, localRho, weight_alpha);
//			Update_SubM_isotropic(Mww, qp, phi, n_u_dofs, JxW, localRho, weight_alpha);
//
//			// Internal tension
//			Update_SubK_isotropic(Kuu, qp, 0, 0, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
//			Update_SubK_isotropic(Kuv, qp, 0, 1, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
//			Update_SubK_isotropic(Kuw, qp, 0, 2, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
//
//			Update_SubK_isotropic(Kvu, qp, 1, 0, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
//			Update_SubK_isotropic(Kvv, qp, 1, 1, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
//			Update_SubK_isotropic(Kvw, qp, 1, 2, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
//
//			Update_SubK_isotropic(Kwu, qp, 2, 0, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
//			Update_SubK_isotropic(Kwv, qp, 2, 1, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
//			Update_SubK_isotropic(Kww, qp, 2, 2, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
//
//			perf_log.pop("Rigidity","Matrix element calculations");
//		}
//
//		// Apply constraints
//		perf_log.push("Constraints","Matrix element calculations");
//		dof_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
//		perf_log.pop("Constraints","Matrix element calculations");
//
//		perf_log.push("Adding elements");
//		system.matrix->add_matrix (Me, dof_indices);
//		system.matrix->add_matrix (Ke, dof_indices);
//		system.rhs->add_vector    (Fe, dof_indices);
//		perf_log.pop("Adding elements");
//	}
//	
//	system.matrix->close();
//	system.rhs->close();
//};
