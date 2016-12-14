/*
 * assemble_intersection_3D.h
 *
 *  Created on: Nov 30, 2015
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef ASSEMBLE_INTERSECTION_3D_H_
#define ASSEMBLE_INTERSECTION_3D_H_

#include "carl_headers.h"

#include "mpi_carl_tools.h"
#include "weak_formulations.h"
#include "coupled_solver.h"

#include "LATIN_solver.h"
#include "CG_solver.h"

#include "PETSC_matrix_operations.h"
#include "weight_parameter_function.h"
#include "anisotropic_elasticity_cubic_sym.h"

const bool MASTER_bPerfLog_assemble_coupling = false;
const bool MASTER_debug_coupling_assemble = false;

namespace carl
{

class libMesh_fe_addresses_3
{
private:
	// Private default constructor
	libMesh_fe_addresses_3();

public:

	// Constructor
	libMesh_fe_addresses_3(libMesh::System& input_system,
			const std::string u_var_name = "u", const std::string v_var_name =
					"v", const std::string w_var_name = "w") :
			eq_system { input_system },
			mesh { eq_system.get_mesh() },
			dim { mesh.mesh_dimension() },
			u_var { input_system.variable_number(u_var_name) },
			v_var { input_system.variable_number(v_var_name) },
			w_var { input_system.variable_number(w_var_name) },
			dof_map { input_system.get_dof_map() },
			fe_type { dof_map.variable_type(u_var) },
			fe_unique_ptr { libMesh::FEBase::build(dim, fe_type) },
			qrule { dim, fe_type.default_quadrature_order() },
			n_dofs { 0 },
			n_dofs_u { 0 },
			n_dofs_v { 0 },
			n_dofs_w { 0 }

	{
		fe_unique_ptr->attach_quadrature_rule(&qrule);
	};

	// Members set at initialization
	libMesh::System& eq_system;
	const libMesh::MeshBase& mesh;
	const unsigned int dim;
	const unsigned int u_var;
	const unsigned int v_var;
	const unsigned int w_var;
	const libMesh::DofMap& dof_map;
	libMesh::FEType fe_type;
	libMesh::UniquePtr<libMesh::FEBase> fe_unique_ptr;
	libMesh::QGauss qrule;

	// Members set with the set_dofs method
	std::vector<libMesh::dof_id_type> dof_indices;
	std::vector<libMesh::dof_id_type> dof_indices_u;
	std::vector<libMesh::dof_id_type> dof_indices_v;
	std::vector<libMesh::dof_id_type> dof_indices_w;

	unsigned int n_dofs;
	unsigned int n_dofs_u;
	unsigned int n_dofs_v;
	unsigned int n_dofs_w;

	void set_DoFs(int idx = 0)
	{
		const libMesh::Elem* elem = mesh.elem(idx);
		dof_map.dof_indices(elem, dof_indices);
		dof_map.dof_indices(elem, dof_indices_u, u_var);
		dof_map.dof_indices(elem, dof_indices_v, v_var);
		dof_map.dof_indices(elem, dof_indices_w, w_var);
		n_dofs = dof_indices.size();
		n_dofs_u = dof_indices_u.size();
		n_dofs_v = dof_indices_v.size();
		n_dofs_w = dof_indices_w.size();
	};
};

class coupling_matrices_3
{
public:
	libMesh::DenseMatrix<libMesh::Number> Me;
	libMesh::DenseSubMatrix<libMesh::Number> Me_uu, Me_uv, Me_uw, Me_vu, Me_vv,
			Me_vw, Me_wu, Me_wv, Me_ww;

	coupling_matrices_3() :
			Me_uu
			{ Me }, Me_uv
			{ Me }, Me_uw
			{ Me }, Me_vu
			{ Me }, Me_vv
			{ Me }, Me_vw
			{ Me }, Me_wu
			{ Me }, Me_wv
			{ Me }, Me_ww
			{ Me }
	{

	}

	void set_matrices(libMesh_fe_addresses_3& system_type_AAA,
			libMesh_fe_addresses_3& system_type_BBB)
	{
		Me.resize(system_type_AAA.n_dofs, system_type_BBB.n_dofs);

		Me_uu.reposition(system_type_AAA.u_var * system_type_AAA.n_dofs_u,
				system_type_BBB.u_var * system_type_BBB.n_dofs_u,
				system_type_AAA.n_dofs_u, system_type_BBB.n_dofs_u);

		Me_uv.reposition(system_type_AAA.u_var * system_type_AAA.n_dofs_u,
				system_type_BBB.v_var * system_type_BBB.n_dofs_v,
				system_type_AAA.n_dofs_u, system_type_BBB.n_dofs_v);

		Me_uw.reposition(system_type_AAA.u_var * system_type_AAA.n_dofs_u,
				system_type_BBB.w_var * system_type_BBB.n_dofs_w,
				system_type_AAA.n_dofs_u, system_type_BBB.n_dofs_w);

		Me_vu.reposition(system_type_AAA.v_var * system_type_AAA.n_dofs_v,
				system_type_BBB.u_var * system_type_BBB.n_dofs_u,
				system_type_AAA.n_dofs_v, system_type_BBB.n_dofs_u);

		Me_vv.reposition(system_type_AAA.v_var * system_type_AAA.n_dofs_v,
				system_type_BBB.v_var * system_type_BBB.n_dofs_v,
				system_type_AAA.n_dofs_v, system_type_BBB.n_dofs_v);

		Me_vw.reposition(system_type_AAA.v_var * system_type_AAA.n_dofs_v,
				system_type_BBB.w_var * system_type_BBB.n_dofs_w,
				system_type_AAA.n_dofs_v, system_type_BBB.n_dofs_w);

		Me_wu.reposition(system_type_AAA.w_var * system_type_AAA.n_dofs_w,
				system_type_BBB.u_var * system_type_BBB.n_dofs_u,
				system_type_AAA.n_dofs_w, system_type_BBB.n_dofs_u);

		Me_wv.reposition(system_type_AAA.w_var * system_type_AAA.n_dofs_w,
				system_type_BBB.v_var * system_type_BBB.n_dofs_v,
				system_type_AAA.n_dofs_w, system_type_BBB.n_dofs_v);

		Me_ww.reposition(system_type_AAA.w_var * system_type_AAA.n_dofs_w,
				system_type_BBB.w_var * system_type_BBB.n_dofs_w,
				system_type_AAA.n_dofs_w, system_type_BBB.n_dofs_w);
	}

	void build_L2_coupling_matrix(const libMesh_fe_addresses_3& system_type_AAA,
			const libMesh_fe_addresses_3& system_type_BBB, int qp,
			const std::vector<std::vector<libMesh::Real> >& corrected_phi_AAA,
			const std::vector<std::vector<libMesh::Real> >& corrected_phi_BBB,
			const std::vector<libMesh::Real>& JxW, double L2_coupling_const)
	{
		L2_Coupling(Me_uu, qp, corrected_phi_AAA, corrected_phi_BBB,
				system_type_AAA.n_dofs_u, system_type_BBB.n_dofs_u, JxW,
				L2_coupling_const);

		L2_Coupling(Me_vv, qp, corrected_phi_AAA, corrected_phi_BBB,
				system_type_AAA.n_dofs_v, system_type_BBB.n_dofs_v, JxW,
				L2_coupling_const);

		L2_Coupling(Me_ww, qp, corrected_phi_AAA, corrected_phi_BBB,
				system_type_AAA.n_dofs_w, system_type_BBB.n_dofs_w, JxW,
				L2_coupling_const);
	}

	void add_H1_coupling_matrix(const libMesh_fe_addresses_3& system_type_AAA,
			const libMesh_fe_addresses_3& system_type_BBB, int qp,
			const std::vector<std::vector<libMesh::RealGradient> >& corrected_dphi_sysAAA,
			const std::vector<std::vector<libMesh::RealGradient> >& corrected_dphi_sysBBB,
			const std::vector<libMesh::Real>& JxW,
			const libMesh::Number H1_coupling_const)
	{
		unsigned int n_components = 3;
		H1_Coupling_Extra_Term(Me_uu, qp, system_type_AAA.u_var,
				system_type_BBB.u_var, n_components, n_components,
				corrected_dphi_sysAAA, corrected_dphi_sysBBB,
				system_type_AAA.n_dofs_u, system_type_BBB.n_dofs_u, JxW,
				H1_coupling_const);

		H1_Coupling_Extra_Term(Me_uv, qp, system_type_AAA.u_var,
				system_type_BBB.v_var, n_components, n_components,
				corrected_dphi_sysAAA, corrected_dphi_sysBBB,
				system_type_AAA.n_dofs_u, system_type_BBB.n_dofs_v, JxW,
				H1_coupling_const);

		H1_Coupling_Extra_Term(Me_uw, qp, system_type_AAA.u_var,
				system_type_BBB.w_var, n_components, n_components,
				corrected_dphi_sysAAA, corrected_dphi_sysBBB,
				system_type_AAA.n_dofs_u, system_type_BBB.n_dofs_w, JxW,
				H1_coupling_const);

		H1_Coupling_Extra_Term(Me_vu, qp, system_type_AAA.v_var,
				system_type_BBB.u_var, n_components, n_components,
				corrected_dphi_sysAAA, corrected_dphi_sysBBB,
				system_type_AAA.n_dofs_v, system_type_BBB.n_dofs_u, JxW,
				H1_coupling_const);

		H1_Coupling_Extra_Term(Me_vv, qp, system_type_AAA.v_var,
				system_type_BBB.v_var, n_components, n_components,
				corrected_dphi_sysAAA, corrected_dphi_sysBBB,
				system_type_AAA.n_dofs_v, system_type_BBB.n_dofs_v, JxW,
				H1_coupling_const);

		H1_Coupling_Extra_Term(Me_vw, qp, system_type_AAA.v_var,
				system_type_BBB.w_var, n_components, n_components,
				corrected_dphi_sysAAA, corrected_dphi_sysBBB,
				system_type_AAA.n_dofs_v, system_type_BBB.n_dofs_w, JxW,
				H1_coupling_const);

		H1_Coupling_Extra_Term(Me_wu, qp, system_type_AAA.w_var,
				system_type_BBB.u_var, n_components, n_components,
				corrected_dphi_sysAAA, corrected_dphi_sysBBB,
				system_type_AAA.n_dofs_w, system_type_BBB.n_dofs_u, JxW,
				H1_coupling_const);

		H1_Coupling_Extra_Term(Me_wv, qp, system_type_AAA.w_var,
				system_type_BBB.v_var, n_components, n_components,
				corrected_dphi_sysAAA, corrected_dphi_sysBBB,
				system_type_AAA.n_dofs_w, system_type_BBB.n_dofs_v, JxW,
				H1_coupling_const);

		H1_Coupling_Extra_Term(Me_ww, qp, system_type_AAA.w_var,
				system_type_BBB.w_var, n_components, n_components,
				corrected_dphi_sysAAA, corrected_dphi_sysBBB,
				system_type_AAA.n_dofs_w, system_type_BBB.n_dofs_w, JxW,
				H1_coupling_const);
	}
	;

	void zero()
	{
		Me.zero();
	}
};

class coupled_system
{
protected:
	// Members

	// -> Equation system maps
	std::pair<std::string, libMesh::EquationSystems*> m_BIG_EquationSystem;
	std::pair<std::string, libMesh::EquationSystems*> m_R_BIG_EquationSystem;

	std::map<std::string, libMesh::EquationSystems*> m_micro_EquationSystemMap;
	std::map<std::string, libMesh::EquationSystems*> m_R_micro_EquationSystemMap;

	std::map<std::string, libMesh::EquationSystems*> m_inter_EquationSystemMap;
	std::map<std::string, libMesh::EquationSystems*> m_mediator_EquationSystemMap;

	// -> Matrix maps
	std::map<std::string, libMesh::PetscMatrix<libMesh::Number>*> m_couplingMatrixMap_mediator_micro;
	std::map<std::string, libMesh::PetscMatrix<libMesh::Number>*> m_couplingMatrixMap_mediator_BIG;
	std::map<std::string, libMesh::PetscMatrix<libMesh::Number>*> m_couplingMatrixMap_mediator_mediator;

	// -> Coordinates vector map
	std::pair<std::string, libMesh::PetscVector<libMesh::Number>* > m_coord_vect_BIG;
	std::map<std::string, libMesh::PetscVector<libMesh::Number>* > m_coord_vect_microMap;

	// -> Point locators needed for the alpha masks
	std::map<std::string, weight_parameter_function*> m_alpha_masks;

	// -> Bools for assembly of the systems
	std::map<std::string, bool> m_bHasAssembled_micro;
	std::map<std::string, bool> m_bUseH1Coupling;
	bool m_bHasAssembled_BIG;
	bool m_bHasDefinedMeshRestrictions;
	bool m_bHasDefinedCoordVector_BIG;
	bool m_bHasDefinedCoordVector_micro;

	// -> Coupling constant maps
	std::map<std::string, double> m_coupling_constantMap;
	std::map<std::string, double> m_coupling_lengthMap;

	// -> Typedefs of the destructor iterators
	typedef std::map<std::string, libMesh::EquationSystems*>::iterator EqSystem_iterator;
	typedef std::map<std::string, libMesh::PetscMatrix<libMesh::Number>*>::iterator Matrix_iterator;
	typedef std::map<std::string, libMesh::PetscVector<libMesh::Number>*>::iterator Vector_iterator;
	typedef std::map<std::string, weight_parameter_function*>::iterator alpha_mask_iterator;

	// -> LATIN solver
	carl::CoupledSolverType m_solver_type;
	std::shared_ptr<coupled_solver> m_coupled_solver;

private:
	coupled_system();

public:
	// Members

	// Constructors
	coupled_system(const libMesh::Parallel::Communicator& comm, carl::CoupledSolverType solver_type = carl::LATIN_MODIFIED_STIFFNESS) :
			m_bHasAssembled_BIG { false },
			m_bHasDefinedMeshRestrictions { false },
			m_bHasDefinedCoordVector_BIG { true },
			m_bHasDefinedCoordVector_micro { true },
			m_solver_type { solver_type }

	{
		switch (m_solver_type)
		{
			case carl::LATIN_MODIFIED_STIFFNESS:
			case carl::LATIN_ORIGINAL_STIFFNESS:
				m_coupled_solver =
						std::shared_ptr<coupled_solver>(new PETSC_LATIN_solver(comm,solver_type));
				break;
			case carl::CG:
				m_coupled_solver =
						std::shared_ptr<coupled_solver>(new PETSC_CG_solver(comm));
				break;
		}
	}
	;

	// Destructor
	~coupled_system()
	{
		this->clear();
	}

	// Methods - adding BIG / micro / inter /mediator systems
	libMesh::EquationSystems& set_BIG_EquationSystem(const std::string& name,
			libMesh::MeshBase& BIGMesh)
	{
		libMesh::EquationSystems* EqSystemPtr = new libMesh::EquationSystems(
				BIGMesh);

		m_BIG_EquationSystem.first = name;
		m_BIG_EquationSystem.second = EqSystemPtr;

		libMesh::PetscVector<libMesh::Number>* coord_vector_Ptr =
				new libMesh::PetscVector<libMesh::Number>(EqSystemPtr->comm());

		m_coord_vect_BIG.first = name;
		m_coord_vect_BIG.second = coord_vector_Ptr;

		return *EqSystemPtr;

		m_bHasDefinedCoordVector_BIG = true;
	}

	libMesh::EquationSystems& set_Restricted_BIG_EquationSystem(const std::string& name,
			libMesh::MeshBase& R_BIGMesh)
	{
		m_bHasDefinedMeshRestrictions = true;
		libMesh::EquationSystems* EqSystemPtr = new libMesh::EquationSystems(
				R_BIGMesh);

		m_R_BIG_EquationSystem.first = name;
		m_R_BIG_EquationSystem.second = EqSystemPtr;

		return *EqSystemPtr;
	}

	template<typename libMesh_MatrixType>
	libMesh::EquationSystems& add_micro_EquationSystem(const std::string& name,
			libMesh::MeshBase& microMesh)
	{
		libMesh::EquationSystems* EqSystemPtr = NULL;
		libMesh::PetscMatrix<libMesh::Number>* Matrix_mediator_micro_Ptr = NULL;
		libMesh::PetscMatrix<libMesh::Number>* Matrix_mediator_BIG_Ptr = NULL;
		libMesh::PetscMatrix<libMesh::Number>* Matrix_mediator_mediator_Ptr =
		NULL;
		libMesh::PetscVector<libMesh::Number>* coord_vector_Ptr = NULL;

		if (!m_micro_EquationSystemMap.count(name))
		{
			// Then add a new system
			EqSystemPtr = new libMesh::EquationSystems(microMesh);

			Matrix_mediator_micro_Ptr = new libMesh_MatrixType(
					microMesh.comm());
			Matrix_mediator_BIG_Ptr = new libMesh_MatrixType(microMesh.comm());
			Matrix_mediator_mediator_Ptr = new libMesh_MatrixType(
					microMesh.comm());
			coord_vector_Ptr = new libMesh::PetscVector<libMesh::Number>(microMesh.comm());

			m_micro_EquationSystemMap.insert(std::make_pair(name, EqSystemPtr));
			m_couplingMatrixMap_mediator_micro.insert(
					std::make_pair(name, Matrix_mediator_micro_Ptr));
			m_couplingMatrixMap_mediator_BIG.insert(
					std::make_pair(name, Matrix_mediator_BIG_Ptr));
			m_couplingMatrixMap_mediator_mediator.insert(
					std::make_pair(name, Matrix_mediator_mediator_Ptr));
			m_coord_vect_microMap.insert(
					std::make_pair(name, coord_vector_Ptr));

			m_bHasAssembled_micro.insert(std::make_pair(name, false));
			m_bUseH1Coupling.insert(std::make_pair(name, false));
		}
		else
		{
			// System already exits, return the system pointer
			std::cerr << " *** Warning: micro system " << name
					<< " already exists!" << std::endl;
			EqSystemPtr = m_micro_EquationSystemMap[name];
		}

		m_bHasDefinedCoordVector_micro = true;

		return *EqSystemPtr;
	}

	libMesh::EquationSystems& add_Restricted_micro_EquationSystem(const std::string& name,
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
	}


	libMesh::EquationSystems& add_inter_EquationSystem(const std::string& name,
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
	}
	;

	libMesh::EquationSystems& add_mediator_EquationSystem(
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
	}
	;

	weight_parameter_function& add_alpha_mask(const std::string& name,
			libMesh::Mesh& alpha_Mesh)
	{
		weight_parameter_function * alpha_msk_Ptr = NULL;

		if (!m_alpha_masks.count(name))
		{
			alpha_msk_Ptr = new weight_parameter_function(alpha_Mesh);
			m_alpha_masks.insert(std::make_pair(name, alpha_msk_Ptr));
		}
		else
		{
			std::cerr << " *** Warning: micro system " << name
					<< " already has an alpha mask!" << std::endl;
			alpha_msk_Ptr = m_alpha_masks[name];
		}
		return *alpha_msk_Ptr;
	}

	weight_parameter_function& get_alpha_mask(const std::string& name)
	{
		return *m_alpha_masks[name];
	}

	void set_alpha_mask_parameters(const std::string& name, double alpha_eps,
			double alpha_coupling_BIG, int subdomain_idx_BIG,
			int subdomain_idx_micro, int subdomain_idx_coupling)
	{
		(m_alpha_masks[name])->set_parameters(alpha_eps, alpha_coupling_BIG,
				subdomain_idx_BIG, subdomain_idx_micro, subdomain_idx_coupling);
	}

	void set_alpha_mask_parameters(const std::string& name,
			int subdomain_idx_BIG, int subdomain_idx_micro,
			int subdomain_idx_coupling)
	{
		(m_alpha_masks[name])->set_parameters(1E-2, 0.5, subdomain_idx_BIG,
				subdomain_idx_micro, subdomain_idx_coupling);
	}

	void set_coupling_parameters(const std::string& name,
			double coupling_constant, double coupling_length)
	{
		m_coupling_constantMap[name] = coupling_constant;
		m_coupling_lengthMap[name] = coupling_length;
	}

	void clear();

	void prepare_coupling_preallocation(
			libMesh::PetscMatrix<libMesh::Number>& coupling_matrix,
			libMesh_fe_addresses_3& row_addresses,
			libMesh_fe_addresses_3& col_addresses,
			const std::unordered_multimap<int,int>&  inter_table
			);

	void assemble_coupling_matrices(const std::string BIG_name,
			const std::string micro_name, const std::string inter_name,
			const std::string mediator_name,
			std::unordered_map<int, int>& equivalence_table_mediator_BIG,
			std::vector<std::pair<int, int> >& intersection_table_mediator_micro,
			std::unordered_multimap<int, int>& intersection_table_inter,
			double coupling_const = 1., bool using_same_mesh_mediator_A = false,
			bool bSameElemsType = true);

	void set_rigid_body_mode(libMesh::ImplicitSystem&  input_system,
			libMesh::PetscVector<libMesh::Number>* coord_vec,
			const std::string& sys_type
			);

	void set_rigid_body_modes_BIG(const std::string& sys_name);

	void set_rigid_body_modes_micro(const std::string micro_name, const std::string& sys_name);

	void assemble_coupling_elasticity_3D(const std::string BIG_name,
			const std::string micro_name, const std::string inter_name,
			const std::string mediator_name,
			std::unordered_map<int, int>& equivalence_table_mediator_BIG,
			std::vector<std::pair<int, int> >& intersection_table_mediator_micro,
			std::unordered_multimap<int, int>& intersection_table_inter,
			bool using_same_mesh_mediator_A = false,
			bool bSameElemsType = true);

	void assemble_coupling_elasticity_3D_parallel(
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
			const std::string BIG_type = "Elasticity",
			const std::string micro_type = "Elasticity",
			bool bSameElemsType = true);

	void assemble_coupling_elasticity_3D_parallel(
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
			const std::string BIG_type = "Elasticity",
			const std::string micro_type = "Elasticity",
			bool bSameElemsType = true);

	void set_restart( bool bUseRestart,
							bool bPrintRestart,
							const std::string restart_base_filename = "restart_",
							bool bPrintMatrix = false);

	void set_micro_system(
			const std::string micro_name,
			const std::string type_name,
			void fptr_assemble(		libMesh::EquationSystems& es,
									const std::string& name, weight_parameter_function& weight_mask) = nullptr );

	void set_macro_system(
			const std::string micro_name,
			const std::string type_name,
			void fptr_assemble(		libMesh::EquationSystems& es,
									const std::string& name, weight_parameter_function& weight_mask) = nullptr);

	void set_micro_system(
			const std::string micro_name,
			const std::string type_name,
			void fptr_assemble(		libMesh::EquationSystems& es,
									const std::string& name, weight_parameter_function& weight_mask, anisotropic_elasticity_tensor_cubic_sym& anisotropy_obj_input),
			anisotropic_elasticity_tensor_cubic_sym& anisotropy_obj);


	void set_LATIN_solver(	const std::string micro_name,const std::string type_name,
							double k_dA = 2.5, double k_dB = 2.5, double k_cA = 2.5, double k_cB = 2.5,
							double eps = 1E-2, int convIter = 10000, double relax = 0.8);

	void set_LATIN_solver(	const std::string micro_name, const std::string type_name,
							void fptr_BIG(		libMesh::EquationSystems& es, const std::string& name,
												weight_parameter_function& alpha_mask),
							void fptr_micro(	libMesh::EquationSystems& es, const std::string& name,
												weight_parameter_function& alpha_mask),
							double k_dA = 2.5, double k_dB = 2.5, double k_cA = 2.5, double k_cB = 2.5,
							double eps = 1E-2, int convIter = 10000, double relax = 0.8);

	void set_LATIN_solver(	const std::string micro_name, const std::string type_name,
							void fptr_BIG(		libMesh::EquationSystems& es,
												const std::string& name, weight_parameter_function& weight_mask),
							void fptr_micro(	libMesh::EquationSystems& es,
												const std::string& name, weight_parameter_function& weight_mask, anisotropic_elasticity_tensor_cubic_sym& anisotropy_obj_input),
							anisotropic_elasticity_tensor_cubic_sym& anisotropy_obj,
							double k_dA = 2.5, double k_dB = 2.5, double k_cA = 2.5, double k_cB = 2.5,
							double eps = 1E-2, int convIter = 10000, double relax = 0.8);

	void set_CG_solver(	const std::string micro_name,const std::string type_name,
							double eps_abs = 1E-50, double eps_rel = 1E-8, int convIter = 100, double div_tol = 1E4);

	void set_CG_solver(	const std::string micro_name, const std::string type_name,
							void fptr_BIG(		libMesh::EquationSystems& es, const std::string& name,
												weight_parameter_function& alpha_mask),
							void fptr_micro(	libMesh::EquationSystems& es, const std::string& name,
												weight_parameter_function& alpha_mask),
							double eps_abs = 1E-50, double eps_rel = 1E-8, int convIter = 100, double div_tol = 1E4);

	void set_CG_solver(	const std::string micro_name, const std::string type_name,
							void fptr_BIG(		libMesh::EquationSystems& es,
												const std::string& name, weight_parameter_function& weight_mask),
							void fptr_micro(	libMesh::EquationSystems& es,
												const std::string& name, weight_parameter_function& weight_mask, anisotropic_elasticity_tensor_cubic_sym& anisotropy_obj_input),
							anisotropic_elasticity_tensor_cubic_sym& anisotropy_obj,
							double eps_abs = 1E-50, double eps_rel = 1E-8, int convIter = 100, double div_tol = 1E4);

	void set_LATIN_nonlinear_solver(const std::string micro_name,
			const std::string type_name_BIG,
			const std::string type_name_micro,
			void fptr_BIG(libMesh::EquationSystems& es, const std::string& name,
					weight_parameter_function& alpha_mask),
			libMesh::EquationSystems& eq_sys_nonlinear, double k_dA = 2.5,
			double k_dB = 2.5, double k_cA = 2.5, double k_cB = 2.5,
			double eps = 1E-2, int convIter = 10000, double relax = 0.8);

	void print_convergence(const std::string& filename);

	void solve(const std::string micro_name, const std::string type_name,
			const std::string conv_name);

	void solve_LATIN_nonlinear(const std::string micro_name, const std::string type_name_micro, const std::string type_name_BIG, const std::string conv_name);

	libMesh::PetscMatrix<libMesh::Number>& get_micro_coupling_matrix(
			const std::string& name);

	libMesh::PetscMatrix<libMesh::Number>& get_BIG_coupling_matrix(
			const std::string& name);

	libMesh::PetscMatrix<libMesh::Number>& get_mediator_coupling_matrix(
			const std::string& name);

	void set_corrected_shapes(
			const std::vector<std::vector<libMesh::Real> >& lambda_weights,
			const std::vector<std::vector<libMesh::Real> >& phi_inter,
			std::vector<std::vector<libMesh::Real> >& phi_corrected);

	void set_corrected_shape_gradients(
			const std::vector<std::vector<libMesh::Real> >& lambda_weights,
			const std::vector<std::vector<libMesh::RealGradient> >& phi_inter,
			std::vector<std::vector<libMesh::RealGradient> >& phi_corrected);

	void get_lambdas(const unsigned int dim, const libMesh::FEType& fe_t,
			const libMesh::Elem* base_elem,
			const std::vector<libMesh::Point>& phys_points,
			std::vector<libMesh::Point>& ref_points,
			std::vector<std::vector<libMesh::Real> >& lambda_weights);

	void print_matrix_micro_info(const std::string& name);

	void print_matrix_BIG_info(const std::string& name);

	void print_matrix_mediator_info(const std::string& name);

	void print_matrices_matlab(const std::string& name,
			const std::string& outputRoot = "coupling");

	void set_BIG_assemble_flag(bool bFlag)
	{
		m_bHasAssembled_BIG = bFlag;
	}
	;

	void set_micro_assemble_flag(std::string name, bool bFlag)
	{
		m_bHasAssembled_micro[name] = bFlag;
	}
	;

	void use_H1_coupling(std::string name)
	{
		m_bUseH1Coupling[name] = true;
	}
	;

	void use_L2_coupling(std::string name)
	{
		m_bUseH1Coupling[name] = false;
	}
};

}
;

#endif /* ASSEMBLE_INTERSECTION_3D_H_ */
