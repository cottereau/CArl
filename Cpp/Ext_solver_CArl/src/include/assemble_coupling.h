/*
 * assemble_coupling.h
 *
 *  Created on: Apr 13, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef ASSEMBLE_COUPLING_H_
#define ASSEMBLE_COUPLING_H_

#include "carl_headers.h"

#include "mpi_carl_tools.h"
#include "weak_formulations.h"
#include "PETSC_matrix_operations.h"

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
	const libMesh::Elem* elem;

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
		const libMesh::Elem* elem_dummy = mesh.elem(idx);
		elem = mesh.elem(idx);
		dof_map.dof_indices(elem_dummy, dof_indices);
		dof_map.dof_indices(elem_dummy, dof_indices_u, u_var);
		dof_map.dof_indices(elem_dummy, dof_indices_v, v_var);
		dof_map.dof_indices(elem_dummy, dof_indices_w, w_var);
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

class assemble_coupling_matrices
{
protected:
	// Members

	// Communicator
	const libMesh::Parallel::Communicator * m_comm;

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

	// -> Bools for assembly of the system
	std::map<std::string, bool> m_bUseH1Coupling;
	std::map<std::string, bool> m_bHasAssembled_micro;
	bool m_bHasAssembled_BIG;
	bool m_bUseNullSpace_BIG;
	bool m_bHasDefinedMeshRestrictions;

	// -> Coupling constant maps
	std::map<std::string, double> m_coupling_rigidityMap;
	std::map<std::string, double> m_coupling_widthMap;

	// -> Typedefs of the destructor iterators
	typedef std::map<std::string, libMesh::EquationSystems*>::iterator EqSystem_iterator;
	typedef std::map<std::string, libMesh::PetscMatrix<libMesh::Number>*>::iterator Matrix_iterator;
	typedef std::map<std::string, libMesh::PetscVector<libMesh::Number>*>::iterator Vector_iterator;

private:
	assemble_coupling_matrices();

public:
	// Members

	// Constructors
	assemble_coupling_matrices(const libMesh::Parallel::Communicator& comm) :
			m_comm { &comm },
			m_bHasAssembled_BIG { false },
			m_bHasDefinedMeshRestrictions { false }
	{
	};

	// Destructor
	~assemble_coupling_matrices()
	{
		this->clear();
	}

	void clear();

	// Methods implemented in assemble_coupling.cpp
	// -> Adding BIG / micro / inter /mediator systems
	libMesh::EquationSystems& set_BIG_EquationSystem(const std::string& name,
			libMesh::MeshBase& BIGMesh);

	libMesh::EquationSystems& set_Restricted_BIG_EquationSystem(const std::string& name,
			libMesh::MeshBase& R_BIGMesh);

	template<typename libMesh_MatrixType>
	libMesh::EquationSystems& add_micro_EquationSystem(const std::string& name,
			libMesh::MeshBase& microMesh)
	{
		libMesh::EquationSystems* EqSystemPtr = NULL;
		libMesh::PetscMatrix<libMesh::Number>* Matrix_mediator_micro_Ptr = NULL;
		libMesh::PetscMatrix<libMesh::Number>* Matrix_mediator_BIG_Ptr = NULL;
		libMesh::PetscMatrix<libMesh::Number>* Matrix_mediator_mediator_Ptr =
		NULL;

		if (!m_micro_EquationSystemMap.count(name))
		{
			// Then add a new system
			EqSystemPtr = new libMesh::EquationSystems(microMesh);

			Matrix_mediator_micro_Ptr = new libMesh_MatrixType(
					microMesh.comm());
			Matrix_mediator_BIG_Ptr = new libMesh_MatrixType(microMesh.comm());
			Matrix_mediator_mediator_Ptr = new libMesh_MatrixType(
					microMesh.comm());

			m_micro_EquationSystemMap.insert(std::make_pair(name, EqSystemPtr));
			m_couplingMatrixMap_mediator_micro.insert(
					std::make_pair(name, Matrix_mediator_micro_Ptr));

			m_couplingMatrixMap_mediator_BIG.insert(
					std::make_pair(name, Matrix_mediator_BIG_Ptr));

			m_couplingMatrixMap_mediator_mediator.insert(
					std::make_pair(name, Matrix_mediator_mediator_Ptr));

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

		return *EqSystemPtr;
	}

	libMesh::EquationSystems& add_Restricted_micro_EquationSystem(const std::string& name,
			libMesh::MeshBase& R_microMesh);

	libMesh::EquationSystems& add_inter_EquationSystem(const std::string& name,
			libMesh::MeshBase& interMesh);

	libMesh::EquationSystems& add_mediator_EquationSystem(
			const std::string& name, libMesh::MeshBase& mediatorMesh);

	/// Set physical parameters
	void set_coupling_parameters(const std::string& name,
			double coupling_rigidity, double coupling_width);

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

	libMesh::PetscMatrix<libMesh::Number>& get_micro_coupling_matrix(
			const std::string& name);

	libMesh::PetscMatrix<libMesh::Number>& get_BIG_coupling_matrix(
			const std::string& name);

	libMesh::PetscMatrix<libMesh::Number>& get_mediator_coupling_matrix(
			const std::string& name);
			
	void print_matrix_micro_info(const std::string& name);

	void print_matrix_BIG_info(const std::string& name);

	void print_matrix_mediator_info(const std::string& name);

	void print_matrices_matlab(const std::string& name,
			const std::string& outputRoot = "coupling");

	void print_PETSC_matrices(const std::string& name, const std::string& outputRoot = "coupling");
	void use_H1_coupling(std::string name);

	void use_L2_coupling(std::string name);

	// Methods implemented in assemble_functions_coupling.cpp
	void prepare_coupling_preallocation(
			libMesh::PetscMatrix<libMesh::Number>& coupling_matrix,
			libMesh_fe_addresses_3& row_addresses,
			libMesh_fe_addresses_3& col_addresses,
			const std::unordered_multimap<int,int>&  inter_table
			);

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

	void check_coupling_construction_3D_parallel(
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
};
};

#endif /* ASSEMBLE_COUPLING_H_ */
