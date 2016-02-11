/*
 * assemble_intersection_3D.h
 *
 *  Created on: Nov 30, 2015
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef ASSEMBLE_INTERSECTION_3D_H_
#define ASSEMBLE_INTERSECTION_3D_H_

#include "common_header.h"
#include "common_header_libmesh.h"
#include "common_functions.h"

#include "weak_formulations.h"
#include "LATIN_solver.h"
#include "PETSC_matrix_operations.h"
#include "weight_parameter_function.h"

const bool MASTER_bPerfLog_assemble_coupling = false;

namespace carl
{

class coupled_system
{
protected:
	// Members

	// -> Equation system maps
	std::pair<std::string, libMesh::EquationSystems* > m_BIG_EquationSystem;
	std::map<std::string, libMesh::EquationSystems* > m_micro_EquationSystemMap;
	std::map<std::string, libMesh::EquationSystems* > m_inter_EquationSystemMap;
	std::map<std::string, libMesh::EquationSystems* > m_restrict_EquationSystemMap;

	// -> Matrix maps
	std::map<std::string, libMesh::PetscMatrix<libMesh::Number>* > m_couplingMatrixMap_restrict_micro;
	std::map<std::string, libMesh::PetscMatrix<libMesh::Number>* > m_couplingMatrixMap_restrict_BIG;
	std::map<std::string, libMesh::PetscMatrix<libMesh::Number>* > m_couplingMatrixMap_restrict_restrict;

	// -> Point locators needed for the alpha masks
	std::map<std::string, weight_parameter_function* > m_alpha_masks;

	// -> Bools for assembly of the systems
	std::map<std::string, bool > m_bHasAssembled_micro;
	bool m_bHasAssembled_BIG;

	typedef std::map<std::string, libMesh::EquationSystems* >::iterator EqSystem_iterator;
	typedef std::map<std::string, libMesh::PetscMatrix<libMesh::Number>* >::iterator Matrix_iterator;
	typedef std::map<std::string, weight_parameter_function* >::iterator alpha_mask_iterator;

	// -> LATIN solver
	PETSC_LATIN_solver m_LATIN_solver;


private:
	coupled_system();

public:
	// Members

	// Constructors
	coupled_system(const libMesh::Parallel::Communicator& comm) :
		m_LATIN_solver { PETSC_LATIN_solver(comm) },
		m_bHasAssembled_BIG { false }
	{
	};

	// Destructor
	~coupled_system()
	{
		this->clear();
	}

	// Methods - adding BIG / micro / inter /restrict systems
	libMesh::EquationSystems& set_BIG_EquationSystem(const std::string& name, libMesh::MeshBase& BIGMesh)
	{
		libMesh::EquationSystems* EqSystemPtr = new libMesh::EquationSystems(BIGMesh);

		m_BIG_EquationSystem.first = name;
		m_BIG_EquationSystem.second = EqSystemPtr;

		return *EqSystemPtr;
	}

	template <typename libMesh_MatrixType>
	libMesh::EquationSystems& add_micro_EquationSystem(const std::string& name, libMesh::MeshBase& microMesh)
	{
		libMesh::EquationSystems* EqSystemPtr = NULL;
		libMesh::PetscMatrix<libMesh::Number>* Matrix_restrict_micro_Ptr = NULL;
		libMesh::PetscMatrix<libMesh::Number>* Matrix_restrict_BIG_Ptr = NULL;
		libMesh::PetscMatrix<libMesh::Number>* Matrix_restrict_restrict_Ptr = NULL;

		if(!m_micro_EquationSystemMap.count(name))
		{
			// Then add a new system
			EqSystemPtr = new libMesh::EquationSystems(microMesh);

			Matrix_restrict_micro_Ptr = new libMesh_MatrixType(microMesh.comm());
			Matrix_restrict_BIG_Ptr = new libMesh_MatrixType(microMesh.comm());
			Matrix_restrict_restrict_Ptr = new libMesh_MatrixType(microMesh.comm());

			m_micro_EquationSystemMap.insert (std::make_pair(name, EqSystemPtr));
			m_couplingMatrixMap_restrict_micro.insert (std::make_pair(name, Matrix_restrict_micro_Ptr));
			m_couplingMatrixMap_restrict_BIG.insert (std::make_pair(name, Matrix_restrict_BIG_Ptr));
			m_couplingMatrixMap_restrict_restrict.insert (std::make_pair(name, Matrix_restrict_restrict_Ptr));
			m_bHasAssembled_micro.insert (std::make_pair(name,false));
		}
		else
		{
			// System already exits, return the system pointer
			std::cerr << " *** Warning: micro system " << name << " already exists!" << std::endl;
			EqSystemPtr = m_micro_EquationSystemMap[name];
		}

		return *EqSystemPtr;
	};

	libMesh::EquationSystems& add_inter_EquationSystem(const std::string& name, libMesh::MeshBase& interMesh)
	{
		libMesh::EquationSystems* EqSystemPtr = NULL;

		if(!m_inter_EquationSystemMap.count(name))
		{
			// Then add a new system
			EqSystemPtr = new libMesh::EquationSystems(interMesh);
			m_inter_EquationSystemMap.insert (std::make_pair(name, EqSystemPtr));
		}
		else
		{
			// System already exits, return the system pointer
			std::cerr << " *** Warning: inter system " << name << " already exists!" << std::endl;
			EqSystemPtr = m_inter_EquationSystemMap[name];
		}

		return *EqSystemPtr;
	};

	libMesh::EquationSystems& add_restrict_EquationSystem(const std::string& name, libMesh::MeshBase& restrictMesh)
	{
		libMesh::EquationSystems* EqSystemPtr = NULL;

		if(!m_restrict_EquationSystemMap.count(name))
		{
			// Then add a new system
			EqSystemPtr = new libMesh::EquationSystems(restrictMesh);
			m_restrict_EquationSystemMap.insert (std::make_pair(name, EqSystemPtr));
		}
		else
		{
			// System already exits, return the system pointer
			std::cerr << " *** Warning: restrict system " << name << " already exists!" << std::endl;
			EqSystemPtr = m_restrict_EquationSystemMap[name];
		}

		return *EqSystemPtr;
	};

	weight_parameter_function& add_alpha_mask(const std::string& name, libMesh::Mesh& alpha_Mesh)
	{
		weight_parameter_function * alpha_msk_Ptr = NULL;

		if(!m_alpha_masks.count(name))
		{
			alpha_msk_Ptr = new weight_parameter_function(alpha_Mesh);
			m_alpha_masks.insert (std::make_pair(name, alpha_msk_Ptr));
		}
		else
		{
			std::cerr << " *** Warning: micro system " << name << " already has an alpha mask!" << std::endl;
			alpha_msk_Ptr = m_alpha_masks[name];
		}
		return *alpha_msk_Ptr;
	}

	weight_parameter_function& get_alpha_mask(const std::string& name)
	{
		return * m_alpha_masks[name];
	}

	void set_alpha_mask_parameters(const std::string& name, double alpha_eps, double alpha_coupling_BIG,
			int subdomain_idx_BIG, int subdomain_idx_micro, int subdomain_idx_coupling)
	{
		(m_alpha_masks[name])->set_parameters(alpha_eps,alpha_coupling_BIG,subdomain_idx_BIG,subdomain_idx_micro,subdomain_idx_coupling);
	}

	void set_alpha_mask_parameters(const std::string& name, int subdomain_idx_BIG, int subdomain_idx_micro, int subdomain_idx_coupling)
	{
		(m_alpha_masks[name])->set_parameters(1E-2,0.5,subdomain_idx_BIG,subdomain_idx_micro,subdomain_idx_coupling);
	}

	void clear();

	void assemble_coupling_matrices(	const std::string BIG_name,
										const std::string micro_name,
										const std::string inter_name,
										const std::string restrict_name,
										std::unordered_map<int,int>& equivalence_table_restrict_BIG,
										std::vector<std::pair<int,int> >& intersection_table_restrict_micro,
										std::unordered_multimap<int,int>& intersection_table_inter,
										double coupling_const = 1.,
										bool using_same_mesh_restrict_A = false,
										bool bSameElemsType = true);

	void assemble_coupling_elasticity_3D(	const std::string BIG_name,
											const std::string micro_name,
											const std::string inter_name,
											const std::string restrict_name,
											std::unordered_map<int,int>& equivalence_table_restrict_BIG,
											std::vector<std::pair<int,int> >& intersection_table_restrict_micro,
											std::unordered_multimap<int,int>& intersection_table_inter,
											double coupling_const = 1.,
											bool using_same_mesh_restrict_A = false,
											bool bSameElemsType = true);

	void assemble_coupling_matrices_old(	const std::string BIG_name,
										const std::string micro_name,
										const std::string inter_name,
										const std::string restrict_name,
										std::unordered_map<int,int>& equivalence_table_restrict_BIG,
										std::vector<std::pair<int,int> >& intersection_table_restrict_micro,
										std::unordered_multimap<int,int>& intersection_table_inter,
										double coupling_const = 1.,
										bool using_same_mesh_restrict_A = false,
										bool bSameElemsType = true);

	void assemble_coupling_elasticity_3D_old(	const std::string BIG_name,
											const std::string micro_name,
											const std::string inter_name,
											const std::string restrict_name,
											std::unordered_map<int,int>& equivalence_table_restrict_BIG,
											std::vector<std::pair<int,int> >& intersection_table_restrict_micro,
											std::unordered_multimap<int,int>& intersection_table_inter,
											double coupling_const = 1.,
											bool using_same_mesh_restrict_A = false,
											bool bSameElemsType = true);

	void set_LATIN_solver(const std::string micro_name, const std::string type_name);

	void set_LATIN_solver(	const std::string micro_name, const std::string type_name,
							void fptr_BIG(		libMesh::EquationSystems& es,
												const std::string& name, weight_parameter_function& alpha_mask),
							void fptr_micro(	libMesh::EquationSystems& es,
												const std::string& name, weight_parameter_function& alpha_mask));

	void solve_LATIN(const std::string micro_name, const std::string type_name);

	libMesh::PetscMatrix<libMesh::Number>& get_micro_coupling_matrix(const std::string& name);

	libMesh::PetscMatrix<libMesh::Number>& get_BIG_coupling_matrix(const std::string& name);

	libMesh::PetscMatrix<libMesh::Number>& get_restrict_coupling_matrix(const std::string& name);

	void set_corrected_shapes(	const std::vector<std::vector<libMesh::Real> >& 	lambda_weights,
								const std::vector<std::vector<libMesh::Real> >& 	phi_inter,
								std::vector<std::vector<libMesh::Real> >& 			phi_corrected);

	void get_lambdas(	const unsigned int 							dim,
						const libMesh::FEType& 						fe_t,
						const libMesh::Elem* 						base_elem,
						const std::vector<libMesh::Point>& 			phys_points,
						std::vector<libMesh::Point>& 				ref_points,
						std::vector<std::vector<libMesh::Real> >& 	lambda_weights);

	void print_matrix_micro_info(const std::string& name);

	void print_matrix_BIG_info(const std::string& name);

	void print_matrix_restrict_info(const std::string& name);

	void print_matrices_matlab(const std::string& name, const std::string& outputRoot = "coupling");

	void set_BIG_assemble_flag(bool bFlag)
	{
		m_bHasAssembled_BIG = bFlag;
	};

	void set_micro_assemble_flag(std::string name, bool bFlag)
	{
		m_bHasAssembled_micro[name] = bFlag;
	};
};

}

#endif /* ASSEMBLE_INTERSECTION_3D_H_ */
