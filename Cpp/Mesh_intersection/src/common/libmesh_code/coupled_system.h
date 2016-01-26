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

namespace carl
{

void lump_matrix(		libMesh::PetscMatrix<libMesh::Number>& matrixInput,
						libMesh::PetscMatrix<libMesh::Number>& matrixOutput);

void print_matrix(libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix);

class coupled_system
{
private:
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

	typedef std::map<std::string, libMesh::EquationSystems* >::iterator EqSystem_iterator;
	typedef std::map<std::string, libMesh::PetscMatrix<libMesh::Number>* >::iterator Matrix_iterator;

	// Methods

public:
	// Members

	// Constructors
	coupled_system()
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

	void clear();

	void assemble_coupling_matrices(	const std::string BIG_name,
										const std::string micro_name,
										const std::string inter_name,
										const std::string restrict_name,
										std::unordered_map<int,int>& equivalence_table_restrict_BIG,
										std::vector<std::pair<int,int> >& intersection_table_restrict_micro,
										std::unordered_multimap<int,int>& intersection_table_inter,
										bool using_same_mesh_restrict_A = false,
										bool bSameElemsType = true);

	libMesh::PetscMatrix<libMesh::Number>& get_micro_coupling_matrix(const std::string& name);

	libMesh::PetscMatrix<libMesh::Number>& get_BIG_coupling_matrix(const std::string& name);

	libMesh::PetscMatrix<libMesh::Number>& get_restrict_coupling_matrix(const std::string& name);

	void print_matrix_micro_info(const std::string& name);

	void print_matrix_BIG_info(const std::string& name);

	void print_matrix_restrict_info(const std::string& name);
};

}

#endif /* ASSEMBLE_INTERSECTION_3D_H_ */
