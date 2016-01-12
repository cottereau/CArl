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

#include "weak_formulations.h"

namespace carl
{

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
	std::map<std::string, libMesh::SparseMatrix<libMesh::Number>* > m_couplingMatrixMap_restrict_micro;
	std::map<std::string, libMesh::SparseMatrix<libMesh::Number>* > m_couplingMatrixMap_restrict_BIG;

	//	std::vector<std::pair<int,int> > m_intersection_table_BIGmicro;
	//	std::unordered_multimap<int,int> m_intersection_table_inter;

	typedef std::map<std::string, libMesh::EquationSystems* >::iterator EqSystem_iterator;
	typedef std::map<std::string, libMesh::SparseMatrix<libMesh::Number>* >::iterator Matrix_iterator;

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
		libMesh::SparseMatrix<libMesh::Number>* Matrix_BIG_micro_Ptr = NULL;
		libMesh::SparseMatrix<libMesh::Number>* Matrix_BIG_BIG_Ptr = NULL;

		if(!m_micro_EquationSystemMap.count(name))
		{
			// Then add a new system
			EqSystemPtr = new libMesh::EquationSystems(microMesh);

			Matrix_BIG_micro_Ptr = new libMesh_MatrixType(microMesh.comm());
			Matrix_BIG_BIG_Ptr = new libMesh_MatrixType(microMesh.comm());

			m_micro_EquationSystemMap.insert (std::make_pair(name, EqSystemPtr));
			m_couplingMatrixMap_restrict_micro.insert (std::make_pair(name, Matrix_BIG_micro_Ptr));
			m_couplingMatrixMap_restrict_BIG.insert (std::make_pair(name, Matrix_BIG_BIG_Ptr));
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

	void clear()
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

			libMesh::SparseMatrix<libMesh::Number> *Mat = toClean->second;
			Mat->clear();
			delete Mat;
			Mat = NULL;

			m_couplingMatrixMap_restrict_micro.erase(toClean);
		}

		while(!m_couplingMatrixMap_restrict_BIG.empty())
		{
			Matrix_iterator toClean = m_couplingMatrixMap_restrict_BIG.begin();

			libMesh::SparseMatrix<libMesh::Number> *Mat = toClean->second;
			Mat->clear();
			delete Mat;
			Mat = NULL;

			m_couplingMatrixMap_restrict_BIG.erase(toClean);
		}

	};

	void assemble_coupling_matrices(	const std::string BIG_name,
										const std::string micro_name,
										const std::string inter_name,
										const std::string restrict_name,
										std::vector<std::pair<int,int> >& intersection_table_restrict_micro,
										std::unordered_multimap<int,int>& intersection_table_inter,
//										libMesh::SparseMatrix<libMesh::Number>& couplingMatrix,
										bool using_same_mesh_restrict_A = false,
										bool bSameElemsType = true);

	libMesh::SparseMatrix<libMesh::Number>& get_micro_coupling_matrix(const std::string& name)
	{
		return * m_couplingMatrixMap_restrict_micro[name];
	}

	libMesh::SparseMatrix<libMesh::Number>& get_BIG_coupling_matrix(const std::string& name)
	{
		return * m_couplingMatrixMap_restrict_BIG[name];
	}

	void print_matrix_micro_info(const std::string& name)
	{
		libMesh::SparseMatrix<libMesh::Number>& CouplingTestMatrix =
							* m_couplingMatrixMap_restrict_micro[name];
		libMesh::Real accumulator = 0;
		std::cout << "| C_(R,B) i,j :" << std::endl;
		for(unsigned int iii = 0; iii < CouplingTestMatrix.m(); ++iii)
		{
			std::cout << "|    ";
			for(unsigned int jjj = 0; jjj < CouplingTestMatrix.n(); ++jjj)
			{
				std::cout <<  CouplingTestMatrix(iii,jjj) << " ";
				accumulator += CouplingTestMatrix(iii,jjj);
			}
			std::cout << std::endl;
		}
		std::cout << "|" << std::endl;
		std::cout << "| Sum( C_i,j ) = " << accumulator << std::endl << std::endl;
	}

	void print_matrix_BIG_info(const std::string& name)
	{
		libMesh::SparseMatrix<libMesh::Number>& CouplingTestMatrix =
							* m_couplingMatrixMap_restrict_BIG[name];
		libMesh::Real accumulator = 0;
		std::cout << "| C_(R,A) i,j :" << std::endl;
		for(unsigned int iii = 0; iii < CouplingTestMatrix.m(); ++iii)
		{
			std::cout << "|    ";
			for(unsigned int jjj = 0; jjj < CouplingTestMatrix.n(); ++jjj)
			{
				std::cout <<  CouplingTestMatrix(iii,jjj) << " ";
				accumulator += CouplingTestMatrix(iii,jjj);
			}
			std::cout << std::endl;
		}
		std::cout << "|" << std::endl;
		std::cout << "| Sum( C_i,j ) = " << accumulator << std::endl;
	}
};
}

// Index tables reading and mapping
void create_mesh_map(std::string &filename, std::unordered_map<int,int> &node_map, std::unordered_map<int,int> &element_map);

void build_mesh_map_Gmsh(std::string &filename, std::unordered_map<int,int> &node_map, std::unordered_map<int,int> &element_map);

#endif /* ASSEMBLE_INTERSECTION_3D_H_ */
