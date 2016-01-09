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
	std::pair<std::string, libMesh::EquationSystems* > m_BIG_EquationSystem;
	std::map<std::string, libMesh::EquationSystems* > m_micro_EquationSystemMap;
	std::map<std::string, libMesh::EquationSystems* > m_inter_EquationSystemMap;
	std::map<std::string, libMesh::EquationSystems* > m_restrict_EquationSystemMap;

	//	std::vector<std::pair<int,int> > m_intersection_table_BIGmicro;
	//	std::unordered_multimap<int,int> m_intersection_table_inter;

	typedef std::map<std::string, libMesh::EquationSystems* >::iterator EqSystem_iterator;

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

	libMesh::EquationSystems& add_micro_EquationSystem(const std::string& name, libMesh::MeshBase& microMesh)
	{
		libMesh::EquationSystems* EqSystemPtr = NULL;

		if(!m_micro_EquationSystemMap.count(name))
		{
			// Then add a new system
			EqSystemPtr = new libMesh::EquationSystems(microMesh);
			m_micro_EquationSystemMap.insert (std::make_pair(name, EqSystemPtr));
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

	};

};
}

void assemble_coupling_matrix(	libMesh::EquationSystems& BIG_eq_system,
								libMesh::EquationSystems& micro_eq_system,
								libMesh::EquationSystems& inter_eq_system,
								libMesh::SparseMatrix<libMesh::Number>& couplingMatrix,
								bool bSameElemsType = true);

void assemble_coupling_matrix(	libMesh::EquationSystems& BIG_eq_system,
								libMesh::EquationSystems& micro_eq_system,
								libMesh::EquationSystems& inter_eq_system,
								std::vector<std::pair<int,int> >& intersection_table_BIGmicro,
								std::unordered_multimap<int,int>& intersection_table_inter,
								libMesh::SparseMatrix<libMesh::Number>& couplingMatrix,
								bool bSameElemsType = true);

void create_mesh_map(std::string &filename, std::unordered_map<int,int> &node_map, std::unordered_map<int,int> &element_map);

void build_mesh_map_Gmsh(std::string &filename, std::unordered_map<int,int> &node_map, std::unordered_map<int,int> &element_map);

#endif /* ASSEMBLE_INTERSECTION_3D_H_ */
