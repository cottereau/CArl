/*
 * mesh_tables.h
 *
 *  Created on: Jan 26, 2016
 *      Author: Thiago Milanetto Schlittler
 *
 *      File containing the tools used to read and generate the mesh tables
 */

#ifndef COMMON_LIBMESH_CODE_MESH_TABLES_H_
#define COMMON_LIBMESH_CODE_MESH_TABLES_H_

#include "common_header.h"
#include "common_functions.h"
#include "common_header_libmesh.h"
#include "mpi_carl_tools.h"
#include <petscsys.h>

namespace carl
{

void set_weight_function_domain_idx(	std::string &filename,
										int& domain_Idx_BIG,
										int& nb_of_domain_Idx_micro,
										std::vector<int>& domain_Idx_micro,
										std::vector<int>& domain_Idx_coupling
										);

void set_mesh_Gmsh(
		libMesh::Mesh& mesh,
		const std::string& mesh_file,
		std::unordered_map<int,int>& mesh_NodeMap,
		std::unordered_map<int,int>& mesh_ElemMap
		);

void set_mesh_Gmsh(
		libMesh::Mesh& mesh,
		const std::string& mesh_file,
		std::unordered_map<int,int> &node_gmsh_to_libmesh_map,
		std::unordered_map<int,int> &node_libmesh_to_gmsh_map,
		std::unordered_map<int,int> &element_gmsh_to_libmesh_map,
		std::unordered_map<int,int> &element_libmesh_to_gmsh_map
		);

void set_mesh_Gmsh(
		libMesh::Mesh& mesh,
		const std::string& mesh_file
		);

void create_mesh_map(
		const std::string &filename,
		std::unordered_map<int,int> &node_map,
		std::unordered_map<int,int> &element_map
		);

void create_mesh_map(
		const std::string &filename,
		std::unordered_map<int,int> &node_gmsh_to_libmesh_map,
		std::unordered_map<int,int> &node_libmesh_to_gmsh_map,
		std::unordered_map<int,int> &element_gmsh_to_libmesh_map,
		std::unordered_map<int,int> &element_libmesh_to_gmsh_map
		);

void create_mesh_map(
		const std::string &filename,
		std::unordered_map<int,int> &node_map,
		std::unordered_map<int,int> &element_map,
		const libMesh::Parallel::Communicator& WorldComm
		);

void create_mesh_map(
		const std::string &filename,
		std::unordered_map<int,int> &node_gmsh_to_libmesh_map,
		std::unordered_map<int,int> &node_libmesh_to_gmsh_map,
		std::unordered_map<int,int> &element_gmsh_to_libmesh_map,
		std::unordered_map<int,int> &element_libmesh_to_gmsh_map,
		const libMesh::Parallel::Communicator& WorldComm
		);

void build_mesh_map_Gmsh(
		const std::string &filename,
		std::unordered_map<int,int> &node_map,
		std::unordered_map<int,int> &element_map
		);

void build_mesh_map_Gmsh(
		const std::string &filename,
		std::unordered_map<int,int> &node_gmsh_to_libmesh_map,
		std::unordered_map<int,int> &node_libmesh_to_gmsh_map,
		std::unordered_map<int,int> &element_gmsh_to_libmesh_map,
		std::unordered_map<int,int> &element_libmesh_to_gmsh_map
		);
;
void build_intersection_and_restriction_tables(
		const libMesh::Parallel::Communicator& WorldComm,
		const std::string& intersection_full_table_Filename,
		const std::string& equivalence_table_A_Filename,
		const std::string& equivalence_table_B_Filename,
		std::vector<carl::IntersectionData>& intersection_full_table,
		std::unordered_map<int,int>& equivalence_table_A_to_R_A,
		std::unordered_map<int,int>& equivalence_table_B_to_R_B,
		std::unordered_map<int,int>& equivalence_table_R_A_to_A,
		std::unordered_map<int,int>& equivalence_table_R_B_to_B
		);

void generate_intersection_tables_partial(	std::string& intersection_table_restrict_B_Filename,
		std::string& intersection_table_I_Filename,
		std::unordered_map<int,int>& mesh_restrict_ElemMap,
		std::unordered_map<int,int>& mesh_micro_ElemMap,
		std::unordered_map<int,int>& mesh_inter_ElemMap,
		std::vector<std::pair<int,int> >& intersection_table_restrict_B,
		std::unordered_multimap<int,int>& intersection_table_I
		);

void generate_intersection_tables_full(		std::string& equivalence_table_restrict_A_Filename,
		std::string& intersection_table_restrict_B_Filename,
		std::string& intersection_table_I_Filename,
		std::unordered_map<int,int>& mesh_restrict_ElemMap,
		std::unordered_map<int,int>& mesh_micro_ElemMap,
		std::unordered_map<int,int>& mesh_BIG_ElemMap,
		std::unordered_map<int,int>& mesh_inter_ElemMap,
		std::unordered_map<int,int>& equivalence_table_restrict_A,
		std::vector<std::pair<int,int> >& intersection_table_restrict_B,
		std::unordered_multimap<int,int>& intersection_table_I
		);

void convert_intersection_table_to_mediator(
		const libMesh::Parallel::Communicator& WorldComm,
		const libMesh::Mesh& mesh_intersection,
		const std::unordered_map<int,int>& mesh_intersection_ElemMap,
		const std::vector<carl::IntersectionData>& intersection_full_table,
		const std::unordered_map<int,int>& equivalence_table_A_to_R_A,
		const std::unordered_map<int,int>& equivalence_table_B_to_R_B,
		const std::unordered_map<int,int>& equivalence_table_mediator,
		std::vector<carl::IntersectionData>& local_intersection_table
		);

void set_equivalence_tables(
		const libMesh::Parallel::Communicator& WorldComm,
		const std::string& equivalence_table_A_Filename,
		const std::string& equivalence_table_B_Filename,

		const std::unordered_map<int,int>& mesh_A_ElemMap,
		const std::unordered_map<int,int>& mesh_RA_ElemMap,

		const std::unordered_map<int,int>& mesh_B_ElemMap,
		const std::unordered_map<int,int>& mesh_RB_ElemMap,

		std::unordered_map<int,int>& equivalence_table_A_to_R_A,
		std::unordered_map<int,int>& equivalence_table_B_to_R_B,
		std::unordered_map<int,int>& equivalence_table_R_A_to_A,
		std::unordered_map<int,int>& equivalence_table_R_B_to_B );

void set_restricted_intersection_pairs_table(
		const std::unordered_map<int,std::pair<int,int> >&  full_intersection_pairs_map,
		const std::unordered_map<int,int>& equivalence_table_A_to_R_A,
		const std::unordered_map<int,int>& equivalence_table_B_to_R_B,
		std::unordered_map<int,std::pair<int,int> >& full_intersection_restricted_pairs_map);

void set_full_intersection_tables(
		const libMesh::Parallel::Communicator& WorldComm,
		const std::string& intersection_full_table_Filename,

		const std::unordered_map<int,int>& mesh_A_ElemMap,
		const std::unordered_map<int,int>& mesh_B_ElemMap,

		std::unordered_map<int,std::pair<int,int> >& full_intersection_pairs_map,
		std::unordered_map<int,int>& full_intersection_meshI_to_inter_map);

void set_intersection_tables(
		const libMesh::Parallel::Communicator& WorldComm,
		const libMesh::Mesh& mesh_intersection,
		const std::string& intersection_full_table_Filename,
		const std::string& equivalence_table_A_Filename,
		const std::string& equivalence_table_B_Filename,
		const std::unordered_map<int,int>& mesh_I_libmesh_to_table_ElemMap,

		const std::unordered_map<int,int>& equivalence_table_A_to_R_A,
		const std::unordered_map<int,int>& equivalence_table_B_to_R_B,

		const std::unordered_map<int,int>& mesh_A_ElemMap,
		const std::unordered_map<int,int>& mesh_B_ElemMap,

		std::unordered_map<int,std::pair<int,int> >& full_intersection_pairs_map,
		std::unordered_map<int,std::pair<int,int> >& full_intersection_restricted_pairs_map,
		std::unordered_map<int,int>& local_intersection_meshI_to_inter_map

		);


//void convert_intersection_table_to_mediator(
//		const libMesh::Mesh& mesh_mediator,
//		const std::string& equivalence_table_mediator_original_Filename,
//		const std::string& intersection_pairs_mediator_other_system_Filename,
//		const std::unordered_map<int,int>& mesh_mediator_ElemMap,
//		const std::unordered_map<int,int>& mesh_other_system_ElemMap,
//		std::vector<std::pair<int,int> >& local_corrected_intersection_pairs );
};





#endif /* COMMON_LIBMESH_CODE_MESH_TABLES_H_ */
