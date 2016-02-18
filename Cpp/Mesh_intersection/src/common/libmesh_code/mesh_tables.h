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
#include "common_header_libmesh.h"

namespace carl
{

void set_weight_function_domain_idx(	std::string &filename,
										int& domain_Idx_BIG,
										int& nb_of_domain_Idx_micro,
										std::vector<int>& domain_Idx_micro,
										std::vector<int>& domain_Idx_coupling
										);

void create_mesh_map(const std::string &filename, std::unordered_map<int,int> &node_map, std::unordered_map<int,int> &element_map);

void build_mesh_map_Gmsh(const std::string &filename, std::unordered_map<int,int> &node_map, std::unordered_map<int,int> &element_map);

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

void set_mesh_Gmsh(	libMesh::Mesh& mesh, const std::string& mesh_file,
					std::unordered_map<int,int>& mesh_NodeMap, std::unordered_map<int,int>& mesh_ElemMap);

void set_mesh_Gmsh(	libMesh::Mesh& mesh, const std::string& mesh_file);

};





#endif /* COMMON_LIBMESH_CODE_MESH_TABLES_H_ */
