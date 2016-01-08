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
