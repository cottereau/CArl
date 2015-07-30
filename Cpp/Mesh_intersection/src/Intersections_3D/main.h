/*
 * main.h
 *
 *  Created on: Jul 30, 2015
 *      Author: Thiago Milanetto Schlittler
 *
 *      	The CGAL typedefs are defined here temporarily, to avoid "breaking"
 *      the working 2D code, they should be moved to the "CGAL_typedefs.h"
 *      header.
 *
 *      TODO: move the CGAL typedefs
 *
 */

#ifndef INTERSECTIONS_3D_MAIN_H_
#define INTERSECTIONS_3D_MAIN_H_

//	--- CGAL kernel and common headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Intersections.h>
#include <CGAL/result_of.h>
#include "CGAL_typedefs.h"
#include "common_header.h"

// --- CGAL 3D mesh headers
// 		Base for 3D domain tetrahedrical meshes ("3D triangulations", in
//		CGAL's notation)
#include <CGAL/Mesh_triangulation_3.h>

// 		3D triangulation data structure
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_vertex_base_3.h>
#include <CGAL/Mesh_cell_base_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>

//		Meshing criteria
#include <CGAL/Mesh_criteria_3.h>

//		Polyhedron domain
#include <CGAL/Polyhedral_mesh_domain_3.h>

//		Mesh creation and refinement
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

//		IO
#include <CGAL/IO/Polyhedron_iostream.h>

// --- CGAL typedefs

//		Polyhedron domain typedefs
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, Kernel> Polyhedron_Mesh_domain;

//		Kernel traits
typedef CGAL::Kernel_traits<Polyhedron_Mesh_domain>::Kernel Kt_3;
typedef CGAL::Regular_triangulation_euclidean_traits_3< Kernel> GeomKt_3;

//		Vertex and cells
typedef CGAL::Mesh_vertex_base_3<Kt_3,Polyhedron_Mesh_domain> 	MeshVb_3;

typedef CGAL::Triangulation_cell_base_3<Kernel> Cb_3;
typedef CGAL::Triangulation_cell_base_with_info_3<FaceInfo,Kernel,Cb_3> CbInfo_3;
typedef CGAL::Regular_triangulation_cell_base_3<GeomKt_3,CbInfo_3> RegularCbInfo_3;
typedef CGAL::Mesh_cell_base_3<Kt_3,Polyhedron_Mesh_domain,RegularCbInfo_3> MeshCbInfo_3;

//		Triangulation data structure
typedef CGAL::Mesh_triangulation_3<
									Polyhedron_Mesh_domain,
									Kt_3,
									CGAL::Sequential_tag,
									MeshVb_3,
									MeshCbInfo_3
								  >::type MeshTr_3;

typedef CGAL::Mesh_complex_3_in_triangulation_3<MeshTr_3> CmplxMesh_3;

typedef CGAL::Mesh_criteria_3<MeshTr_3> Mesh_3_criteria;

#endif /* INTERSECTIONS_3D_MAIN_H_ */
