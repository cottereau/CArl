/*
 * defines_3D.h
 *
 *  Created on: Jul 31, 2015
 *      Author: breubreubreu
 *
 *         	The CGAL typedefs are defined here temporarily, to avoid "breaking"
 *      the working 2D code, they should be moved to the "CGAL_typedefs.h"
 *      header.
 *
 *      TODO: move the CGAL typedefs
 */

#ifndef DEFINES_3D_H_
#define DEFINES_3D_H_

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

//		3D convex hull
#include <CGAL/convex_hull_3.h>

// --- CGAL typedefs

////		Polyhedron domain typedefs
//typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
//typedef Polyhedron::Vertex_iterator  Vertex_iterator;
//typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, Kernel> Polyhedron_Mesh_domain;
//
////		Kernel traits
//typedef CGAL::Kernel_traits<Polyhedron_Mesh_domain>::Kernel Kt_3;
//typedef CGAL::Regular_triangulation_euclidean_traits_3< Kernel> GeomKt_3;
//
////		Vertex and cells
//typedef CGAL::Mesh_vertex_base_3<GeomKt_3,Polyhedron_Mesh_domain> 	MeshVb_3;
//
//typedef CGAL::Triangulation_cell_base_3<Kernel> Cb_3;
//typedef CGAL::Regular_triangulation_cell_base_3<GeomKt_3> RegularCb_3;
//// typedef CGAL::Triangulation_cell_base_with_info_3<CellInfo,Kernel,RegularCb_3> RegularCbInfo_3;;
//typedef CGAL::Mesh_cell_base_3<GeomKt_3,Polyhedron_Mesh_domain,RegularCb_3> MeshCbInfo_3;
//
////		Triangulation data structure
//typedef CGAL::Mesh_triangulation_3<
//									Polyhedron_Mesh_domain
////									GeomKt_3,
////									CGAL::Sequential_tag,
////									MeshVb_3,
////									MeshCbInfo_3
//								  >::type MeshTr_3;
//
//typedef CGAL::Mesh_complex_3_in_triangulation_3<MeshTr_3> CmplxMesh_3;
//
//typedef CGAL::Mesh_criteria_3<MeshTr_3> Mesh_3_criteria;

//		Polyhedron typedefs
typedef CGAL::Polyhedron_3<Kernel> 		Polyhedron;
typedef CGAL::Tetrahedron_3<Kernel>		Tetrahedron;

// Cell base
typedef CGAL::Triangulation_cell_base_with_info_3<CellInfo,Kernel>	CbInfo_3;

// Vertex base
typedef CGAL::Triangulation_vertex_base_3<Kernel>					Vb_3;

// Triangulation DS
typedef CGAL::Triangulation_data_structure_3<Vb_3,CbInfo_3>			TdsInfo_3;

// Triangulation
typedef CGAL::Delaunay_triangulation_3<Kernel,TdsInfo_3>			DT_3;

// --- Constrained triangulation iterator typedefs
typedef DT_3::Vertex 					Vertex_3;
typedef DT_3::Edge						Edge_3;
typedef DT_3::Facet 					Facet_3;
typedef DT_3::Cell	 					Cell_3;

typedef DT_3::Vertex_handle				Vertex_handle_3;
typedef DT_3::Cell_handle				Cell_handle_3;

typedef DT_3::All_vertices_iterator 	All_vertices_iterator_3;
typedef DT_3::All_edges_iterator 		All_edges_iterator_3;
typedef DT_3::All_facets_iterator 		All_facets_iterator_3;
typedef DT_3::All_cells_iterator 		All_cells_iterator_3;

typedef DT_3::Finite_vertices_iterator 	Finite_vertices_iterator_3;
typedef DT_3::Finite_edges_iterator 	Finite_edges_iterator_3;
typedef DT_3::Finite_facets_iterator 	Finite_facets_iterator_3;
typedef DT_3::Finite_cells_iterator 	Finite_cells_iterator_3;

typedef DT_3::Facet_circulator			Facet_circulator_3;
typedef DT_3::Cell_circulator			Cell_circulator_3;

typedef DT_3::Locate_type				Locate_type_3;

typedef CGAL::cpp11::result_of<Kernel::Intersect_2(Triangle_3, Triangle_3)>::type
		Triangle_3_Intersection_Variant;

#endif /* MESH_INTERSECTION_SRC_COMMON_DEFINES_3D_H_ */
