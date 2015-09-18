/*
 * defines_intersection_2D.h
 *
 *  Created on: May 6, 2015
 *      Author: Thiago Milanetto Schlittler
 *
 *  	This file contains the CGAL typedefs used by the
 */

#ifndef CGAL_DEFINES_H_
#define CGAL_DEFINES_H_

//	--- CGAL kernel and common headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Intersections.h>
#include <CGAL/result_of.h>
#include "common_header.h"

// --- CGAL polygon and polyhedrons
#include <CGAL/Polygon_2.h>
#include <CGAL/Polyhedron_3.h>

// --- CGAL 3D convex hull
#include <CGAL/convex_hull_3.h>

// --- CGAL 2D triangulation headers
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

// --- CGAL 3D triangulation headers
#include <CGAL/Triangulation_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>

// Structure used to add extra information to the 2D triangles
struct VertexInfo_2
{
	int ExtIndex;
};

struct VertexInfo_3
{
	int ExtIndex;
};

struct FaceInfo
{
	int ExtIndex;
	int ExtImportedIndex;
};

struct CellInfo
{
	int ExtIndex;
	int ExtImportedIndex;
};

/*
 * 		CGAL typedefs
 */

// --- Kernel type
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

// --- Geometry typedefs
typedef Kernel::Point_2	 	Point_2;
typedef Kernel::Triangle_2	Triangle_2;
typedef Kernel::Segment_2	Segment_2;
typedef CGAL::Bbox_2		Bbox_2;

typedef Kernel::Point_3		Point_3;
typedef Kernel::Triangle_3  Triangle_3;
typedef Kernel::Segment_3	Segment_3;
typedef CGAL::Bbox_3		Bbox_3;

typedef CGAL::Polygon_2<Kernel>			Polygon_2;
typedef CGAL::Polyhedron_3<Kernel>		Polyhedron;
typedef CGAL::Tetrahedron_3<Kernel>		Tetrahedron;

// --- 2D Delaunay mesh with info
typedef CGAL::Triangulation_vertex_base_with_info_2<VertexInfo_2,Kernel> VbInfo_2;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo,Kernel> FbInfo_2;

typedef CGAL::Triangulation_data_structure_2<VbInfo_2, FbInfo_2> Tds_2;
typedef CGAL::Delaunay_triangulation_2<Kernel,Tds_2> DT_2;

// --- 2D Triangulation iterator typedefs
typedef DT_2::Vertex 					Vertex_2;
typedef DT_2::Edge						Edge_2;
typedef DT_2::Face 						Face_2;

typedef DT_2::Vertex_handle				Vertex_handle_2;
typedef DT_2::Face_handle				Face_handle_2;

typedef DT_2::Face_iterator 			Face_iterator_2;

typedef DT_2::All_vertices_iterator 	All_vertices_iterator_2;
typedef DT_2::All_faces_iterator 		All_faces_iterator_2;
typedef DT_2::Finite_faces_iterator		Finite_faces_iterator_2;
typedef DT_2::Face_circulator			Face_circulator_2;
typedef DT_2::Finite_vertices_iterator	Finite_vertices_iterator_2;

// --- 3D Delaunay mesh with info
typedef CGAL::Triangulation_vertex_base_with_info_3<VertexInfo_3,Kernel> VbInfo_3;
typedef CGAL::Triangulation_cell_base_with_info_3<CellInfo,Kernel>	CbInfo_3;

typedef CGAL::Triangulation_data_structure_3<VbInfo_3,CbInfo_3> Tds_3;
typedef CGAL::Delaunay_triangulation_3<Kernel,Tds_3> DT_3;

// --- 3D Triangulation iterator typedefs
typedef DT_3::Vertex 					Vertex_3;
typedef DT_3::Edge						Edge_3;
typedef DT_3::Facet 					Facet_3;
typedef DT_3::Cell	 					Cell_3;

typedef DT_3::Vertex_handle				Vertex_handle_3;
typedef DT_3::Cell_handle				Cell_handle_3;

typedef DT_3::Cell_iterator				Cell_iterator_3;

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

// Intersection type typedef (boost::variant)
typedef CGAL::cpp11::result_of<Kernel::Intersect_2(Triangle_2, Triangle_2)>::type
		Triangle_2_Intersection_Variant;

typedef CGAL::cpp11::result_of<Kernel::Intersect_2(Triangle_3, Triangle_3)>::type
		Triangle_3_Intersection_Variant;

// --- Geomview (for visualization, only works on Linux / MacOSX)
#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Triangulation_geomview_ostream_2.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>

#endif /* CGAL_DEFINES_H_ */
