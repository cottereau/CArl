/*
 * defines_intersection_2D.h
 *
 *  Created on: May 6, 2015
 *      Author: Thiago Milanetto Schlittler
 *
 *  	This file contains the CGAL typedefs used by the
 */

#ifndef CGAL_DEFINES_2_H_
#define CGAL_DEFINES_2_H_

//	--- CGAL kernel and common headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Intersections.h>
#include <CGAL/result_of.h>
#include "common_header.h"

// --- CGAL 2D triangulation headers
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

// Structure used to add extra information to the 2D triangles
struct FaceInfo
{
	int ExtIndex;
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
//
// --- Constrained Delaunay mesh with info
typedef CGAL::Triangulation_vertex_base_2<Kernel> Vb_2;

typedef CGAL::Delaunay_mesh_face_base_2<Kernel> Fb_2;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo,Kernel,Fb_2> FbInfo_2;

typedef CGAL::Triangulation_data_structure_2<Vb_2, FbInfo_2> TdsInfo_2;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel,TdsInfo_2> CDT;

typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria_2;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria_2> Mesher_2;

////// --- Delaunay mesh with info
//typedef CGAL::Triangulation_vertex_base_2<Kernel> Vb;
//typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo,Kernel> FbInfo;
//
//typedef CGAL::Triangulation_data_structure_2<Vb, FbInfo> Tds;
//typedef CGAL::Delaunay_triangulation_2<Kernel,Tds> CDT;
//
// --- Constrained triangulation iterator typedefs
typedef CDT::Face 						Face_2;
typedef CDT::Vertex_handle				Vertex_handle_2;
typedef CDT::Face_handle				Face_handle_2;
typedef CDT::Face_iterator 				Face_iterator_2;
typedef CDT::All_faces_iterator 		All_Face_iterator_2;
typedef CDT::Finite_faces_iterator		Finite_face_iterator_2;
typedef CDT::Face_circulator			Face_circulator_2;
typedef CDT::Finite_vertices_iterator	Finite_vertices_iterator_2;


// Intersection type typedef (boost::variant)
typedef CGAL::cpp11::result_of<Kernel::Intersect_2(Triangle_2, Triangle_2)>::type
		Triangle_2_Intersection_Variant;

// Boost Random defines
extern boost::random::lagged_fibonacci607 m_rng;

#endif /* CGAL_DEFINES_2_H_ */
