/*
 * intersection_tools.cpp
 *
 *  Created on: Apr 17, 2016
 *      Author: Thiago Milanetto Schlittler
 *
 *  \brief **STAT/DYN-DI/DYN-CG** Intersection step
 */

#include "intersection_tools.h"

namespace carl
{

void Intersection_Tools::set_element_indexes()
{
  m_TET_tetrahedrons.resize(1,std::vector<unsigned int>(4,0));
  m_TET_triangles.resize(4,std::vector<unsigned int>(3,0));
  m_TET_edges.resize(6,std::vector<unsigned int>(2,0));

  m_HEX_tetrahedrons.resize(5,std::vector<unsigned int>(4,0));
  m_HEX_triangles.resize(12,std::vector<unsigned int>(3,0));
  m_HEX_edges.resize(12,std::vector<unsigned int>(2,0));

  // Set up tetra tetrahedrons (...)
  m_TET_tetrahedrons[0] = {0, 1, 2, 3};

  // Set up tetra triangles
  m_TET_triangles[0] = {0, 1, 2};
  m_TET_triangles[1] = {1, 2, 3};
  m_TET_triangles[2] = {0, 1, 3};
  m_TET_triangles[3] = {0, 2, 3};

  // Set up tetra edges
  m_TET_edges[0] = {0, 1};
  m_TET_edges[1] = {0, 2};
  m_TET_edges[2] = {0, 3};
  m_TET_edges[3] = {1, 3};
  m_TET_edges[4] = {3, 2};
  m_TET_edges[5] = {2, 1};

  // Set up hex tetrahedrons
  m_HEX_tetrahedrons[0] = {0, 1, 3, 4};
  m_HEX_tetrahedrons[1] = {5, 1, 4, 6};
  m_HEX_tetrahedrons[2] = {7, 3, 4, 6};
  m_HEX_tetrahedrons[3] = {2, 6, 1, 3};

  m_HEX_tetrahedrons[4] = {1, 4, 6, 3};

  // Set up hex triangles
  m_HEX_triangles[0] = {0, 1, 2};
  m_HEX_triangles[1] = {0, 3, 2};

  m_HEX_triangles[2] = {0, 1, 5};
  m_HEX_triangles[3] = {0, 4, 5};

  m_HEX_triangles[4] = {0, 3, 7};
  m_HEX_triangles[5] = {0, 4, 7};

  m_HEX_triangles[6] = {6, 5, 4};
  m_HEX_triangles[7] = {6, 7, 4};

  m_HEX_triangles[8]  = {6, 5, 1};
  m_HEX_triangles[9]  = {6, 2, 1};

  m_HEX_triangles[10] = {6, 7, 3};
  m_HEX_triangles[11] = {6, 2, 3};

  // Set up hex edges
  m_HEX_edges[0] = {0, 1};
  m_HEX_edges[1] = {1, 2};
  m_HEX_edges[2] = {2, 3};
  m_HEX_edges[3] = {3, 0};

  m_HEX_edges[4] = {4, 5};
  m_HEX_edges[5] = {5, 6};
  m_HEX_edges[6] = {6, 7};
  m_HEX_edges[7] = {7, 4};

  m_HEX_edges[8]  = {0, 4};
  m_HEX_edges[9]  = {1, 5};
  m_HEX_edges[10] = {2, 6};
  m_HEX_edges[11] = {3, 7};
}

bool Intersection_Tools::elements_do_intersect(
    std::vector<ExactPoint_3> & elem_C_points,
    std::vector<std::vector<unsigned int> > & elem_C_tetras,
    std::vector<std::vector<unsigned int> > & elem_C_triangles,
    std::vector<ExactPoint_3> & elem_D_points,
    std::vector<std::vector<unsigned int> > & elem_D_tetras,
    std::vector<std::vector<unsigned int> > & elem_D_triangles)
{
  bool bElemIntersect = false;

  // Test intersections between C's tetrahedrons and D's triangles
  for(unsigned int jjj = 0; jjj < elem_D_triangles.size(); ++jjj)
  {
    std::vector<unsigned int> & work_triangle = elem_D_triangles[jjj];
    m_test_triangle = ExactTriangle_3(elem_D_points[work_triangle[0]],
        elem_D_points[work_triangle[1]],
        elem_D_points[work_triangle[2]]);

    for(unsigned int iii = 0; iii < elem_C_tetras.size(); ++iii)
    {
      std::vector<unsigned int> & work_tetra = elem_C_tetras[iii];
      m_test_tetra = ExactTetrahedron(  elem_C_points[work_tetra[0]],
          elem_C_points[work_tetra[1]],
          elem_C_points[work_tetra[2]],
          elem_C_points[work_tetra[3]]);

      bElemIntersect = CGAL::do_intersect(m_test_triangle,m_test_tetra);

      if(bElemIntersect)
      {
        // Found intersection!
        break;
      }
    }

    if(bElemIntersect)
    {
      // Found intersection!
      break;
    }
  }

  if(!bElemIntersect)
  {
    // Test intersections between D's tetrahedrons and C's triangles
    for(unsigned int jjj = 0; jjj < elem_C_triangles.size(); ++jjj)
    {
      std::vector<unsigned int> & work_triangle = elem_C_triangles[jjj];
      m_test_triangle = ExactTriangle_3(elem_C_points[work_triangle[0]],
          elem_C_points[work_triangle[1]],
          elem_C_points[work_triangle[2]]);

      for(unsigned int iii = 0; iii < elem_D_tetras.size(); ++iii)
      {
        std::vector<unsigned int> & work_tetra = elem_D_tetras.at(iii);
        m_test_tetra = ExactTetrahedron(  elem_D_points[work_tetra[0]],
            elem_D_points[work_tetra[1]],
            elem_D_points[work_tetra[2]],
            elem_D_points[work_tetra[3]]);

        bElemIntersect = CGAL::do_intersect(m_test_triangle,m_test_tetra);

        if(bElemIntersect)
        {
          // Found intersection!
          break;
        }
      }

      if(bElemIntersect)
      {
        // Found intersection!
        break;
      }
    }
  }
  return bElemIntersect;
}

const libMesh::Elem * Intersection_Tools::FindFirstIntersection(
    const libMesh::Elem * Query_elem,
    std::unique_ptr<libMesh::PointLocatorBase> & point_locator,
    bool        bGuaranteeQueryIsInMesh)
{
  libMesh::PointLocatorBase& locator = *point_locator.get();
  if(!bGuaranteeQueryIsInMesh)
  {
    // Then we are sure that the query element is inside the mesh, only
    // one search needed
  }
  else
  {
    // Better check all the vertices ...
    unsigned int elem_nb_nodes = Query_elem->n_nodes();
    libMesh::Point dummyPoint;
    bool bInsideTheMesh = true;

    // Just to be sure, check if one of the points intersect the mesh
    for(unsigned int iii = 0; iii < elem_nb_nodes; ++iii)
    {
      dummyPoint = Query_elem->point(iii);
      const libMesh::Elem * Patch_elem = locator(Query_elem->point(iii));

      if(Patch_elem == NULL)
      {
        bInsideTheMesh = false;
        break;
      }
    }

    homemade_assert_msg(bInsideTheMesh, "Query element is not fully inside tested mesh!\n");
  }

  return locator(Query_elem->point(0));
};

bool Intersection_Tools::FindAllIntersection(
    const libMesh::Elem * Query_elem,
    std::unique_ptr<libMesh::PointLocatorBase> & point_locator,
    std::set<unsigned int>  & Intersecting_elems)
{
  libMesh::PointLocatorBase& locator = *point_locator.get();

  // Clear output
  Intersecting_elems.clear();

  // Search each vertex
  unsigned int elem_nb_nodes = Query_elem->n_nodes();
  libMesh::Point dummyPoint;

  int nbOfInters = 0;

  // Just to be sure, check if one of the points intersect the mesh
  for(unsigned int iii = 0; iii < elem_nb_nodes; ++iii)
  {
    dummyPoint = Query_elem->point(iii);
    const libMesh::Elem * Patch_elem = locator(Query_elem->point(iii));

    if(Patch_elem != NULL)
    {
      Intersecting_elems.insert(Patch_elem->id());
      ++nbOfInters;
    }
  }

  return nbOfInters != 0;
};

bool Intersection_Tools::libMesh_exact_do_intersect(
    const libMesh::Elem * elem_A,
    const libMesh::Elem * elem_B)
{
  // The booleans
  bool bBboxIntersect = false;
  bool bElemIntersect = false;

  unsigned int n_nodes_A = elem_A->n_nodes();
  unsigned int n_nodes_B = elem_B->n_nodes();

  // First, convert both elements to CGAL exact point vectors
  convert_elem_to_exact_points(elem_A,m_exact_points_A);
  convert_elem_to_exact_points(elem_B,m_exact_points_B);

  // Fast check: bounding boxes
  std::vector<ExactPoint_3>::const_iterator exact_points_A_begin = m_exact_points_A.begin();
  std::vector<ExactPoint_3>::const_iterator exact_points_B_begin = m_exact_points_B.begin();

  Bbox_3 exact_bbox_A = CGAL::bbox_3(exact_points_A_begin,exact_points_A_begin + n_nodes_A);
  Bbox_3 exact_bbox_B = CGAL::bbox_3(exact_points_B_begin,exact_points_B_begin + n_nodes_B);

  bBboxIntersect = CGAL::do_intersect(exact_bbox_A,exact_bbox_B);

  if(bBboxIntersect)
  {
    // Bbox intersect, test intersection between tetrahedrons and triangles

    // Pointers that will depend on the element type;
    std::vector<std::vector<unsigned int> > * elem_A_tetras    = NULL;
    std::vector<std::vector<unsigned int> > * elem_A_triangles = NULL;

    std::vector<std::vector<unsigned int> > * elem_B_tetras    = NULL;
    std::vector<std::vector<unsigned int> > * elem_B_triangles = NULL;

    if(elem_A->type() == libMesh::TET4)
    {
      // Use tetrahedron geometry
      elem_A_tetras = &m_TET_tetrahedrons;
      elem_A_triangles = &m_TET_triangles;
    }
    else if(elem_A->type() == libMesh::HEX8)
    {
      // Use hexaedron geometry
      elem_A_tetras = &m_HEX_tetrahedrons;
      elem_A_triangles = &m_HEX_triangles;
    }
    else
    {
      homemade_error_msg("Unsupported element type! Must be either TET4 or HEX8");
    }

    if(elem_B->type() == libMesh::TET4)
    {
      // Use tetrahedron geometry
      elem_B_tetras = &m_TET_tetrahedrons;
      elem_B_triangles = &m_TET_triangles;
    }
    else if(elem_B->type() == libMesh::HEX8)
    {
      // Use hexaedron geometry
      elem_B_tetras = &m_HEX_tetrahedrons;
      elem_B_triangles = &m_HEX_triangles;
    }
    else
    {
      homemade_error_msg("Unsupported element type! Must be either TET4 or HEX8");
    }

    bElemIntersect = this->elements_do_intersect(m_exact_points_A, *elem_A_tetras, *elem_A_triangles,
        m_exact_points_B, *elem_B_tetras, *elem_B_triangles);
  }
  else
  {
    bElemIntersect = false;
  }

  return bElemIntersect;
}

bool Intersection_Tools::libMesh_exact_intersection(const libMesh::Elem * elem_A,
                const libMesh::Elem * elem_B,
                std::set<Point_3> & points_out,
                bool bCreateNewNefForA,
                bool bConvertPoints,
                bool bTestNeeded)
{
  bool bElemIntersect = true;

  if(bTestNeeded)
  {
    // Test the intersection beforehand
    m_perf_log.push("Test intersection","Exact intersection construction");
    bElemIntersect = libMesh_exact_do_intersect(elem_A,elem_B);
    m_perf_log.pop("Test intersection","Exact intersection construction");
  }
  else if(bConvertPoints)
  {
    // Test already made somewhere else, but we need to set up the exact
    // points.
    m_perf_log.push("Point conversion to exact","Exact intersection construction");
    if(bCreateNewNefForA)
    {
      convert_elem_to_exact_points(elem_A,m_exact_points_A);
    }
    convert_elem_to_exact_points(elem_B,m_exact_points_B);
    m_perf_log.pop("Point conversion to exact","Exact intersection construction");
  }

  if(bElemIntersect)
  {
    // Generate the Nef polyhedrons
    unsigned int n_nodes_A = elem_A->n_nodes();
    unsigned int n_nodes_B = elem_B->n_nodes();

    std::vector<ExactPoint_3>::const_iterator exact_points_A_begin = m_exact_points_A.begin();
    std::vector<ExactPoint_3>::const_iterator exact_points_B_begin = m_exact_points_B.begin();

    m_perf_log.push("Nef polyhedron construction","Exact intersection construction");
    if(bCreateNewNefForA)
    {
      convert_exact_points_to_Nef(  exact_points_A_begin,
                      exact_points_A_begin + n_nodes_A,
                      m_nef_A);
    }
    else
    {
      homemade_assert_msg(!m_nef_A.is_empty(), "Intersection base is empty!\n");
    }

    convert_exact_points_to_Nef(  exact_points_B_begin,
                    exact_points_B_begin + n_nodes_B,
                    m_nef_B);

    m_perf_log.pop("Nef polyhedron construction","Exact intersection construction");

    // Intersect them
    m_perf_log.push("Nef polyhedron intersection","Exact intersection construction");
    m_nef_I = m_nef_A*m_nef_B;
    m_perf_log.pop("Nef polyhedron intersection","Exact intersection construction");

    m_perf_log.push("Point set output","Exact intersection construction");
    if(!m_nef_I.is_empty())
    {
      // Intersection exists! Create output
      bElemIntersect = true;

      points_out.clear();
      for(Nef_Polyhedron::Vertex_const_iterator it_vertex = m_nef_I.vertices_begin();
          it_vertex != m_nef_I.vertices_end();
          ++it_vertex)
      {
        points_out.insert(ConvertExactToInexact(it_vertex->point()));
      }
    }
    else
    {
      bElemIntersect = false;
    }
    m_perf_log.pop("Point set output","Exact intersection construction");
  }

  return bElemIntersect;
}

bool Intersection_Tools::libMesh_exact_do_intersect_inside_coupling(const libMesh::Elem * elem_A,
                        const libMesh::Elem * elem_B,
                        bool bCreateNewNefForA)
{
  bool bElemIntersect = true;

  // Assert if C was built
  homemade_assert_msg(!m_nef_C.is_empty(), "Coupling restriction element was not set yet!\n");

  // Test the intersection beforehand. This also generates the exact
  // points needed.
  bElemIntersect = libMesh_exact_do_intersect(elem_A,elem_B);

  if(bElemIntersect)
  {
    // Generate the Nef polyhedrons
    unsigned int n_nodes_A = elem_A->n_nodes();
    unsigned int n_nodes_B = elem_B->n_nodes();

    std::vector<ExactPoint_3>::const_iterator exact_points_A_begin = m_exact_points_A.begin();
    std::vector<ExactPoint_3>::const_iterator exact_points_B_begin = m_exact_points_B.begin();

    if(bCreateNewNefForA)
    {
      // A was updated, must also update A*C
      convert_exact_points_to_Nef(  exact_points_A_begin,
                      exact_points_A_begin + n_nodes_A,
                      m_nef_A);
      m_nef_I_AC = m_nef_A*m_nef_C;
    }
    else
    {
      homemade_assert_msg(!m_nef_I_AC.is_empty(), "Intersection base is empty!\n");
    }
    convert_exact_points_to_Nef(  exact_points_B_begin,
                    exact_points_B_begin + n_nodes_B,
                    m_nef_B);

    // Intersect them
    m_nef_I = m_nef_I_AC*m_nef_B;

    if(!m_nef_I.is_empty() && m_nef_I.number_of_volumes() > 1)
    {
      // Intersection exists!
      bElemIntersect = true;
    }
    else
    {
      bElemIntersect = false;
    }
  }

  return bElemIntersect;
}

bool Intersection_Tools::libMesh_exact_intersection_inside_coupling(const libMesh::Elem * elem_A,
                        const libMesh::Elem * elem_B,
                        std::set<libMesh::Point> & points_out,
                        bool bCreateNewNefForA,
                        bool bConvertPoints,
                        bool bTestNeeded)
{
  bool bElemIntersect = true;

  // Assert if C was built
  homemade_assert_msg(!m_nef_C.is_empty(), "Coupling restriction element was not set yet!\n");

  // Dummy CGAL point
  Point_3   dummy_CGAL_point;
  libMesh::Point dummy_libMesh_point;

  std::vector<double> bbox_dims(6,0);
  double inter_vol = 0;

  if(bTestNeeded)
  {
    // Test the intersection beforehand
    m_perf_log.push("Test intersection","Exact intersection construction inside coupling");
    bElemIntersect = libMesh_exact_do_intersect(elem_A,elem_B);
    m_perf_log.pop("Test intersection","Exact intersection construction inside coupling");
  }
  else if(bConvertPoints)
  {
    // Test already made somewhere else, but we need to set up the exact
    // points.
    m_perf_log.push("Point conversion to exact","Exact intersection construction inside coupling");
    convert_elem_to_exact_points(elem_A,m_exact_points_A);
    convert_elem_to_exact_points(elem_B,m_exact_points_B);
    m_perf_log.pop("Point conversion to exact","Exact intersection construction inside coupling");
  }

  if(bElemIntersect)
  {
    // Generate the Nef polyhedrons
    unsigned int n_nodes_A = elem_A->n_nodes();
    unsigned int n_nodes_B = elem_B->n_nodes();

    std::vector<ExactPoint_3>::const_iterator exact_points_A_begin = m_exact_points_A.begin();
    std::vector<ExactPoint_3>::const_iterator exact_points_B_begin = m_exact_points_B.begin();

    m_perf_log.push("Nef polyhedron construction","Exact intersection construction inside coupling");
    if(bCreateNewNefForA)
    {
      // Then we need to build a new nef polyhedron
      convert_exact_points_to_Nef(  exact_points_A_begin,
                      exact_points_A_begin + n_nodes_A,
                      m_nef_A);
      m_nef_I_AC = m_nef_A*m_nef_C;
    }
    else
    {
      homemade_assert_msg(!m_nef_I_AC.is_empty(), "Intersection base is empty!\n");
    }

    convert_exact_points_to_Nef(  exact_points_B_begin,
                    exact_points_B_begin + n_nodes_B,
                    m_nef_B);
    m_perf_log.pop("Nef polyhedron construction","Exact intersection construction inside coupling");

    // Intersect them
    m_perf_log.push("Nef polyhedron intersection","Exact intersection construction inside coupling");

    m_nef_I = m_nef_B*m_nef_I_AC;

    m_perf_log.pop("Nef polyhedron intersection","Exact intersection construction inside coupling");

    m_perf_log.push("Point set output","Exact intersection construction inside coupling");

    if(   !m_nef_I.is_empty() &&
         m_nef_I.number_of_volumes() > 1 &&
         m_nef_I.number_of_vertices() > 3 &&
         m_nef_I.number_of_facets() > 1)
    {
      // Intersection exists! Extract points!
      points_out.clear();
      bElemIntersect = true;

      for(Nef_Polyhedron::Vertex_const_iterator it_vertex = m_nef_I.vertices_begin();
          it_vertex != m_nef_I.vertices_end();
          ++it_vertex)
      {
        dummy_CGAL_point = ConvertExactToInexact(it_vertex->point());
        points_out.insert(libMesh::Point(dummy_CGAL_point.x(),dummy_CGAL_point.y(),dummy_CGAL_point.z()));
      }

      // Sometimes, CGAL will generate a very small volume, which follow
      // the "if" conditions above, but results in a point or facet when
      // converted to inexact values. The checks below test for it.
      std::set<libMesh::Point>::const_iterator it_set = points_out.begin();
      dummy_libMesh_point = *it_set;

      bbox_dims[0] = dummy_libMesh_point(0);
      bbox_dims[1] = dummy_libMesh_point(0);
      bbox_dims[2] = dummy_libMesh_point(1);
      bbox_dims[3] = dummy_libMesh_point(1);
      bbox_dims[4] = dummy_libMesh_point(2);
      bbox_dims[5] = dummy_libMesh_point(2);
      ++it_set;

      for(; it_set != points_out.end(); ++it_set)
      {
        dummy_libMesh_point = *it_set;
        if(bbox_dims[0] > dummy_libMesh_point(0)) // Min x
          bbox_dims[0] = dummy_libMesh_point(0);
        if(bbox_dims[1] < dummy_libMesh_point(0)) // Max x
          bbox_dims[1] = dummy_libMesh_point(0);
        if(bbox_dims[2] > dummy_libMesh_point(1)) // Min y
          bbox_dims[2] = dummy_libMesh_point(1);
        if(bbox_dims[3] < dummy_libMesh_point(1)) // Max y
          bbox_dims[3] = dummy_libMesh_point(1);
        if(bbox_dims[4] > dummy_libMesh_point(2)) // Min z
          bbox_dims[4] = dummy_libMesh_point(2);
        if(bbox_dims[5] < dummy_libMesh_point(2)) // Max z
          bbox_dims[5] = dummy_libMesh_point(2);

//        it_set->print();
//        std::cout << std::endl;
      }

      inter_vol = std::abs((bbox_dims[1] - bbox_dims[0]) * (bbox_dims[3] - bbox_dims[2]) * (bbox_dims[5] - bbox_dims[4]));
//      std::cout << " -> Volume = " << inter_vol << std::endl;
      if(points_out.size() < 4 || inter_vol < m_Min_Inter_Volume )
      {
        // Invalid intersection!
        bElemIntersect = false;

        std::set<libMesh::Point>::const_iterator it_set = points_out.begin();
//        std::cout << " -> Volume = " << inter_vol << std::endl;
        for(; it_set != points_out.end(); ++it_set)
        {
//          it_set->print();
//          std::cout << std::endl;
        }

      }
    }
    else
    {
      bElemIntersect = false;
    }

    m_perf_log.pop("Point set output","Exact intersection construction inside coupling");
  }

  return bElemIntersect;
}

void Intersection_Tools::libmesh_set_coupling_nef_polyhedron(const libMesh::Elem * elem_C)
{
  unsigned int n_nodes_C = elem_C->n_nodes();

  // First, convert the element to a CGAL exact point vector
  convert_elem_to_exact_points(elem_C,m_exact_points_C);

  // Second, convert to Nef
  std::vector<ExactPoint_3>::const_iterator exact_points_C_begin = m_exact_points_C.begin();
  convert_exact_points_to_Nef(  exact_points_C_begin,
                  exact_points_C_begin + n_nodes_C,
                  m_nef_C);

  // And, finally, associate the tetra and triangle tables to it
  if(elem_C->type() == libMesh::TET4)
  {
    // Use tetrahedron geometry
    m_elem_C_tetrahedrons = &m_TET_tetrahedrons;
    m_elem_C_triangles = &m_TET_triangles;
  }
  else if(elem_C->type() == libMesh::HEX8)
  {
    // Use hexaedron geometry
    m_elem_C_tetrahedrons = &m_HEX_tetrahedrons;
    m_elem_C_triangles = &m_HEX_triangles;
  }
  else
  {
    homemade_error_msg("Unsupported element type! Must be either TET4 or HEX8");
  }
}

void Intersection_Tools::convert_elem_to_exact_points(  const libMesh::Elem       * elem_input,
                  std::vector<ExactPoint_3> & points_output)
{
  libMesh::Point dummyPoint;
  for(unsigned int iii = 0; iii < elem_input->n_nodes(); ++iii)
  {
    dummyPoint = elem_input->point(iii);
    points_output[iii] = ExactPoint_3(  dummyPoint(0),
                      dummyPoint(1),
                      dummyPoint(2));
  }
}

void Intersection_Tools::convert_exact_points_to_Nef(std::vector<ExactPoint_3>::const_iterator it_begin,
                 std::vector<ExactPoint_3>::const_iterator it_end,
                 Nef_Polyhedron & nef_out)
{
  CGAL::convex_hull_3(it_begin,it_end,m_dummyPoly,ExactHullTraits);
  nef_out = Nef_Polyhedron(m_dummyPoly);
}

}
