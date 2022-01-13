#include "mesh_intersection_methods.h"

namespace carl
{

void Mesh_Intersection::triangulate_intersection(const std::set<libMesh::Point> & input_points)
{
  m_libMesh_PolyhedronMesh.clear();

  for(std::set<libMesh::Point>::const_iterator it_set = input_points.begin();
      it_set != input_points.end();
      ++it_set)
  {
    libMesh::Point dummy = *it_set;
    m_libMesh_PolyhedronMesh.add_point(*it_set);
  }

  switch (m_MeshingMethod)
  {
  case IntersectionMeshingMethod::LIBMESH_TETGEN:
  {
    libMesh::TetGenMeshInterface      temp_tetgen_interface(m_libMesh_PolyhedronMesh);
    temp_tetgen_interface.triangulate_pointset();
    break;
  }
  case IntersectionMeshingMethod::CGAL:
  {
    // Clear CGAL's mesh
    m_CGAL_PolyhedronMesh.clear();

    // This is a superfluous extra point conversion:
    // (CGAL Exact -> CGAL Inexact -> LibMesh -> CGAL Inexact)
    // But changing this would need a re-structuring of the code (TODO !)
    Point_3 dummy_CGAL_point;
    const libMesh::Point * dummy_libMesh_point;
    for(int iii = 0; iii < input_points.size(); ++iii)
    {
      dummy_libMesh_point = &m_libMesh_PolyhedronMesh.point(iii);
      dummy_CGAL_point = Point_3( (*dummy_libMesh_point)(0),(*dummy_libMesh_point)(1),(*dummy_libMesh_point)(2));

      // Insert vertex and save the corresponding
      Vertex_handle_3 dummy_vertex_handle = m_CGAL_PolyhedronMesh.insert(dummy_CGAL_point);
      dummy_vertex_handle->info().ExtIndex = iii;
    }

    // At this point, CGAL generated the mesh, now we have to convert it to libMesh
    unsigned int dummy_vertex_idx = 0;
    for(  Finite_cells_iterator_3 cell_it = m_CGAL_PolyhedronMesh.finite_cells_begin();
        cell_it != m_CGAL_PolyhedronMesh.finite_cells_end();
        ++cell_it)
    {
      libMesh::Elem * elem_pointer =  m_libMesh_PolyhedronMesh.add_elem(new libMesh::Tet4);
      for(int iii = 0; iii < 4; ++iii)
      {
        dummy_vertex_idx = cell_it->vertex(iii)->info().ExtIndex;
        elem_pointer->set_node(iii) = &(m_libMesh_PolyhedronMesh.node_ref(dummy_vertex_idx));
      }
    }
    m_libMesh_PolyhedronMesh.prepare_for_use();

    break;
  }
  }
}

void Mesh_Intersection::update_intersection_mesh()
{
  // Test which elements must be added, and insert their vertices
  libMesh::ReplicatedMesh::element_iterator   it_poly_mesh = m_libMesh_PolyhedronMesh.elements_begin();

  for(  ; it_poly_mesh != m_libMesh_PolyhedronMesh.elements_end();
      ++it_poly_mesh)
  {
    libMesh::Elem * poly_elem = * it_poly_mesh;

    if(std::abs(poly_elem->volume()) > m_vol_tol)
    {
      // Add a new element!
      libMesh::Elem * mesh_elem = m_libMesh_Mesh.add_elem(new libMesh::Tet4);
      
      // -> First, add the nodes
      for(unsigned int iii = 0; iii < 4; ++iii)
      {
        convert_to_discrete(poly_elem->point(iii),m_dummy_discrete_point);

        if(m_discrete_vertices.find(m_dummy_discrete_point) == m_discrete_vertices.end())
        {
          // New vertex! Add it to the mesh

          m_discrete_vertices[m_dummy_discrete_point] = m_nb_of_vertices;
          m_libMesh_Mesh.add_point(poly_elem->point(iii),m_nb_of_vertices);
          mesh_elem->set_node(iii) = m_libMesh_Mesh.node_ptr(m_nb_of_vertices);
          ++m_nb_of_vertices;
        }

        // Associate vertex to the new element
        mesh_elem->set_node(iii) = m_libMesh_Mesh.node_ptr(m_discrete_vertices[m_dummy_discrete_point]);
      }
      // Increase the number of elements
      ++m_nb_of_elements;
    }
  }
}

void Mesh_Intersection::update_intersection_pairs(unsigned int elem_idx_A, unsigned int elem_idx_B, unsigned int inter_id)
{
  m_intersection_pairs[inter_id] = std::make_pair(elem_idx_A, elem_idx_B);
}

void Mesh_Intersection::update_intersection_element_range(unsigned int range_start, unsigned int range_end, unsigned int inter_id)
{
  m_intersection_element_range[inter_id] = std::make_pair(range_start, range_end);
}

const libMesh::ReplicatedMesh & Mesh_Intersection::mesh()
{
  return m_libMesh_Mesh;
}

libMesh::Point & Mesh_Intersection::min_point()
{
  return m_Grid_MinPoint;
}

libMesh::Point & Mesh_Intersection::max_point()
{
  return m_Grid_MaxPoint;
}

double Mesh_Intersection::eps()
{
  return m_eps;
}

double Mesh_Intersection::min_vol()
{
  return m_vol_tol;
}

std::vector<long> & Mesh_Intersection::grid_sizes()
{
  return m_GridN;
}

long Mesh_Intersection::grid_min_size()
{
  return m_GridN_min;
}

void Mesh_Intersection::initialize()
{
  m_libMesh_Mesh.clear();
  m_discrete_vertices.clear();
  m_bMeshFinalized = false;
  m_intersection_pairs.clear();
  m_intersection_couplings.clear();
  m_intersection_element_range.clear();
  m_nb_of_intersections = 0;
  m_nb_of_elements = 0 ;
  m_nb_of_vertices = 0 ;
  m_nb_of_points = 0 ;
}

void Mesh_Intersection::preallocate_grid( int map_preallocation )
{
//  m_Grid_to_mesh_vertex_idx.reserve(map_preallocation);
  m_discrete_vertices.reserve(map_preallocation);
}

void Mesh_Intersection::set_grid_constraints(const libMesh::Mesh & mesh_A, const libMesh::Mesh & mesh_B, double vol_tol)
{
  libMesh::BoundingBox bbox_A = libMesh::MeshTools::create_bounding_box(mesh_A);
  libMesh::BoundingBox bbox_B = libMesh::MeshTools::create_bounding_box(mesh_B);

  // Just to be sure, test if the bboxes intersect!
  homemade_assert_msg(bbox_A.intersects(bbox_B),"Meshes' bounding boxes do not intersect!\n");

  // Set the (future) intersection bbox corners
  std::vector<double> eps_candidates(3,0);
  for(unsigned int iii = 0; iii < 3; ++iii)
  {
    m_Grid_MinPoint(iii) = std::max(bbox_A.min()(iii),bbox_B.min()(iii));
    m_Grid_MaxPoint(iii) = std::min(bbox_A.max()(iii),bbox_B.max()(iii));
    eps_candidates[iii] = (m_Grid_MaxPoint(iii) - m_Grid_MinPoint(iii))/m_GridN_min;
  }

  m_eps = *std::min_element(eps_candidates.begin(),eps_candidates.end());

  for(unsigned int iii = 0; iii < 3; ++iii)
  {
    m_Grid_MinPoint(iii) -= 2*m_eps;
    m_Grid_MaxPoint(iii) += 2*m_eps;
  }

  if( vol_tol < 0 )
  {
    // Grossily estimate the volume of A's and B's elements using the bbox
    double grid_volume    =   (m_Grid_MaxPoint(0) - m_Grid_MinPoint(0)) *
                  (m_Grid_MaxPoint(1) - m_Grid_MinPoint(1)) *
                  (m_Grid_MaxPoint(2) - m_Grid_MinPoint(2));


    double fraction_vol_A =   (bbox_A.max()(0) - bbox_A.min()(0)) *
                  (bbox_A.max()(1) - bbox_A.min()(1)) *
                  (bbox_A.max()(2) - bbox_A.min()(2)) /
                  grid_volume;

    double fraction_vol_B =   (bbox_B.max()(0) - bbox_B.min()(0)) *
                  (bbox_B.max()(1) - bbox_B.min()(1)) *
                  (bbox_B.max()(2) - bbox_B.min()(2)) /
                  grid_volume;

    unsigned int est_elem =   std::max(fraction_vol_A * mesh_A.n_elem(),fraction_vol_B * mesh_B.n_elem());

    m_vol_tol = 1E-6 * grid_volume / est_elem;
  }
  else
  {
    m_vol_tol = vol_tol;
  }

  for(unsigned int iii = 0; iii < 3; ++iii)
  {
    m_GridN[iii] = (m_Grid_MaxPoint(iii) - m_Grid_MinPoint(iii)) / m_eps + 1;
  }

  if(m_bPrintDebug)
  {
    std::cout << "    DEBUG: discrete grid" << std::endl;
    std::cout << " -> eps             : " << m_eps << std::endl;
    std::cout << " -> volume          : " << m_vol_tol << std::endl;
    std::cout << " -> Grid dimensions : " << m_GridN[0] << " " << m_GridN[1] << " " << m_GridN[2] << " " << std::endl  << std::endl;
  }
}

void Mesh_Intersection::increase_intersection_mesh( const std::set<libMesh::Point> & input_points,
                  unsigned int elem_idx_A,
                  unsigned int elem_idx_B,
                  unsigned int elem_idx_C)
{
  m_bMeshFinalized = false;
  unsigned int intersection_range_start = m_nb_of_elements;

  //  First, triangulate the point set!
  triangulate_intersection(input_points);

  //  Second, add the points to the grid-to-mesh map and update the intersection mesh,
  // taking into account if the associated elements are big enough.
  update_intersection_mesh();

  //  Third, update the intersection pairs, if needed
  unsigned int intersection_range_end = m_nb_of_elements;
  if(intersection_range_end != intersection_range_start)
  {
    update_intersection_pairs(elem_idx_A, elem_idx_B,m_nb_of_intersections);
    m_intersection_couplings[m_nb_of_intersections] = elem_idx_C;
    update_intersection_element_range(intersection_range_start,intersection_range_end,m_nb_of_intersections);
    ++m_nb_of_intersections;
  }
}

void Mesh_Intersection::export_intersection_data(const std::string & filename_base, const std::string & mesh_format)
{
  homemade_assert_msg(m_bMeshFinalized, "Mesh not prepared for use yet!\n");

  // Set the filenames depending on the processor rank
  std::string mesh_file_out = filename_base + mesh_format;
  std::string table_file_out = filename_base + "_inter_table.dat";

  // Print the mesh
  libMesh::NameBasedIO output_mesh(m_libMesh_Mesh);
  output_mesh.write(mesh_file_out);

  // Print the intersection table
  std::ofstream table_out(table_file_out);
  unsigned int nb_of_inter_elements = 0;
  table_out << m_nb_of_intersections << " " << m_libMesh_Mesh.n_elem() << " " << m_libMesh_Mesh.n_nodes() << std::endl;

  for(unsigned int iii = 0; iii < m_nb_of_intersections; ++iii)
  {
    nb_of_inter_elements = m_intersection_element_range[iii].second
               - m_intersection_element_range[iii].first;
    table_out << iii << " " << m_intersection_pairs[iii].first << " "
                        << m_intersection_pairs[iii].second << " "
                << nb_of_inter_elements;
    for(unsigned int jjj = 0; jjj < nb_of_inter_elements; ++jjj)
    {
      table_out << " " << m_intersection_element_range[iii].first + jjj;
    }
    table_out << std::endl;
  }
  table_out.close();
}

void Mesh_Intersection::prepare_for_use()
{
  if(m_bPrintDebug)
  {
    // Print information about the number of collisions
    size_t collisions = 0;
    for (size_t bucket = 0; bucket != m_discrete_vertices.bucket_count(); ++bucket)
    {
        if (m_discrete_vertices.bucket_size(bucket) > 1)
        {
            collisions += m_discrete_vertices.bucket_size(bucket) - 1;
        }
    }

    std::cout   << "    DEBUG: discrete grid hash collisions" << std::endl;
    std::cout   << " -> Nb. of collisions / size : " << collisions << " / " << m_discrete_vertices.size()
          << " (" << 100.*collisions/m_discrete_vertices.size() << "%)" << std::endl << std::endl;
  }

  m_libMesh_Mesh.prepare_for_use();
  m_bMeshFinalized = true;
}

double Mesh_Intersection::get_total_volume()
{
  double volume = 0;

  homemade_assert_msg(m_bMeshFinalized,"Intersection mesh was not prepared for use!\n");

  for(libMesh::ReplicatedMesh::element_iterator it_elem = m_libMesh_Mesh.elements_begin();
      it_elem != m_libMesh_Mesh.elements_end();
      ++it_elem)
  {
    const libMesh::Elem * elem = * it_elem;
    volume += elem->volume();
  }

  return volume;
}

double Mesh_Intersection::get_intersection_volume(std::set<libMesh::Point> & input_points)
{
  double volume = 0;

  triangulate_intersection(input_points);

  for(libMesh::ReplicatedMesh::element_iterator it_elem = m_libMesh_PolyhedronMesh.elements_begin();
      it_elem != m_libMesh_PolyhedronMesh.elements_end();
      ++it_elem)
  {
    const libMesh::Elem * elem = * it_elem;
    volume += elem->volume();
  }

  return volume;
}

void Mesh_Intersection::convert_to_discrete(const libMesh::Point& iPoint, std::vector<long>& oPoint)
{
  oPoint[0] = lround( (iPoint(0) -  m_Grid_MinPoint(0) )/m_eps);
  oPoint[1] = lround( (iPoint(1) -  m_Grid_MinPoint(1) )/m_eps);
  oPoint[2] = lround( (iPoint(2) -  m_Grid_MinPoint(2) )/m_eps);
}
}
