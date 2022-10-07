/*
 * stitch_meshes.cpp
 *
 *  Created on: Jun 2, 2016
 *      Author: Thiago Milanetto Schlittler
 *
 *  \brief **STAT/DYN-DI/DYN-CG** Intersection step
 */

#include "stitch_meshes.h"

namespace carl
{
const libMesh::ReplicatedMesh & Stitch_Meshes::mesh()
{
  return m_Stitched_mesh;
}

void Stitch_Meshes::set_base_filenames(std::vector<std::string> & mesh_filenames, std::vector<std::string> & table_filenames)
{
  homemade_assert_msg(mesh_filenames.size() == table_filenames.size(), "File lists have different sizes!\n");

  m_nb_files = table_filenames.size();

  // Resize the string vectors
  m_mesh_filenames.resize(m_nb_files);
  m_table_filenames.resize(m_nb_files);

  for(unsigned int iii = 0; iii < m_nb_files; ++iii)
  {
    m_mesh_filenames[iii] = mesh_filenames[iii];
    m_table_filenames[iii] = table_filenames[iii];
  }

  m_bFilenamesSet = true;
}

void Stitch_Meshes::set_base_filenames(const std::string & filename_base, const std::string & mesh_format, unsigned int nb_of_files )
{
  // Set the number of files
  if(nb_of_files == 0)
  {
    m_nb_files = m_nodes;
  }
  else
  {
    m_nb_files = nb_of_files;
  }

  // Resize the string vectors
  m_mesh_filenames.resize(m_nb_files);
  m_table_filenames.resize(m_nb_files);

  for(unsigned int iii = 0; iii < m_nb_files; ++iii)
  {
    m_mesh_filenames[iii] = filename_base + "_r_" + std::to_string(iii) + "_n_" + std::to_string(m_nb_files) + mesh_format;
    m_table_filenames[iii] = filename_base + "_r_" + std::to_string(iii) + "_n_" + std::to_string(m_nb_files) + "_inter_table.dat";
  }

  m_mesh_output = m_base_output + ".e";
  m_table_output= m_base_output + "_inter_pairs.dat";

  m_bFilenamesSet = true;
}

void Stitch_Meshes::preallocate_grid(int map_preallocation)
{
  m_discrete_vertices.reserve(map_preallocation);
//  m_Grid_to_mesh_vertex_idx.reserve(map_preallocation);
  m_bGridPreallocated = true;
}

void Stitch_Meshes::set_grid_constraints(const libMesh::Mesh & mesh_A, const libMesh::Mesh & mesh_B, double vol_tol )
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

  for(unsigned int iii = 0; iii < 3; ++iii)
  {
    m_GridN[iii] = (m_Grid_MaxPoint(iii) - m_Grid_MinPoint(iii)) / m_eps + 1;
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

  // Mark grid as ready to use
  m_bGridDefined = true;

  if(m_bPrintDebug)
  {
    std::cout << "    DEBUG: discrete grid" << std::endl;
    std::cout << " -> eps             : " << m_eps << std::endl;
    std::cout << " -> volume          : " << m_vol_tol << std::endl;
    std::cout << " -> Grid dimensions : " << m_GridN[0] << " " << m_GridN[1] << " " << m_GridN[2] << " " << std::endl  << std::endl;
  }
}

void Stitch_Meshes::set_grid_constraints(Mesh_Intersection & mesh_inter_obj)
{
  // Copy everything from the argument
  m_Grid_MinPoint = mesh_inter_obj.min_point();
  m_Grid_MaxPoint = mesh_inter_obj.max_point();
  m_eps = mesh_inter_obj.eps();
  m_vol_tol = mesh_inter_obj.min_vol();
  m_GridN = mesh_inter_obj.grid_sizes();
  m_GridN_min = mesh_inter_obj.grid_min_size();

  // Mark grid as ready to use
  m_bGridDefined = true;

  if(m_bPrintDebug)
  {
    std::cout << "    DEBUG: discrete grid" << std::endl;
    std::cout << " -> eps             : " << m_eps << std::endl;
    std::cout << " -> volume          : " << m_vol_tol << std::endl;
    std::cout << " -> Grid dimensions : " << m_GridN[0] << " " << m_GridN[1] << " " << m_GridN[2] << " " << std::endl  << std::endl;
  }
}

void Stitch_Meshes::join_tables()
{
  // Check if the grids and files were set, and if the grid was preallocated
  homemade_assert_msg(m_bGridDefined,"Grid not set!\n");
  homemade_assert_msg(m_bFilenamesSet,"Filenames not set!\n");

  // -> First, read the intersection tables files, to get the grid hash
  //    function preallocation.

  std::ifstream table_file;
  m_nb_of_intersections = 0;
  m_nb_of_elements = 0;
  m_maximum_nb_of_nodes = 0;

  unsigned int temp_nb_of_intersections = 0;
  unsigned int temp_nb_of_elements = 0;
  unsigned int temp_nb_of_nodes = 0;

  for(unsigned int iii = 0; iii < m_nb_files; ++iii)
  {
    table_file.open(m_table_filenames[iii]);
    table_file >> temp_nb_of_intersections;
    table_file >> temp_nb_of_elements;
    table_file >> temp_nb_of_nodes;

    table_file.close();

    m_nb_of_intersections += temp_nb_of_intersections;
    m_nb_of_elements += temp_nb_of_elements;
    m_maximum_nb_of_nodes += temp_nb_of_nodes;
  }

  // -> Preallocate the data structures
  preallocate_grid(100*m_maximum_nb_of_nodes);

  m_Stitched_mesh.reserve_elem(m_nb_of_elements);
  m_Stitched_mesh.reserve_nodes(m_maximum_nb_of_nodes);

  m_intersection_pairs.resize(m_nb_of_intersections);
  m_intersection_nb_of_elements.resize(m_nb_of_intersections);
  m_restriction_set_first.reserve(m_nb_of_intersections);
  m_restriction_set_second.reserve(m_nb_of_intersections);

  // -> Second, read the data that will be used to reconstruct the
  //    intersection tables in the end: intersection pairs, and number of
  //    elements per intersection.
  unsigned int intersection_idx = 0;
  unsigned int dummy_uint = 0;
  unsigned int nb_of_elems = 0;
  for(unsigned int iii = 0; iii < m_nb_files; ++iii)
  {
    table_file.open(m_table_filenames[iii]);
    table_file >> temp_nb_of_intersections;
    table_file >> temp_nb_of_elements;
    table_file >> temp_nb_of_nodes;

    for(unsigned int jjj = 0; jjj < temp_nb_of_intersections; ++jjj)
    {
      table_file  >> dummy_uint;
      table_file  >> m_intersection_pairs[intersection_idx].first
            >> m_intersection_pairs[intersection_idx].second;
      table_file  >> m_intersection_nb_of_elements[intersection_idx];
      carl::jump_lines(table_file);

      // Prepare data for the restriction mesh
      nb_of_elems += m_intersection_nb_of_elements[intersection_idx];
      m_restriction_set_first.insert(m_intersection_pairs[intersection_idx].first);
      m_restriction_set_second.insert(m_intersection_pairs[intersection_idx].second);
      ++intersection_idx;
    }

    table_file.close();
  }

  // -> Fourth, re-build the intersection tables
  unsigned int intersection_elem_idx = 0;

  std::ofstream joined_tables_file(m_table_output);
  joined_tables_file  << m_nb_of_intersections << " "
            << nb_of_elems << std::endl;

  for(unsigned int iii = 0; iii < m_nb_of_intersections; ++iii)
  {
    joined_tables_file  << m_intersection_pairs[iii].first << " "
              << m_intersection_pairs[iii].second << std::endl;
    // joined_tables_file   << iii << " "
    //          << m_intersection_pairs[iii].first << " "
    //          << m_intersection_pairs[iii].second << " "
    //          << m_intersection_nb_of_elements[iii] << " ";
    // for(unsigned jjj = 0; jjj < m_intersection_nb_of_elements[iii]; ++jjj)
    // {
    //  joined_tables_file  << intersection_elem_idx << " ";
    //  ++intersection_elem_idx;
    // }
    // joined_tables_file   << std::endl;
  }
  joined_tables_file.close();
}

void Stitch_Meshes::stitch_meshes()
{
  homemade_assert_msg(m_bGridPreallocated,"Grid not preallocated!\n");

  // -> Third, stitch the meshes
  libMesh::ReplicatedMesh temp_mesh(m_world_comm,3);
  temp_mesh.allow_renumbering(false);

  unsigned int full_mesh_nb_elems = 0;
  unsigned int full_mesh_nb_nodes = 0;

  libMesh::Elem * copy_elem = NULL;
  libMesh::Elem * mesh_elem = NULL;
  libMesh::Node * mesh_node = NULL;

  double dummy_volume = 0;

  // -> Stitch!
  for(unsigned int iii = 0; iii < m_nb_files; ++iii)
  {
    // -> Open mesh file
    temp_mesh.read(m_mesh_filenames[iii]);

    // -> Insert nodes
    libMesh::ReplicatedMesh::element_iterator it_mesh = temp_mesh.elements_begin();

    for( ; it_mesh != temp_mesh.elements_end(); ++it_mesh)
    {
      copy_elem = * it_mesh;

      // -> Each element is unique, so no tests for insertion
      mesh_elem = libMesh::Elem::build(libMesh::TET4).release();
      mesh_elem->set_id(full_mesh_nb_elems);
      mesh_elem->processor_id(0);

      // -> First, add the nodes
      for(unsigned int jjj = 0; jjj < 4; ++jjj)
      {
        convert_to_discrete(copy_elem->point(jjj),m_dummy_discrete_point);
        if(m_discrete_vertices.find(m_dummy_discrete_point) == m_discrete_vertices.end())
        {
          // New vertex! Add it to the mesh
          m_discrete_vertices[m_dummy_discrete_point] = full_mesh_nb_nodes;
          mesh_node = m_Stitched_mesh.add_point(copy_elem->point(jjj),full_mesh_nb_nodes,0);
          ++full_mesh_nb_nodes;
        }
        else
        {
          mesh_node = m_Stitched_mesh.node_ptr(m_discrete_vertices[m_dummy_discrete_point]);
        }

        // Associate vertex to the new element
        mesh_elem->set_node(jjj) = mesh_node;

      }

      m_Stitched_mesh.add_elem(mesh_elem);
      dummy_volume += mesh_elem->volume();

      ++full_mesh_nb_elems;
    }
  }

  // Print information about the number of collisions
  if(m_bPrintDebug)
  {
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
  m_Stitched_mesh.prepare_for_use();

  // Print mesh
  libMesh::NameBasedIO output_mesh(m_Stitched_mesh);
  output_mesh.write(m_mesh_output);

  if(m_bPrintDebug)
  {
    std::cout << "    DEBUG: stitched mesh" << std::endl;
    std::cout << " -> Volume : " << dummy_volume << std::endl  << std::endl;
  }

  int wrong_volume = 0;
  libMesh::ReplicatedMesh::element_iterator elem_begin = m_Stitched_mesh.local_elements_begin();
  libMesh::ReplicatedMesh::element_iterator elem_end = m_Stitched_mesh.local_elements_end();

  for( ; elem_begin != elem_end; ++elem_begin)
  {
    libMesh::Elem * dummy_elem = * elem_begin;
    if(std::abs(dummy_elem->volume()) < m_vol_tol)
    {
      ++wrong_volume;
    }
  }

  std::cout << " -> bad volumes : " << wrong_volume << " ( " <<  m_vol_tol << " ) " << std::endl;
};

void Stitch_Meshes::convert_to_discrete(const libMesh::Point& iPoint, std::vector<long>& oPoint)
{
  oPoint[0] = lround( (iPoint(0) -  m_Grid_MinPoint(0) )/m_eps);
  oPoint[1] = lround( (iPoint(1) -  m_Grid_MinPoint(1) )/m_eps);
  oPoint[2] = lround( (iPoint(2) -  m_Grid_MinPoint(2) )/m_eps);
}

const std::unordered_set<unsigned int>* Stitch_Meshes::get_restricted_set_pointer_first()
{
  return &m_restriction_set_first;
};

const std::unordered_set<unsigned int>* Stitch_Meshes::get_restricted_set_pointer_second()
{
  return &m_restriction_set_second;
};

}




