/*
 * stitch_meshes.cpp
 *
 *  Created on: Jun 2, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#include "stitch_meshes.h"

namespace carl
{
const libMesh::SerialMesh & Stitch_Intersection_Meshes::mesh()
{
	return m_Stitched_mesh;
}

void Stitch_Intersection_Meshes::set_base_filenames(std::vector<std::string> & mesh_filenames, std::vector<std::string> & table_filenames)
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

void Stitch_Intersection_Meshes::set_base_filenames(const std::string & filename_base, const std::string & mesh_format, unsigned int nb_of_files )
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
		m_table_filenames[iii] = filename_base + "_r_" + std::to_string(iii) + "_n_" + std::to_string(m_nb_files) + "_inter_table_Full.dat";
	}

	m_mesh_output = m_base_output + ".msh";
	m_table_output= m_base_output + "_inter_table_Full.dat";

	m_bFilenamesSet = true;
}

void Stitch_Intersection_Meshes::preallocate_grid(int map_preallocation)
{
	m_Grid_to_mesh_vertex_idx.reserve(map_preallocation);
	m_bGridPreallocated = true;
}

void Stitch_Intersection_Meshes::set_grid_constraints(const libMesh::Mesh & mesh_A, const libMesh::Mesh & mesh_B)
{
	libMesh::MeshTools::BoundingBox bbox_A = libMesh::MeshTools::bounding_box(mesh_A);
	libMesh::MeshTools::BoundingBox bbox_B = libMesh::MeshTools::bounding_box(mesh_B);

	// Just to be sure, test if the bboxes intersect!
	homemade_assert_msg(bbox_A.intersect(bbox_B),"Meshes' bounding boxes do not intersect!\n");

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

	// Mark grid as ready to use
	m_bGridDefined = true;

	if(m_bPrintDebug)
	{
		std::cout << "    DEBUG: discrete grid" << std::endl;
		std::cout << " -> eps             : " << m_eps << std::endl;
		std::cout << " -> Grid dimensions : " << m_GridN[0] << " " << m_GridN[1] << " " << m_GridN[2] << " " << std::endl  << std::endl;
	}
}

void Stitch_Intersection_Meshes::set_grid_constraints(Mesh_Intersection & mesh_inter_obj)
{
	// Copy everything from the argument
	m_Grid_MinPoint = mesh_inter_obj.min_point();
	m_Grid_MaxPoint = mesh_inter_obj.max_point();
	m_eps = mesh_inter_obj.eps();
	m_GridN = mesh_inter_obj.grid_sizes();
	m_GridN_min = mesh_inter_obj.grid_min_size();

	// Mark grid as ready to use
	m_bGridDefined = true;

	if(m_bPrintDebug)
	{
		std::cout << "    DEBUG: discrete grid" << std::endl;
		std::cout << " -> eps             : " << m_eps << std::endl;
		std::cout << " -> Grid dimensions : " << m_GridN[0] << " " << m_GridN[1] << " " << m_GridN[2] << " " << std::endl  << std::endl;
	}
}

void Stitch_Intersection_Meshes::stitch_meshes()
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
	preallocate_grid(2*m_maximum_nb_of_nodes);

	m_Stitched_mesh.reserve_elem(m_nb_of_elements);
	m_Stitched_mesh.reserve_nodes(m_maximum_nb_of_nodes);

	m_intersection_pairs.resize(m_nb_of_intersections);
	m_intersection_nb_of_elements.resize(m_nb_of_intersections);

	// -> Second, read the data that will be used to reconstruct the
	//    intersection tables in the end: intersection pairs, and number of
	//    elements per intersection.
	unsigned int intersection_idx = 0;
	unsigned int dummy_uint = 0;
	for(unsigned int iii = 0; iii < m_nb_files; ++iii)
	{
		table_file.open(m_table_filenames[iii]);
		table_file >> temp_nb_of_intersections;
		table_file >> temp_nb_of_elements;
		table_file >> temp_nb_of_nodes;

		for(unsigned int jjj = 0; jjj < temp_nb_of_intersections; ++jjj)
		{
			table_file 	>> dummy_uint;
			table_file 	>> m_intersection_pairs[intersection_idx].first
						>> m_intersection_pairs[intersection_idx].second;
			table_file 	>> m_intersection_nb_of_elements[intersection_idx];
			carl::jump_lines(table_file);
			++intersection_idx;
		}

		table_file.close();
	}

	homemade_assert_msg(m_bGridPreallocated,"Grid not preallocated!\n");

	// -> Third, stitch the meshes
	libMesh::SerialMesh temp_mesh(m_world_comm,3);
	temp_mesh.allow_renumbering(false);
	long grid_value = 0;

	unsigned int full_mesh_nb_elems = 0;
	unsigned int full_mesh_nb_nodes = 0;
	unsigned int partial_mesh_nb_elems = 0;
	unsigned int partial_mesh_nb_nodes = 0;

	libMesh::Elem * copy_elem = NULL;
	libMesh::Elem * mesh_elem = NULL;
	libMesh::Node * mesh_node = NULL;

	double dummy_volume = 0;

	// -> Stitch!
	for(unsigned int iii = 0; iii < m_nb_files; ++iii)
	{
		// -> Open mesh file
		temp_mesh.read(m_mesh_filenames[iii]);

		partial_mesh_nb_elems = temp_mesh.n_elem();
		partial_mesh_nb_nodes = temp_mesh.n_nodes();

		// -> Insert nodes
		libMesh::SerialMesh::element_iterator it_mesh = temp_mesh.elements_begin();
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
				grid_value = convert_to_grid(copy_elem->point(jjj));
				if(m_Grid_to_mesh_vertex_idx.find(grid_value)==m_Grid_to_mesh_vertex_idx.end())
				{
					// New vertex! Add it to the mesh
					m_Grid_to_mesh_vertex_idx[grid_value] = full_mesh_nb_nodes;
					mesh_node = m_Stitched_mesh.add_point(copy_elem->point(jjj),full_mesh_nb_nodes,0);
					++full_mesh_nb_nodes;
				}
				else
				{
					mesh_node = m_Stitched_mesh.node_ptr(m_Grid_to_mesh_vertex_idx[grid_value]);
				}

				// Associate vertex to the new element
				mesh_elem->set_node(jjj) = mesh_node;

			}

			m_Stitched_mesh.add_elem(mesh_elem);
			dummy_volume += mesh_elem->volume();

			++full_mesh_nb_elems;
		}
	}

	m_Stitched_mesh.prepare_for_use();
	m_Stitched_mesh.write(m_mesh_output);

	// -> Fourth, re-build the intersection tables
	unsigned int intersection_elem_idx = 0;

	std::ofstream joined_tables_file(m_table_output);
	joined_tables_file 	<< m_nb_of_intersections << " "
						<< m_Stitched_mesh.n_elem() << std::endl;
	for(unsigned int iii = 0; iii < m_nb_of_intersections; ++iii)
	{
		joined_tables_file 	<< iii << " "
							<< m_intersection_pairs[iii].first << " "
							<< m_intersection_pairs[iii].second << " "
							<< m_intersection_nb_of_elements[iii] << " ";
		for(unsigned jjj = 0; jjj < m_intersection_nb_of_elements[iii]; ++jjj)
		{
			joined_tables_file 	<< intersection_elem_idx + 1 << " ";
			++intersection_elem_idx;
		}
		joined_tables_file 	<< std::endl;
	}
	joined_tables_file.close();

	if(m_bPrintDebug)
	{
		std::cout << "    DEBUG: stitched mesh" << std::endl;
		std::cout << " -> Volume : " << dummy_volume << std::endl  << std::endl;
	}
};

long Stitch_Intersection_Meshes::convert_to_grid(const libMesh::Point iPoint)
{
	long dummy =  lround( (iPoint(0) -  m_Grid_MinPoint(0) )/m_eps) * m_GridN[1]*m_GridN[2]
				+ lround( (iPoint(1) -  m_Grid_MinPoint(1) )/m_eps) * m_GridN[1]
				+ lround( (iPoint(2) -  m_Grid_MinPoint(2) )/m_eps);
	homemade_assert_msg(dummy > -1, "Negative grid index!\n");

	return dummy;
}
}


