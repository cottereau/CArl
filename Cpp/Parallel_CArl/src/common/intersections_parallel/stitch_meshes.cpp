/*
 * stitch_meshes.cpp
 *
 *  Created on: Jun 2, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#include "stitch_meshes.h"

namespace carl
{
void Stitch_Intersection_Meshes::set_base_filenames(const std::string & filename_base, const std::string & mesh_format, int nb_of_files )
{
	// Set the number of files
	if(nb_of_files < 0)
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

	// Create the filenames!
	std::ifstream table_file;
	m_nb_of_intersections = 0;
	m_nb_of_elements = 0;
	m_nb_of_vertices = 0;

	unsigned int temp_nb_of_intersections = 0;
	unsigned int temp_nb_of_elements = 0;
	unsigned int temp_nb_of_vertices = 0;
	for(int iii = 0; iii < m_nb_files; ++iii)
	{
		m_mesh_filenames[iii] = filename_base + std::to_string(iii) + "_n_" + std::to_string(m_nb_files) + mesh_format;
		m_table_filenames[iii] = filename_base + std::to_string(iii) + "_n_" + std::to_string(m_nb_files) + "_inter_table_Full.dat";

		// Read the first line of each intersection table file, to get the preallocations
		if(m_rank == 0)
		{
			table_file.open(m_table_filenames[iii]);
			table_file >> temp_nb_of_intersections;
			table_file >> temp_nb_of_elements;
			table_file >> temp_nb_of_vertices;
			table_file.close();

			m_nb_of_intersections += temp_nb_of_intersections;
			m_nb_of_elements += temp_nb_of_elements;
			m_nb_of_vertices += temp_nb_of_vertices;
		}
	}

	// Filenames have been set
	m_bFilenamesSet = true;

	// Broadcast the data
	m_comm.broadcast(m_nb_of_intersections);
	m_comm.broadcast(m_nb_of_elements);
	m_comm.broadcast(m_nb_of_vertices);

	// Preallocate the data - for now, only in proc 0!
	if(m_rank == 0)
	{
		preallocate_grid(2*m_nb_of_vertices);
	}
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

}


