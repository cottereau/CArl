#include "mesh_intersection_methods.h"

namespace carl
{

void Mesh_Intersection::update_intersection_vertices(	const std::set<libMesh::Point> & input_points)
{
	long dummy_long_int;
	m_nb_of_points = 0;

	for(	std::set<libMesh::Point>::const_iterator it_points = input_points.begin();
			it_points != input_points.end();
			++it_points)
	{
		dummy_long_int = convert_to_grid(*it_points);
		if(m_Grid_to_mesh_vertex_idx.find(dummy_long_int)==m_Grid_to_mesh_vertex_idx.end())
		{
			// New vertex! Add it to the mesh
			m_Grid_to_mesh_vertex_idx[dummy_long_int] = m_nb_of_vertices;
			m_intersection_point_indexes[m_nb_of_points] = m_nb_of_vertices;
			m_libMesh_Mesh.add_point(*it_points,m_nb_of_vertices);
			++m_nb_of_vertices;
			++m_nb_of_points;
		}
		else
		{
			// Recover the index corresponding to the point
			m_intersection_point_indexes[m_nb_of_points] = m_Grid_to_mesh_vertex_idx[dummy_long_int];
			++m_nb_of_points;
		}
	}
}

void Mesh_Intersection::triangulate_intersection(const std::set<libMesh::Point> & input_points)
{
	m_libMesh_PolyhedronMesh.clear();

	for(std::set<libMesh::Point>::const_iterator it_set = input_points.begin();
			it_set != input_points.end();
			++it_set)
	{
		m_libMesh_PolyhedronMesh.add_point(*it_set);
	}
	m_TetGenInterface.triangulate_pointset();
}

void Mesh_Intersection::update_intersection_mesh()
{
	// Test which elements must be added, and insert their vertices
	libMesh::SerialMesh::element_iterator   it_poly_mesh = m_libMesh_PolyhedronMesh.elements_begin();
	long dummy_long_int = 0;
	for(	; it_poly_mesh != m_libMesh_PolyhedronMesh.elements_end();
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
				dummy_long_int = convert_to_grid(poly_elem->point(iii));
				if(m_Grid_to_mesh_vertex_idx.find(dummy_long_int)==m_Grid_to_mesh_vertex_idx.end())
				{
					// New vertex! Add it to the mesh
					m_Grid_to_mesh_vertex_idx[dummy_long_int] = m_nb_of_vertices;
					m_libMesh_Mesh.add_point(poly_elem->point(iii),m_nb_of_vertices);
					mesh_elem->set_node(iii) = m_libMesh_Mesh.node_ptr(m_nb_of_vertices);
					++m_nb_of_vertices;
				}

				// Associate vertex to the new element
				mesh_elem->set_node(iii) = m_libMesh_Mesh.node_ptr(m_Grid_to_mesh_vertex_idx[dummy_long_int]);
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

const libMesh::SerialMesh & Mesh_Intersection::mesh()
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
	m_Grid_to_mesh_vertex_idx.clear();
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
	m_Grid_to_mesh_vertex_idx.reserve(map_preallocation);
}

void Mesh_Intersection::set_grid_constraints(const libMesh::Mesh & mesh_A, const libMesh::Mesh & mesh_B)
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

	m_vol_tol = std::pow(m_eps,2);

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

void Mesh_Intersection::increase_intersection_mesh(	const std::set<libMesh::Point> & input_points,
									unsigned int elem_idx_A,
									unsigned int elem_idx_B,
									unsigned int elem_idx_C)
{
	m_bMeshFinalized = false;
	unsigned int intersection_range_start = m_nb_of_elements;

	//	First, triangulate the point set!
	triangulate_intersection(input_points);

	//	Second, add the points to the grid-to-mesh map and update the intersection mesh,
	// taking into account if the associated elements are big enough.
	update_intersection_mesh();

/*	/i/ 	First, add the points to the grid-to-mesh map, and create new
	// vertices if needed.
	update_intersection_vertices(input_points);

	//	Second, triagulate the point set
	triangulate_intersection(input_points);

	//  Third, update the intersection mesh
	// 	The sets conserve the order, so there is no need of an equivalence
	// table between the m_libMesh_PolyhedronMesh and m_libMesh_Mesh
	// vertices.
	update_intersection_mesh();
*/
	//	Third, update the intersection pairs, if needed
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
	std::string table_file_out = filename_base + "_inter_table_Full.dat";

	// Print the mesh
	m_libMesh_Mesh.write(mesh_file_out);

	// Print the intersection table
	std::ofstream table_out(table_file_out);
	unsigned int nb_of_inter_elements = 0;
	table_out << m_nb_of_intersections << " " << m_libMesh_Mesh.n_elem() << " " << m_libMesh_Mesh.n_nodes() << std::endl;

	for(unsigned int iii = 0; iii < m_nb_of_intersections; ++iii)
	{
		nb_of_inter_elements = m_intersection_element_range[iii].second
							 - m_intersection_element_range[iii].first;
		table_out << iii << " " << m_intersection_pairs[iii].first << " "
				                << m_intersection_pairs[iii].second  << " "
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
	m_libMesh_Mesh.prepare_for_use();
	m_bMeshFinalized = true;
}

double Mesh_Intersection::get_total_volume()
{
	double volume = 0;

	homemade_assert_msg(m_bMeshFinalized,"Intersection mesh was not prepared for use!\n");

	for(libMesh::SerialMesh::element_iterator it_elem = m_libMesh_Mesh.elements_begin();
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

	for(libMesh::SerialMesh::element_iterator it_elem = m_libMesh_PolyhedronMesh.elements_begin();
			it_elem != m_libMesh_PolyhedronMesh.elements_end();
			++it_elem)
	{
		const libMesh::Elem * elem = * it_elem;
		volume += elem->volume();
	}

	return volume;
}

long Mesh_Intersection::convert_to_grid(const libMesh::Point iPoint)
{
	long dummy =  lround( (iPoint(0) -  m_Grid_MinPoint(0) )/m_eps) * m_GridN[1]*m_GridN[2]
				+ lround( (iPoint(1) -  m_Grid_MinPoint(1) )/m_eps) * m_GridN[1]
				+ lround( (iPoint(2) -  m_Grid_MinPoint(2) )/m_eps);
	homemade_assert_msg(dummy > -1, "Negative grid index!\n");

	return dummy;
}
}
