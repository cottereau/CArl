/*
 * raw_data_parser.cpp
 *
 *  Created on: Jan 29, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#include "raw_data_parser.h"

Raw_data_parser::Raw_data_parser(long grid_n_min) :
	nb_of_vertices_upper_limit { 0 },
	nb_of_faces { 0 },
	nb_of_edges { 0 },
	nb_of_cells { 0 },
	nb_of_vertices { 0 },
	nb_of_grains { 0 },
	m_GridN_min { grid_n_min },
	m_vol_tol { -1 },
	m_eps { -1 }
{
	m_Grid_MinPoint.resize(3,-1);
	m_Grid_MaxPoint.resize(3,-1);
	m_GridN.resize(3,-1);
}

Raw_data_parser::~Raw_data_parser()
{
	// TODO Auto-generated destructor stub
}

void Raw_data_parser::set_grid_constraints(std::vector<double> Grid_MinPoint, std::vector<double> Grid_MaxPoint)
{
	// Set the (future) intersection bbox corners
	std::vector<double> eps_candidates(3,0);

	m_Grid_MinPoint = Grid_MinPoint;
	m_Grid_MaxPoint = Grid_MaxPoint;
	for(unsigned int iii = 0; iii < 3; ++iii)
	{
		eps_candidates[iii] = (m_Grid_MaxPoint[iii] - m_Grid_MinPoint[iii])/m_GridN_min;
	}

	m_eps = *std::min_element(eps_candidates.begin(),eps_candidates.end());

	for(unsigned int iii = 0; iii < 3; ++iii)
	{
		m_Grid_MinPoint[iii] -= 2*m_eps;
		m_Grid_MaxPoint[iii] += 2*m_eps;
	}

	for(unsigned int iii = 0; iii < 3; ++iii)
	{
		m_GridN[iii] = (m_Grid_MaxPoint[iii] - m_Grid_MinPoint[iii]) / m_eps + 1;
	}

		std::cout << "    DEBUG: discrete grid" << std::endl;
		std::cout << " -> eps             : " << m_eps << std::endl;
		std::cout << " -> volume          : " << m_vol_tol << std::endl;
		std::cout << " -> Grid dimensions : " << m_GridN[0] << " " << m_GridN[1] << " " << m_GridN[2] << " " << std::endl  << std::endl;
}

void Raw_data_parser::convert_to_discrete(std::vector<double>& iPoint, std::vector<long>& oPoint)
{
	oPoint[0] = lround( (iPoint[0] -  m_Grid_MinPoint[0] )/m_eps);
	oPoint[1] = lround( (iPoint[1] -  m_Grid_MinPoint[1] )/m_eps);
	oPoint[2] = lround( (iPoint[2] -  m_Grid_MinPoint[2] )/m_eps);
}

void Raw_data_parser::prepare_table_dimensions(std::string& filename)
{
	// Second, build the point hash
	std::ifstream input_stream(filename);

	// Data
	int local_nb_of_vertices = 0;
	int local_nb_of_faces = 0;
	int local_nb_of_edges = 0;

	int dummy_data = 0;
	while(true)
	{

		/* File format:
		 * 		cell_id nb_of_vertices nb_of_faces
		 * 		[list of vertices]
		 * 		[nb of vertices per face]
		 * 		[list of faces]
		 *
		 * 		cell_id nb_of_vertices nb_of_faces
		 * 		[list of vertices]
		 * 		[nb of vertices per face]
		 * 		[list of faces]
		 *
		 * 		...
		 */

		input_stream >> dummy_data;

		if( input_stream.eof() ) break;

		input_stream >> local_nb_of_vertices;
		input_stream >> local_nb_of_faces;

		jump_lines(input_stream,local_nb_of_vertices + 1);

		local_nb_of_edges = 0;
		for(int iii = 0; iii < local_nb_of_faces; ++iii)
		{
			input_stream 	>> dummy_data;
			local_nb_of_edges += dummy_data;
		}

		jump_lines(input_stream,local_nb_of_faces + 2);

		nb_of_vertices_upper_limit += local_nb_of_vertices;
	}
	input_stream.close();
};

void Raw_data_parser::prepare_tables()
{
	m_discrete_vertices.reserve(nb_of_vertices_upper_limit);
	m_tess_vertices.reserve(nb_of_vertices_upper_limit);
}

void Raw_data_parser::add_file_data(std::string& filename)
{
	// First, check that the read is correct:
	//    only read the data and print it -> OK!
	//
	// Second, build the point hash
	std::ifstream input_stream(filename);
	std::ofstream debug_stream("debug.dat");

	// Data
	int cell_phys_id = 0;
	int cell_nb_of_vertices;
	int cell_nb_of_faces;

	std::vector<int> cell_vertex_idx;
	std::vector<std::vector<double> > cell_vertex;
	std::vector<long> vertex_discrete_point(3,0);

	std::vector<int> face_nb_of_vertices;
	std::vector<std::vector<int> > face_vertex_list;

	//
	while(true)
	{

		/* File format:
		 * 		cell_id nb_of_vertices nb_of_faces
		 * 		[list of vertices]
		 * 		[nb of vertices per face]
		 * 		[list of faces]
		 *
		 * 		cell_id nb_of_vertices nb_of_faces
		 * 		[list of vertices]
		 * 		[nb of vertices per face]
		 * 		[list of faces]
		 *
		 * 		...
		 */

		input_stream >> cell_phys_id >> cell_nb_of_vertices >> cell_nb_of_faces;

		if( input_stream.eof() ) break;

		cell_vertex_idx.resize(cell_nb_of_vertices);
		cell_vertex.resize(cell_nb_of_vertices,std::vector<double>(3,0));

		face_nb_of_vertices.resize(cell_nb_of_faces);
		face_vertex_list.resize(cell_nb_of_faces);

		for(int iii = 0; iii < cell_nb_of_vertices; ++iii)
		{
			cell_vertex_idx[iii] = iii;
			input_stream 	>> cell_vertex[iii][0]
							>> cell_vertex[iii][1]
							>> cell_vertex[iii][2];

			this->convert_to_discrete(cell_vertex[iii],vertex_discrete_point);
			if(m_discrete_vertices.find(vertex_discrete_point) == m_discrete_vertices.end())
			{
				// New vertex! Add it to the tesselation
				m_discrete_vertices[vertex_discrete_point] = nb_of_vertices;
				m_tess_vertices[nb_of_vertices] = cell_vertex[iii];
				++nb_of_vertices;
			}
		}

		for(int iii = 0; iii < cell_nb_of_faces; ++iii)
		{
			input_stream 	>> face_nb_of_vertices[iii];
			face_vertex_list[iii].resize(face_nb_of_vertices[iii]);
		}

		for(int iii = 0; iii < cell_nb_of_faces; ++iii)
		{
			for(int jjj = 0; jjj < face_nb_of_vertices[iii]; ++jjj)
			{
				input_stream 	>> face_vertex_list[iii][jjj];
			}
		}

		std::string dummy;
		std::getline(input_stream,dummy);

		// DEBUG check all this
		debug_stream << cell_phys_id << " "
				     << cell_nb_of_vertices << " "
					 << cell_nb_of_faces << std::endl;

		for(int iii = 0; iii < cell_nb_of_vertices; ++iii)
		{
			debug_stream << cell_vertex[iii][0] << " "
						 << cell_vertex[iii][1] << " "
						 << cell_vertex[iii][2] << std::endl;
		}

		for(int iii = 0; iii < cell_nb_of_faces; ++iii)
		{
			debug_stream << face_nb_of_vertices[iii] << " ";
		}
		debug_stream << std::endl;

		for(int iii = 0; iii < cell_nb_of_faces; ++iii)
		{
			for(int jjj = 0; jjj < face_nb_of_vertices[iii]; ++jjj)
			{
				debug_stream << face_vertex_list[iii][jjj] << " ";
			}
			debug_stream << std::endl;
		}
		debug_stream << std::endl;
	}
	input_stream.close();
	debug_stream.close();
}

void Raw_data_parser::print_points(std::string filename, double weight)
{
	std::ofstream gmshOutput(filename, std::ofstream::trunc);

	m_weight = weight;

	gmshOutput << "// Weights" << std::endl;
	gmshOutput << "w1 = " << weight << ";" << std::endl << std::endl;

	gmshOutput << "// Points" << std::endl;
	for(int iii = 0; iii < nb_of_vertices; ++iii)
	{
		gmshOutput 	<< "Point(" << iii + 1 << ") = {";
		for(int jjj = 0; jjj < 2; ++jjj)
		{
			gmshOutput << m_tess_vertices[iii][jjj] << ", ";
		}
		gmshOutput << m_tess_vertices[iii][2] << ", w1};" << std::endl;
	}
	gmshOutput << std::endl;
	gmshOutput.close();
}

void Raw_data_parser::print_cells(std::string& input_data_filename, std::string& output_data_filename)
{
	// Build data!
	std::ifstream input_stream(input_data_filename);
	std::ofstream output_stream(output_data_filename,std::ofstream::app);

	// Data
	int cell_phys_id = 0;
	int cell_nb_of_vertices;
	int cell_nb_of_faces;

	std::vector<int> cell_vertex_idx;
	std::vector<double> vertex_continuous_point(3,0);
	std::vector<long> vertex_discrete_point(3,0);

	std::vector<int> face_nb_of_vertices;
	std::vector<std::vector<int> > face_vertex_list;
	std::vector<int> face_edges;
	std::vector<int> face_idx;

	int edge_start = 0;
	int edge_end = 0;

	while(true)
	{

		/* File format:
		 * 		cell_id nb_of_vertices nb_of_faces
		 * 		[list of vertices]
		 * 		[nb of vertices per face]
		 * 		[list of faces]
		 *
		 * 		cell_id nb_of_vertices nb_of_faces
		 * 		[list of vertices]
		 * 		[nb of vertices per face]
		 * 		[list of faces]
		 *
		 * 		...
		 */

		input_stream >> cell_phys_id >> cell_nb_of_vertices >> cell_nb_of_faces;

		if( input_stream.eof() ) break;

		// Add cell to list that will be used later to build the physical grains
		if(m_grain_set.find(cell_phys_id) == m_grain_set.end())
		{
			// new grain!
			m_grain_set.insert(cell_phys_id);
		}

		++nb_of_cells;
		m_grain_mapping.insert(std::make_pair(cell_phys_id,nb_of_cells));

		cell_vertex_idx.resize(cell_nb_of_vertices);

		face_nb_of_vertices.resize(cell_nb_of_faces);
		face_vertex_list.resize(cell_nb_of_faces);
		face_idx.resize(cell_nb_of_faces);

		for(int iii = 0; iii < cell_nb_of_vertices; ++iii)
		{
			input_stream 	>> vertex_continuous_point[0]
							>> vertex_continuous_point[1]
							>> vertex_continuous_point[2];

			this->convert_to_discrete(vertex_continuous_point,vertex_discrete_point);
			if(m_discrete_vertices.find(vertex_discrete_point) == m_discrete_vertices.end())
			{
				std::cout << "ERROR" << std::endl;
			}
			cell_vertex_idx[iii] = m_discrete_vertices[vertex_discrete_point];
		}

		for(int iii = 0; iii < cell_nb_of_faces; ++iii)
		{
			input_stream 	>> face_nb_of_vertices[iii];
			face_vertex_list[iii].resize(face_nb_of_vertices[iii]);
		}

		for(int iii = 0; iii < cell_nb_of_faces; ++iii)
		{
			for(int jjj = 0; jjj < face_nb_of_vertices[iii]; ++jjj)
			{
				input_stream 	>> face_vertex_list[iii][jjj];
			}

			face_edges.resize(face_nb_of_vertices[iii],0);

			for(int jjj = 0; jjj < face_nb_of_vertices[iii] - 1; ++jjj)
			{
				edge_start = cell_vertex_idx[face_vertex_list[iii][jjj]];
				edge_end   = cell_vertex_idx[face_vertex_list[iii][jjj+1]];

				this->print_edge(edge_start,edge_end,output_stream);
				face_edges[jjj] = nb_of_edges;
			}

			edge_start = cell_vertex_idx[face_vertex_list[iii][face_nb_of_vertices[iii] - 1]];
			edge_end   = cell_vertex_idx[face_vertex_list[iii][0]];

			this->print_edge(edge_start,edge_end,output_stream);
			face_edges[face_nb_of_vertices[iii] - 1] = nb_of_edges;

			this->print_face(face_edges,output_stream);
			face_idx[iii] = nb_of_faces;
		}

		print_single_cell(face_idx,output_stream);

		std::string dummy;
		std::getline(input_stream,dummy);
	}
	input_stream.close();

	output_stream.close();
}

void Raw_data_parser::print_edge(int edge_start, int edge_end, std::ofstream& output_stream)
{
	double 	edge_length = std::sqrt(std::pow(m_tess_vertices[edge_start][0] - m_tess_vertices[edge_end][0],2.) +
			std::pow(m_tess_vertices[edge_start][1] - m_tess_vertices[edge_end][1],2.) +
			std::pow(m_tess_vertices[edge_start][2] - m_tess_vertices[edge_end][2],2.));

	int nb_of_points_per_edge = 0;

	if(edge_length > m_weight)
	{
		nb_of_points_per_edge = std::ceil(edge_length / m_weight) + 1;
	}
	else
	{
		nb_of_points_per_edge = 2;
	}

	// Add the edge
	++nb_of_edges;
	output_stream	<< "Line(" << nb_of_edges << ") = {"
					<< edge_start + 1 << ", " << edge_end + 1 << "};" << std::endl;
	output_stream  << "Transfinite Line(" << nb_of_edges << ") = "
				   << nb_of_points_per_edge << " Using Bump 0.25;" << std::endl;
}

void Raw_data_parser::print_face(std::vector<int>& face_edges, std::ofstream& output_stream)
{
	++nb_of_faces;
	output_stream	<< "Line Loop(" << nb_of_faces << ") = {";

	for(unsigned int jjj = 0; jjj < face_edges.size() - 1; ++jjj)
	{
		output_stream	<< face_edges[jjj] << ", ";
	}

	// -> Close the face
	output_stream	<< face_edges[face_edges.size() - 1] << "};" << std::endl;
	output_stream	<< "Plane Surface(" <<  nb_of_faces << ") = {" <<  nb_of_faces << "};" << std::endl;
}

void Raw_data_parser::print_single_cell(std::vector<int>& face_idx, std::ofstream& output_stream)
{
	output_stream << "Surface Loop(" << nb_of_cells << ") = {";
	for(unsigned int jjj = 0; jjj < face_idx.size() - 1; ++jjj)
	{
		output_stream	<< face_idx[jjj] << ", ";
	}

	// -> Close the cell
	output_stream	<< face_idx[face_idx.size() - 1] << "};" << std::endl;

	// -> Set up the volumes
	output_stream	<< "Volume(" <<  nb_of_cells
				<< ") = {" <<  nb_of_cells << "};" << std::endl;

}

void Raw_data_parser::print_physical_groups(std::string filename)
{
	std::ofstream output_stream(filename,std::ofstream::app);

	output_stream << "Coherence;" << std::endl;

	// Print the physical grains
	auto grain_set_it = m_grain_set.begin();
	auto grain_set_it_end = m_grain_set.end();

	for( ; grain_set_it != grain_set_it_end; ++grain_set_it)
	{
		auto range_cells = m_grain_mapping.equal_range(*grain_set_it);
		std::multimap<int,int>::iterator it_end = range_cells.second;
		--it_end;

		output_stream	<< "Physical Volume(" <<  *grain_set_it
					<< ") = {";
		for(auto it = range_cells.first; it != it_end; ++it)
		{
			output_stream << it->second << ",";
		}
		output_stream << it_end->second << "};" << std::endl;
	}

	output_stream.close();
}

