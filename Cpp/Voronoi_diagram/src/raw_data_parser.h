/*
 * raw_data_parser.h
 *
 *  Created on: Jan 29, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef RAW_DATA_PARSER_H_
#define RAW_DATA_PARSER_H_

#include "voro_functions.h"
#include <map>

template<typename T>
void jump_lines(T& filestream, unsigned int numberOfLines = 1)
{
	std::string dummy;
	for(unsigned int iii = 0; iii < numberOfLines; ++iii)
		std::getline(filestream,dummy);
};

class Raw_data_parser
{
private:
	int nb_of_vertices_upper_limit;
	int nb_of_faces;
	int nb_of_edges;
	int nb_of_cells;

	int nb_of_vertices;
	int nb_of_grains;

	std::vector<double> m_Grid_MinPoint;
	std::vector<double> m_Grid_MaxPoint;
	long m_GridN_min;
	std::vector<long> m_GridN;
	double m_vol_tol;
	double m_eps;

	double m_weight;

	std::unordered_map<std::vector<long>, unsigned int, PointHash_3D, PointHash_3D_Equal >
		m_discrete_vertices;

	std::unordered_map<unsigned int,std::vector<double> >
		m_tess_vertices;

	std::multimap<int,int> m_grain_mapping;
	std::unordered_set<int> m_grain_set;

	Raw_data_parser()
	{

	};

public:
	Raw_data_parser(long grid_n_min);
	virtual ~Raw_data_parser();

	void set_grid_constraints(std::vector<double> Grid_MinPoint, std::vector<double> Grid_MaxPoint);
	void convert_to_discrete(std::vector<double>& Grid_MinPoint, std::vector<long>& oPoint);

	void prepare_table_dimensions(std::string& filename);
	void prepare_tables();

	void add_file_data(std::string& filename);

	void print_points(std::string filename, double weight);
	void print_cells(std::string& input_data_filename, std::string& output_data_filename);
	void print_edge(int edge_start, int edge_end, std::ofstream& output_stream);
	void print_face(std::vector<int>& face_edges, std::ofstream& output_stream);
	void print_single_cell(std::vector<int>& face_idx, std::ofstream& output_stream);

	void print_physical_groups(std::string filename);


};
#endif /* RAW_DATA_PARSER_H_ */

