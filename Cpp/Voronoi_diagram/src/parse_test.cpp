/*
 * parse_test.cpp
 *
 *  Created on: Jan 29, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#include "voro++.hh"
#include "raw_data_parser.h"

int main(int argc, char *argv[])
{
	GetPot cmd_line(argc,argv);
	std::string input_file_1;
	std::string input_file_2;
	std::string output_file;

	double weight = 0.1;

	if ( cmd_line.search(1, "-i1") )
	{
		input_file_1 = cmd_line.next(input_file_1.c_str());
	}
	if ( cmd_line.search(1, "-i2") )
	{
		input_file_2 = cmd_line.next(input_file_2.c_str());
	}

	if ( cmd_line.search(1, "-o") )
	{
		output_file = cmd_line.next(output_file.c_str());
	}

	if ( cmd_line.search(1, "-w") )
	{
		weight = cmd_line.next(weight);
	}

	Raw_data_parser parser(static_cast<long>(1E9));
	std::vector<double> min_point = {-1,-0.5,-0.5};
	std::vector<double> max_point = {1,0.5,0.5};

	parser.set_grid_constraints(min_point,max_point);

	parser.prepare_table_dimensions(input_file_1);
	parser.prepare_table_dimensions(input_file_2);

	parser.prepare_tables();

	parser.add_file_data(input_file_1);
	parser.add_file_data(input_file_2);

	parser.print_points(output_file,weight);

	parser.print_cells(input_file_1,output_file);
	parser.print_cells(input_file_2,output_file);

	parser.print_physical_groups(output_file);

	return 0;
}


