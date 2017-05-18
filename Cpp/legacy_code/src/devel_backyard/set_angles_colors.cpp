/*
 * set_angles_colors.cpp
 *
 *  Created on: Sep 6, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#include "set_angles_colors.h"

int main(int argc, char** argv)
{
	libMesh::LibMeshInit init(argc, argv);

	GetPot command_line(argc, argv);

	std::string mesh_file;
	std::string angles_file;
	std::string mesh_out_file;

	// Set mesh files
	if (command_line.search(1, "--meshI")) {
		mesh_file = command_line.next(
				mesh_file);
	} else {
		homemade_error_msg("Missing the input mesh file!");
	}

	if (command_line.search(1, "--angles")) {
		angles_file = command_line.next(
				angles_file);
	} else {
		homemade_error_msg("Missing the angles file!");
	}

	if (command_line.search(1, "--meshO")) {
		mesh_out_file = command_line.next(
				mesh_out_file);
	} else {
		homemade_error_msg("Missing the output mesh file!");
	}

	libMesh::Parallel::Communicator& WorldComm = init.comm();

	libMesh::Mesh mesh_input(WorldComm, 3);
	mesh_input.allow_renumbering(false);

	libMesh::ExodusII_IO exodus_output(mesh_input);
	exodus_output.read(mesh_file);
	mesh_input.prepare_for_use();

	libMesh::EquationSystems original_systems(mesh_input);
	libMesh::ExplicitSystem& new_data = original_systems.add_system<libMesh::ExplicitSystem>("PhysicalConstants");

	const std::vector<std::string>& elem_var_names = exodus_output.get_elem_var_names();
	const std::vector<std::string>& node_var_names = exodus_output.get_nodal_var_names();

	for(unsigned int iii = 0; iii < elem_var_names.size(); ++iii)
	{
		new_data.add_variable(elem_var_names[iii],libMesh::CONSTANT, libMesh::MONOMIAL);
	}

	for(unsigned int iii = 0; iii < node_var_names.size(); ++iii)
	{
		new_data.add_variable(node_var_names[iii],libMesh::CONSTANT, libMesh::MONOMIAL);
	}

	new_data.init();

//	new_data.add_variable("Index", libMesh::CONSTANT, libMesh::MONOMIAL);
//	new_data.add_variable("Angle_x", libMesh::CONSTANT, libMesh::MONOMIAL);
//	new_data.add_variable("Angle_y", libMesh::CONSTANT, libMesh::MONOMIAL);
//	new_data.add_variable("Angle_z", libMesh::CONSTANT, libMesh::MONOMIAL);
//	new_data.add_variable("color_r", libMesh::CONSTANT, libMesh::MONOMIAL);
//	new_data.add_variable("color_g", libMesh::CONSTANT, libMesh::MONOMIAL);
//	new_data.add_variable("color_b", libMesh::CONSTANT, libMesh::MONOMIAL);

	for(unsigned int iii = 0; iii < elem_var_names.size(); ++iii)
	{
		exodus_output.copy_elemental_solution(new_data,elem_var_names[iii],elem_var_names[iii]);
	}

	for(unsigned int iii = 0; iii < node_var_names.size(); ++iii)
	{
		exodus_output.copy_nodal_solution(new_data,node_var_names[iii],node_var_names[iii]);
	}

//	double BIG_E = 0,BIG_Mu = 0;
//	carl::anisotropic_elasticity_tensor_cubic_sym anisotropy_data(original_systems,angles_file,BIG_E,BIG_Mu);



//	original_data.add_variable("u",libMesh::FIRST, libMesh::LAGRANGE);
//	original_data.add_variable("v",libMesh::FIRST, libMesh::LAGRANGE);
//	original_data.add_variable("w",libMesh::FIRST, libMesh::LAGRANGE);
//
//	exodus_input.copy_nodal_solution(original_data,"u","u");
//	exodus_input.copy_nodal_solution(original_data,"v","v");
//	exodus_input.copy_nodal_solution(original_data,"w","w");

//	for(unsigned int iii = 0; iii < elem_var_names.size(); ++iii)
//	{
//		original_data.add_variable(elem_var_names[iii],libMesh::CONSTANT, libMesh::MONOMIAL);
//		exodus_input.copy_elemental_solution(original_data,elem_var_names[iii],elem_var_names[iii]);
//	}

//	for(unsigned int iii = 0; iii < node_var_names.size(); ++iii)
//	{
//		original_data.add_variable(node_var_names[iii],libMesh::CONSTANT, libMesh::MONOMIAL);
//		exodus_input.copy_nodal_solution(original_data,node_var_names[iii],node_var_names[iii]);
//	}

	std::set<std::string> system_names;
	system_names.insert("PhysicalConstants");
	exodus_output.write_equation_systems(mesh_out_file,original_systems,&system_names);

	return 0;
}


