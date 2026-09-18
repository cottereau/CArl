/*
 * libmesh_update_init_cond_input_parser.h
 *
 *  Created on: August 5, 2026
 *      Author: Romain Ruyssen
 */

#include "libmesh_update_init_cond_input_parser.h"
void get_input_params(GetPot& field_parser,
		libmesh_update_init_cond_input_params& input_params) {

	// Set time mesh file
	if ( field_parser.search(1, "TimeMeshFile") )
	{
		input_params.time_mesh_file = field_parser.next(input_params.time_mesh_file);
	}
	else
	{
		homemade_error_msg("Missing the Time Mesh file!");
	}

    // Set space mesh file
    if ( field_parser.search(1, "SpaceMeshFile") )
    {
        input_params.space_mesh_file = field_parser.next(input_params.space_mesh_file);
    }
    else
    {
        homemade_error_msg("Missing the Space Mesh file!");
    }

	// Set Newmark parameters
	if ( field_parser.search(1, "NewmarkParameters") )
	{
		input_params.newmark_params_file = field_parser.next(input_params.newmark_params_file);
	}
	else
	{
		homemade_error_msg("Missing the Newmark parameters file!");
	}

	// Set solution file
    if ( field_parser.search(1, "Solution") )
    {
        input_params.solution_file = field_parser.next(input_params.solution_file);
    }
    else
    {
        homemade_error_msg("Missing the solution file!");
    }   

    // Set initial conditions files
    if ( field_parser.search(1, "InitCond_disp") )
    {
        input_params.init_cond_disp_file = field_parser.next(input_params.init_cond_disp_file);
    }
    else
    {
        homemade_error_msg("Missing the initial displacement file!");
    }

    if ( field_parser.search(1, "InitCond_vel") )
    {
        input_params.init_cond_vel_file = field_parser.next(input_params.init_cond_vel_file);
    }
    else
    {
        homemade_error_msg("Missing the initial velocity file!");
    }

    if ( field_parser.search(1, "InitCond_acc") )
    {
        input_params.init_cond_acc_file = field_parser.next(input_params.init_cond_acc_file);
    }
    else
    {
        homemade_error_msg("Missing the initial acceleration file!");
    }
};
