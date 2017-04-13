/*
 * intersection_input_parser.h
 *
 *  Created on: Feb 14, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef INTERSECTION_INPUT_PARSER_H_
#define INTERSECTION_INPUT_PARSER_H_

/**
 @file intersection_input_parser.h
 */
#include "carl_headers.h"

namespace carl
{
/**	\brief Structure containing the parameters for the parallel intersection search test program (source: CArl_build_intersections.cpp)
 *	
 *		Details on the parameters setup are found in the documentation of carl::get_intersection_input_params(GetPot& field_parser,
 *		parallel_intersection_params& input_params).
 */
struct parallel_intersection_params {
	std::string mesh_A;	///< Mesh A path.
	std::string mesh_B; ///< Mesh B path.
	std::string mesh_C; ///< Coupling mesh path.
	std::string output_base; ///< Output filename base.

	bool bStitchInterMeshes;	///< Stitch the intersection meshes?
	bool bVerbose;				///< Print coupling partitioning?
	
	carl::IntersectionMeshingMethod inter_meshing_method; ///< Intersection meshing method. *Values*: carl::IntersectionMeshingMethod::CGAL or carl::IntersectionMeshingMethod::LIBMESH_TETGEN.
};

/**	\brief Parser function for the parallel intersection search program (source: CArl_build_intersections.cpp)
 *	
 *	Required parameters:
 *	- `MeshA`, `-mA` or `--meshA` : path to the mesh A.
 *	- `MeshB`, `-mB` or `--meshB` : path to the mesh B.
 *	- `MeshC`, `-mC` or `--meshC` : path to the coupling mesh C.
 *
 *  Optional parameters:
 *	- `OutputBase`, `-mO` or `--output` : base of the output files (including folders). *Default*: `test_inter`.
 *  - `MeshingMethod` or `--meshingMethodType` : intersection meshing method. *Values*: `CGAL` or `LIBMESH_TETGEN`. *Default*: `CGAL`.
 *  
 *  Boolean flags:
 *  - `StitchInterMeshes` : do not stich together the intersection meshes. 
 *  - `VerboseOutput` or `--verbose` : print some extra information, such as the coupling mesh partitioning.
 */
void get_intersection_input_params(GetPot& field_parser,
		parallel_intersection_params& input_params) {

	// Set mesh files
	if (field_parser.search(3, "--meshA", "-mA", "MeshA")) {
		input_params.mesh_A = field_parser.next(
				input_params.mesh_A);
	} else {
		homemade_error_msg("Missing mesh A!"); 
	}

	if (field_parser.search(3, "--meshB", "-mB", "MeshB")) {
		input_params.mesh_B = field_parser.next(
				input_params.mesh_B);
	} else {
		homemade_error_msg("Missing mesh B!"); 
	}

	if (field_parser.search(3, "--meshC", "-mC", "MeshC")) {
		input_params.mesh_C = field_parser.next(
				input_params.mesh_C);
	} else {
		homemade_error_msg("Missing the coupling mesh C!"); 
	}

	if (field_parser.search(3, "--output", "-mO", "OutputBase")) {
		input_params.output_base = field_parser.next(
				input_params.output_base);
	} else {
		input_params.output_base = "inter_search";
	}

	std::string meshing_method;
	input_params.inter_meshing_method = carl::IntersectionMeshingMethod::CGAL;
	if (field_parser.search(2, "--meshingMethodType", "MeshingMethod")) {
		meshing_method = field_parser.next(
				meshing_method);
		if(meshing_method == "CGAL")
		{
			input_params.inter_meshing_method = carl::IntersectionMeshingMethod::CGAL;
		}
		else if(meshing_method == "TETGEN" )
		{
			input_params.inter_meshing_method = carl::IntersectionMeshingMethod::LIBMESH_TETGEN;
		}
	}

	if(field_parser.search(1,"StitchInterMeshes")) {
		input_params.bStitchInterMeshes = true;
	}
	else
	{
		input_params.bStitchInterMeshes = false;
	}

	if(field_parser.search(2,"VerboseOutput", "--verbose")) {
		input_params.bVerbose = true;
	}
	else
	{
		input_params.bVerbose = false;
	}
};

};
#endif /* INTERSECTION_INPUT_PARSER_H_ */
