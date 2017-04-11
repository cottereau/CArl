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
/**	\brief Structure containing the parameters for the parallel intersection search test program (source: parallel_intersection_test.cpp)
 *	
 *		Details on the parameters setup are found in the documentation of carl::get_input_params(GetPot& field_parser,
 *		parallel_intersection_test_params& input_params).
 */
struct parallel_intersection_test_params {
	std::string mesh_A;	///< Mesh A path.
	std::string mesh_B; ///< Mesh B path.
	std::string mesh_C; ///< Coupling mesh path.
	std::string output_base; ///< Output filename base.

	carl::SearchMethod search_type;	///< Search type (either CGAL or LIBMESH_TETGEN)

	bool bSkipIntersectionConstruction;	///< Skip the construction of the intersections? *Default*: false.
	bool bSkipIntersectionPartitioning;	///< Skip the intersection partitioning? *Default*: false.
	bool bSkipRestriction;	///< Skip the mesh restrictions? *Default*: false.
	bool bSkipMeshStitching;	///< Skip the intersection mesh stitching? *Default*: false.
	bool bExportScalingData;	///< Export the scaling data? *Default*: false.

	carl::IntersectionMeshingMethod inter_meshing_method;
};

/**	\brief Parser function for the parallel intersection search test program (source: parallel_intersection_test.cpp)
 *	
 *	Required parameters:
 *	- `MeshA`, `-mA` or `--meshA` : path to the mesh A.
 *	- `MeshB`, `-mB` or `--meshB` : path to the mesh B.
 *	- `MeshC`, `-mC` or `--meshC` : path to the coupling mesh C.
 *
 *  Optional parameters:
 *	- `OutputBase`, `-mO` or `--output` : base of the output files (including folders). *Default*: `test_inter`.
 * 	- `SearchType` or `--searchType`    : search type for the patch intersections. *Values*: `FRONT`, `BRUTE` or `BOTH`. *Default*: `BRUTE`.
 *  - `MeshingMethod` or `--meshingMethodType` : intersection meshing method. *Values*: `CGAL` or `LIBMESH_TETGEN`. *Default*: `CGAL`.
 *  
 *  Boolean flags:
 *  - `SkipIntersectionConstruction` : skip the construction of the intersections (useful to test the search times).
 *  - `SkipIntersectionPartitioning` : skip the repartitioning of the intersection construction.
 *  - `SkipRestriction` : do not build meshes' A and B restriction to the coupling region.
 *  - `ExportScalingData` : export the data for the algorithm scaling.
 *  - `SkipMeshStitching` : do not stich together the intersection meshes. 
 */
void get_input_params(GetPot& field_parser,
		parallel_intersection_test_params& input_params) {

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
		input_params.output_base = "test_inter";
	}

	std::string search_type;
	if (field_parser.search(2, "--searchType", "SearchType")) {
		search_type = field_parser.next(
				search_type);
		if(search_type == "Front" || search_type == "front" || search_type == "FRONT")
		{
			input_params.search_type = carl::SearchMethod::FRONT;
		}
		else if(search_type == "Brute" || search_type == "brute" || search_type == "BRUTE")
		{
			input_params.search_type = carl::SearchMethod::BRUTE;
		}
		else if(search_type == "Both" || search_type == "both" || search_type == "BOTH")
		{
			input_params.search_type = carl::SearchMethod::BOTH;
		}
		else
		{
			input_params.search_type = carl::SearchMethod::BRUTE;
		}
	}
	else
	{
		input_params.search_type = carl::SearchMethod::BRUTE;
	}

	std::string meshing_method;
	if (field_parser.search(2, "--meshingMethodType", "MeshingMethod")) {
		meshing_method = field_parser.next(
				search_type);
		if(meshing_method == "CGAL")
		{
			input_params.inter_meshing_method = carl::IntersectionMeshingMethod::CGAL;
		}
		else if(meshing_method == "TETGEN" )
		{
			input_params.inter_meshing_method = carl::IntersectionMeshingMethod::LIBMESH_TETGEN;
		}
		else
		{
			input_params.inter_meshing_method = carl::IntersectionMeshingMethod::LIBMESH_TETGEN;
		}
	}
	else
	{
		input_params.inter_meshing_method = carl::IntersectionMeshingMethod::CGAL;
	}

	if(field_parser.search(1,"SkipIntersectionConstruction")) {
		input_params.bSkipIntersectionConstruction = true;
	}
	else
	{
		input_params.bSkipIntersectionConstruction = false;
	}

	if(field_parser.search(1,"SkipIntersectionPartitioning")) {
		input_params.bSkipIntersectionPartitioning = true;
	}
	else
	{
		input_params.bSkipIntersectionPartitioning = false;
	}

	if(field_parser.search(1,"SkipRestriction")) {
		input_params.bSkipRestriction = true;
	}
	else
	{
		input_params.bSkipRestriction = false;
	}

	if(field_parser.search(1,"ExportScalingData")) {
		input_params.bExportScalingData = true;
	}
	else
	{
		input_params.bExportScalingData = false;
	}

	if(field_parser.search(1,"SkipMeshStitching")) {
		input_params.bSkipMeshStitching = true;
	}
	else
	{
		input_params.bSkipMeshStitching = false;
	}
};

};
#endif /* INTERSECTION_INPUT_PARSER_H_ */
