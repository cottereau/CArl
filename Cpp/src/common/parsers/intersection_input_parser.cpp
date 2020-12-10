/*
 * intersection_input_parser.h
 *
 *  Created on: Feb 14, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#include "intersection_input_parser.h"

namespace carl
{
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
