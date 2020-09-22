/*
 * carl_assemble_coupling_input_parser.h
 *
 *  Created on: Feb 14, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#include "carl_assemble_coupling_input_parser.h"

namespace carl
{

void get_assemble_coupling_input_params(GetPot& field_parser,
    coupling_assemble_coupling_input_params& input_params) {

  // Set mesh files
  if (field_parser.search(3, "--meshA", "-mA", "MeshA")) {
    input_params.mesh_BIG_file = field_parser.next(
        input_params.mesh_BIG_file);
  } else {
    homemade_error_msg("Missing the A mesh file!");
  }

  if (field_parser.search(3, "--meshB", "-mB", "MeshB")) {
    input_params.mesh_micro_file = field_parser.next(
        input_params.mesh_micro_file);
  } else {
    homemade_error_msg("Missing the B mesh file!");
  }

  if (field_parser.search(3, "--meshI", "-mI", "InterBase")) {
    input_params.common_inter_file = field_parser.next(
        input_params.common_inter_file);
  } else {
    homemade_error_msg("Missing the path to the intersection files!");
  }

  // Set coupling parameters
  if( field_parser.search(2, "--ce","CouplingWidth") )
  {
    input_params.coupling_width = field_parser.next(input_params.coupling_width);
  } else {
    homemade_error_msg("Missing the coupling region width!");
  }

  if( field_parser.search(2, "--ck","CouplingRigidity") )
  {
    input_params.coupling_rigidity = field_parser.next(input_params.coupling_rigidity);
  } else {
    homemade_error_msg("Missing the coupling rigidity!");
  }

  // Output
  if (field_parser.search(3, "--output", "-mO", "OutputFolder"))
  {
    input_params.output_folder = field_parser.next(
      input_params.output_folder);
  } else {
    input_params.output_folder = "";
  }

  if (field_parser.search(3, "--meshAR", "-mAR", "Mesh_A_Restriction")) {
    input_params.mesh_restrict_BIG_file = field_parser.next(
        input_params.mesh_restrict_BIG_file);
  } else {
    input_params.mesh_restrict_BIG_file = input_params.common_inter_file + "_A_restriction.msh";
  }

  if (field_parser.search(3, "--meshBR", "-mBR", "Mesh_B_Restriction")) {
    input_params.mesh_restrict_micro_file = field_parser.next(
        input_params.mesh_restrict_micro_file);
  } else {
    input_params.mesh_restrict_micro_file = input_params.common_inter_file + "_B_restriction.msh";
  }

  // Set the equivalence tables
  if (field_parser.search(2, "--tableRA", "Mesh_A_RestrictionEquivalenceTable")) {
    input_params.equivalence_table_restrict_BIG_file = field_parser.next(
        input_params.equivalence_table_restrict_BIG_file);
  } else {
    input_params.equivalence_table_restrict_BIG_file = input_params.common_inter_file + "_A_restriction_restrict.dat";
  }

  if (field_parser.search(2, "--tableRB", "Mesh_B_RestrictionEquivalenceTable")) {
    input_params.equivalence_table_restrict_micro_file = field_parser.next(
        input_params.equivalence_table_restrict_micro_file);
  } else {
    input_params.equivalence_table_restrict_micro_file = input_params.common_inter_file + "_B_restriction_restrict.dat";
  }

  // Set the mediator mesh
  input_params.mediator_type = carl::MediatorType::USE_MACRO;
  input_params.mesh_mediator_file = input_params.mesh_restrict_BIG_file;

  std::string mediator_type;
  if (field_parser.search(1,"MediatorMesh"))
  {
    mediator_type = field_parser.next(mediator_type);

    if(mediator_type == "UseRestricted_A")
    {
      input_params.mediator_type = carl::MediatorType::USE_MACRO;
      input_params.mesh_mediator_file = input_params.mesh_restrict_BIG_file;
    }
    else if(mediator_type == "UseRestricted_B")
    {
      input_params.mediator_type = carl::MediatorType::USE_MICRO;
      input_params.mesh_mediator_file = input_params.mesh_restrict_micro_file;
    }
  }
  if (field_parser.search(1,"Dynamic")) {
    input_params.dynamic_analysis = true;
  } else {
    input_params.dynamic_analysis = false;
  }

  if (field_parser.search(1, "deltatA")) {
    input_params.deltatA = field_parser.next(input_params.deltatA);
  } else {
    input_params.deltatA = 0.25;
  }

  if (field_parser.search(1, "betaA")) {
    input_params.betaA = field_parser.next(input_params.betaA);
  } else {
    input_params.betaA = 0.25;
  }

  if (field_parser.search(1, "gammaA")) {
    input_params.gammaA = field_parser.next(input_params.gammaA);
  } else {
    input_params.gammaA = 0.5;
  }

  if (field_parser.search(1, "deltatB")) {
    input_params.deltatB = field_parser.next(input_params.deltatB);
  } else {
    input_params.deltatB = 0.25;
  }

  if (field_parser.search(1, "betaB")) {
    input_params.betaB = field_parser.next(input_params.betaB);
  } else {
    input_params.betaB = 0.25;
  }

  if (field_parser.search(1, "gammaB")) {
    input_params.gammaB = field_parser.next(input_params.gammaB);
  } else {
    input_params.gammaB = 0.5;
  }

};

};
/* Local Variables:                                                        */
/* mode: c                                                                 */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=2 ts=2 et tw=80 smartindent :                               */
