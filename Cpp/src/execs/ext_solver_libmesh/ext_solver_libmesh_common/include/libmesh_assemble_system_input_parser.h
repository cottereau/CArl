/*
 * libmesh_assemble_system_input_parser.h
 *
 *  Created on: Apr 16, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef LIBMESH_ASSEMBLE_SYSTEM_INPUT_PARSER_H_
#define LIBMESH_ASSEMBLE_SYSTEM_INPUT_PARSER_H_

#include "common_header_ext_solver_libmesh.h"
#include "ext_solver_libmesh_enums.h"
#include "newmark_param_parser.h"

struct libmesh_assemble_input_params {
  std::string mesh_file;          ///< Path to the system mesh.
  std::string physical_params_file;   ///< Physical parameters.
  WeightFunctionSystemType system_type; ///< Indicates if the system to be assembled is a micro or a macro system (used to choose the proper weight function).

  std::string mesh_weight_file;   ///< Path to the mesh containing the weight region indices.
  std::string weight_domain_idx_file; ///< Path to the file identifying the weight function regions.

  std::string output_base;  ///< Output filename base.
  bool bCalculateRBVectors; ///< Build and export the rigid body modes vectors?
  bool dynamic_analysis; // dynamic analysis
//};
//
//struct dynamic_params 
//{
  /*dynamic params*/
  carl::NewmarkParams Newmark;
  // double deltat; //step time of calculation
  // double alpha;
  // double beta; 
  // double gamma;

  double CM;            ///< Coefficient for proportional damping, *Defalut*: 0
  double CK;            ///< Coefficient for proportional damping, *Defalut*: 0
  std::string dis_vec_name, vel_vec_name, acc_vec_name;
  //unsigned int n_timesteps; // number of time step
  //unsigned int write_interval;
  //bool solver_quiet;
  //double relative_step_tolerance;
  //double relative_residual_tolerance;
  //unsigned int max_nonlinear_iterations;
  //unsigned int max_linear_iterations;
  //double initial_linear_tolerance;
  //double absolute_residual_tolerance;
};

struct material_params
{
  /* material params */
  double rho = 10000;
  double poisson_ratio = 0.5;
  double young_modulus = 1e6;
};

/** \brief Parser function for the coupled solver test programs.
 *  
 *  Required parameters:
 *    - `Mesh`: path to the mesh.
 *    - `PhysicalParameters` : physical parameters.
 *    - `SystemType` : parameter used to tell the assembler which weight functions must be used. *Values*: `Micro` or `Macro`.
 *    - `MeshWeight` : path to the mesh defining the domains of the Arlequin weight parameters.
 *    - `WeightIndexes` : path to the indices of the domains of the Arlequin weight parameters.
 *    - `NewmarkParameters` : Newmark Parameter File
 *
 *  Optional parameter:
 *    - `OutputBase` or `--output` : base of the output files (including folders). *Default*: `test_system`.
 *
 *  Boolean flags:
 *    - `ExportRBVectors` : build and export the rigid body modes vectors.
 * 
 *  * Optional parameters*:
 *    - `CM` : Coefficient for proportional damping, *Defalut*: 0
 *    - `CK` : Coefficient for proportional damping, *Defalut*: 0
 */
void get_input_params(GetPot& field_parser,
  libmesh_assemble_input_params& input_params);

//void get_input_params(GetPot& field_parser,
//    dynamic_params& input_params);


void get_input_params(GetPot& field_parser,
    material_params& input_params);

void get_file_params(GetPot command_line,
          libmesh_assemble_input_params& carl_params,  
          material_params& input_material_params);


#endif /* LIBMESH_ASSEMBLE_SYSTEM_INPUT_PARSER_H_ */
/* Local Variables:                                                        */
/* mode: c                                                                 */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=2 ts=2 et tw=80 smartindent :                               */