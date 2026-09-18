/*
 * libmesh_init_dyn_input_parser.h
 *
 *  Created on: Apr 3, 2026
 *      Author: Romain Ruyssen
 */

#ifndef LIBMESH_INIT_SYSTEM_DYN_INPUT_PARSER_H_
#define LIBMESH_INIT_SYSTEM_DYN_INPUT_PARSER_H_

#include "common_header_ext_solver_libmesh.h"
#include "ext_solver_libmesh_enums.h"

struct libmesh_init_dyn_input_params {
	std::string mesh_file;				  ///< Path to the system mesh.
	std::string physical_params_file;	  ///< Physical parameters.
	std::string time_discretization_file; ///< Time discretization parameters.
	std::string newmark_params_file;      ///< Newmark parameters.
	WeightFunctionSystemType system_type; ///< Indicates if the system to be assembled is a micro or a macro system (used to choose the proper weight function).

	std::string mesh_weight_file;		///< Path to the mesh containing the weight region indices.
	std::string weight_domain_idx_file; ///< Path to the file identifying the weight function regions.

	std::string output_base_matrix; 	///< Output filename base.
	std::string output_base_init_cond; 	///< Output filename base.
	bool bCalculateRBVectors;	///< Build and export the rigid body modes vectors?
};

/**	\brief Parser function for the coupled solver test programs.
 *	
 *	Required parameters:
 *	  - `Mesh` : path to the mesh.
 *    - `PhysicalParameters` : physical parameters.
 *    - 'NewmarkParemeters' : Newmark scheme parameters.
 *    - `SystemType` : parameter used to tell the assembler which weight functions must be used. *Values*: `Micro` or `Macro`.
 *	  - `MeshWeight` : path to the mesh defining the domains of the Arlequin weight parameters.
 *    - `WeightIndexes` : path to the indices of the domains of the Arlequin weight parameters.
 *
 *  Optional parameter:
 *    - `OutputBase` or `--output` : base of the output files (including folders). *Default*: `test_system`.
 *
 *  Boolean flags:
 *    - `ExportRBVectors` : build and export the rigid body modes vectors.
 */
void get_input_params(GetPot& field_parser,
		libmesh_init_dyn_input_params& input_params);
#endif /* LIBMESH_INIT_SYSTEM_DYN_INPUT_PARSER_H_ */
