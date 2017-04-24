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

struct libmesh_assemble_input_params {
	std::string mesh_file;				  ///< Path to the system mesh.
	std::string physical_params_file;	  ///< Physical parameters.
	WeightFunctionSystemType system_type; ///< Indicates if the system to be assembled is a micro or a macro system (used to choose the proper weight function).

	std::string mesh_weight_file;		///< Path to the mesh containing the weight region indices.
	std::string weight_domain_idx_file; ///< Path to the file identifying the weight function regions.

	std::string output_base; 	///< Output filename base.
	bool bCalculateRBVectors;	///< Build and export the rigid body modes vectors?
};

/**	\brief Parser function for the coupled solver test programs.
 *	
 *	Required parameters:
 *	  - `Mesh` : path to the mesh.
 *    - `PhysicalParameters` : physical parameters.
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
		libmesh_assemble_input_params& input_params);
#endif /* LIBMESH_ASSEMBLE_SYSTEM_INPUT_PARSER_H_ */
