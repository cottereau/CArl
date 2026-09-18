/*
 * carl_time_mono_iterate_finish_input_parser.h
 *
 *  Created on: August 4, 2026
 *      Author: Romain Ruyssen
 */

#ifndef CARL_TIME_MONO_ITERATE_FINISH_INPUT_PARSER_H_
#define CARL_TIME_MONO_ITERATE_FINISH_INPUT_PARSER_H_

#include "carl_headers.h"

namespace carl
{
/// Structure containing the parameters for the setup initialization of the Time Mono solver.
struct time_mono_iterate_finish_params {
	
	// Cluster 
	ClusterSchedulerType scheduler; ///< Cluster scheduler software type. *Values*: PBS, SLURM (code not implemented for the later yet).

	// Path to "scratch" folder
	std::string scratch_folder_path;	///< Path to the folder which will be used to save the temporary files during the solve operation

	// Path for the copying function of the FETI outputs
	std::string feti_solution_path;			///< FETI solution path
    std::string results_folder_path;			///< Results folder path
	std::string results_file_name;			///< Results file name

	// Path to the time mesh file for the time loop completion check
	std::string time_mesh_file;	                ///< Path to the file containing the time mesh.
	
};

/**	\brief Parser function for the Time_Mono solver.

	Required parameters:
    - `ScratchFolderPath` : path to the folder where the temporary files used by the Time_Mono solver will be saved.
    - `TimeMeshFile` : path to the file containing the time mesh.
	- 'ResultsFolderPath' : path to the folder where the results will be saved.
	- `ResultsFileName` : name of the results files. The files will be saved in the `results_folder_path` folder.
	- `ResultsFileType` : type of the results files.
 */
void get_input_params(GetPot& field_parser,
		time_mono_iterate_finish_params& input_params);

};
#endif /* CARL_TIME_MONO_ITERATE_FINISH_INPUT_PARSER_H_ */
