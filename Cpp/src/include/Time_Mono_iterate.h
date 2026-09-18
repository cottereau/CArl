/*
 * Time_Mono_iterate.h
 *
 *  Created on: August 4, 2017
 *      Author: Romain Ruyssen
 */

#ifndef TIME_MONO_ITERATE_H_
#define TIME_MONO_ITERATE_H_

#include "carl_headers.h"
#include "carl_time_mono_iterate_finish_input_parser.h"

namespace carl
{

/**	\brief Class containing the operations needed for the FETI solver.
 *
 *	This class is used by the several `CArl_Time_Mono_iterate` 
 */
class	Time_Mono_Iterate
{
protected:

	//  --- Miscellaneous declarations
	libMesh::Parallel::Communicator& m_comm;	///< Communicator
	
	// --- Params
	time_mono_iterate_finish_params m_input_params;	///< Input parameters

	// --- Lecture de l'état de l'itération précédente
	bool 	    m_bStateReadUpdated;			///< Have the state of the previous iteration been read?
	int         m_ttt;							///< Current iteration index

    // --- Copie des champs sortis de FETI dans le dossier résultats
    bool        m_bCopyFETIoutputs;			    ///< Have the FETI outputs been copied to the results folder?

	// --- Time loop completion check
	bool m_bTimeLoopChecked;					///< Has the time loop been completeness been checked?
	
	/// Default constructor
	Time_Mono_Iterate();

public:
	/// Constructor with scratch folder path, coupling matrices base filename, and libMesh communicator
	Time_Mono_Iterate(libMesh::Parallel::Communicator& comm, 
		const time_mono_iterate_finish_params& input_params) :
		m_comm { comm },
		m_input_params { input_params },
		m_bStateReadUpdated { false },
		m_ttt { 0 },
		m_bCopyFETIoutputs { false },
		m_bTimeLoopChecked { false }
	{
	};

	/// Method to read the state of the previous iteration.
	void read_and_update_state();

	/// Method to copy the FETI outputs to the results folder.
	void copy_FETI_outputs();
	
	/// Method to check if the time loop has been completed.
	bool check_time_loop_completed();

};
}

#endif /* TIME_MONO_ITERATE_H_ */