#include "Time_Mono_iterate.h"

namespace carl
{
    
void Time_Mono_Iterate::read_and_update_state()
{

	// ONLY read in proc 0!
	if(m_comm.rank() == 0)
	{

		// Import the scalar data
		// File format:  ttt
		std::ifstream scalar_data_if;
		scalar_data_if.open(m_input_params.scratch_folder_path + "/Time_Mono_state.dat");
		scalar_data_if >> m_ttt;
		scalar_data_if.close();

		// Incrémentation
		m_ttt += 1;

		// Update the state
		std::ofstream scalar_data_of;
		scalar_data_of.open(m_input_params.scratch_folder_path + "/Time_Mono_state.dat");
		scalar_data_of.precision(15);
		scalar_data_of << m_ttt;
		scalar_data_of << std::endl;
		scalar_data_of.close();

	}
	
	m_bStateReadUpdated = true;
};

void Time_Mono_Iterate::copy_FETI_outputs()
{

	// ONLY read in proc 0!
    if(m_comm.rank() == 0)
	{

		// Copy the FETI outputs
		std::string command_string_A;
		std::string command_string_B;
		
		std::string name_file_in_A = m_input_params.feti_solution_path + "/coupled_sol_A.petscvec";
		std::string name_file_in_B = m_input_params.feti_solution_path + "/coupled_sol_B.petscvec";
		
		std::string name_file_out_result_A = m_input_params.results_folder_path + "/" + m_input_params.results_file_name + "A_" + std::to_string(m_ttt) + ".petscvec";
		std::string name_file_out_result_B = m_input_params.results_folder_path + "/" + m_input_params.results_file_name + "B_" + std::to_string(m_ttt) + ".petscvec";
		
		command_string_A = "cp " + name_file_in_A + " " + name_file_out_result_A;
		command_string_B = "cp " + name_file_in_B + " " + name_file_out_result_B;

		carl::exec_command(command_string_A.c_str());
		carl::exec_command(command_string_B.c_str());

		std::cout << command_string_A << std::endl;
		std::cout << command_string_B << std::endl;

	}
	
	m_bCopyFETIoutputs = true;
};

bool Time_Mono_Iterate::check_time_loop_completed()
{
	homemade_assert_msg(m_bStateReadUpdated,"State not read and updated yet!");

	bool TimeLoopCompleted = false;

	// Récupération du nombre de pas de temps total
  	std::ifstream file_time_steps(m_input_params.time_mesh_file);
  	double dummy_1;
  	double dummy_2;
  	double nb_of_time_steps;
  	file_time_steps >> dummy_1 >> dummy_2 >> nb_of_time_steps;

	// Check if the time loop is completed
	if(m_ttt > nb_of_time_steps)
    {
		TimeLoopCompleted = true;
	} 

	m_bTimeLoopChecked = true;

	return TimeLoopCompleted;

};

};