/*
 * \file FETI_dyn_operations_setup.cpp
 *
 *  Created on: Nov 23,2021
 *      Author: Chensheng Luo
 * 
 * \brief **DYN**   functions responsible for for setup calculations in Solving steps
 */

#include "FETI_dyn_operations.h"

namespace carl
{

	void FETI_Dyn_Operations::init_prepare_force(int prepare_mode, 
      	std::string& input_file, 
      	int innerTimes, 
      	double small_deltat,
      	carl::DynSystemVectorPath* vectorsA,
      	carl::DynSystemVectorPath* vectorsB){
  		// File parser
  		GetPot field_parser;

  		carl::dyn_force_params input_params;
  		// If there is an input file, parse it to get the parameters. Else, parse the command line
        field_parser.parse_input_file(input_file, "#", "\n", " \t\n");
  	    get_input_params(field_parser, prepare_mode,input_params);

		if(prepare_mode == ForcePrepareMethod::MODAL_SINUS){

				this->prepare_force_vector_by_modal_and_sinus(input_params.modal_A,
					vectorsA,
					input_params.amplitude_A,
					input_params.frequency_A,
					input_params.initialPhase_A,
					small_deltat,
					0,
					innerTimes);
				this->prepare_force_vector_by_modal_and_sinus(input_params.modal_B,
					vectorsB,
					input_params.amplitude_B,
					input_params.frequency_B,
					input_params.initialPhase_B,
					small_deltat,
					0,1);
			}
		else if(prepare_mode == ForcePrepareMethod::MODAL_CONSTANT){
				this->prepare_force_vector_by_modal_and_constant(input_params.modal_A,
					vectorsA,
					input_params.amplitude_A);
				this->prepare_force_vector_by_modal_and_constant(input_params.modal_B,
					vectorsB,
					input_params.amplitude_B);
			}
		else if(prepare_mode == ForcePrepareMethod::MODAL_LINEAR){
				this->prepare_force_vector_by_modal_and_slope(input_params.modal_A,
					vectorsA,
					input_params.slope_A,
					input_params.saturation_A,
					input_params.offset_A,
					small_deltat,
					0,
					innerTimes);
				this->prepare_force_vector_by_modal_and_slope(input_params.modal_B,
					vectorsB,
					input_params.slope_B,
					input_params.saturation_B,
					input_params.offset_B,
					small_deltat,
					0,1);
			}
		else if(prepare_mode == ForcePrepareMethod::MODAL_PRODUCT){
				this->prepare_force_vector_by_modal_and_product(input_params.modal_A,
					vectorsA,
					input_params.time_series_A,
					input_params.interpolation_method_A,
					small_deltat,
					0,
					innerTimes);
				this->prepare_force_vector_by_modal_and_product(input_params.modal_B,
					vectorsB,
					input_params.time_series_B,
					input_params.interpolation_method_B,
					small_deltat,
					0,
					1);
		}else{
			homemade_error_msg("Please choose a correct force preparation mode!");
		}

		std::cout << "Force init preparation finish!" <<std::endl;

	}

	void FETI_Dyn_Operations::prepare_force_vector_by_modal_and_constant(std::string force_path,
    	carl::DynSystemVectorPath* vectors,
    	double amplitude){

		Vec modal;

		VecCreate(m_comm.get(),&modal);
		carl::read_PETSC_vector(modal,force_path,m_comm.get());

		VecScale(modal,amplitude);


		carl::write_PETSC_vector(modal,vectors->this_force,0,m_comm.get(),1);
		carl::write_PETSC_vector(modal,vectors->next_force,0,m_comm.get(),1);
		#ifdef PRINT_MATLAB_DEBUG
			libMesh::PetscVector<libMesh::Number> vector_force(modal,m_comm);
   			vector_force.print_matlab(vectors->this_force+".m");
   			vector_force.print_matlab(vectors->next_force+".m");
		#endif

		VecDestroy(&modal);

	}

	void FETI_Dyn_Operations::prepare_force_vector_by_modal_and_sinus(std::string force_path,
    	carl::DynSystemVectorPath* vectors,
    	double amplitude,
    	double frequency,
    	double initialPhase,
    	double small_deltat,
    	int index,
    	int timestep)
	{
		Vec modal;
		Vec forceThis,forceNext;

		VecCreate(m_comm.get(),&modal);
		carl::read_PETSC_vector(modal,force_path,m_comm.get());
		VecDuplicate(modal,&forceThis);
		VecDuplicate(modal,&forceNext);


		VecCopy(modal,forceThis);
		VecScale(forceThis,amplitude*sin(2*M_PI*frequency*small_deltat*index+initialPhase));
		carl::write_PETSC_vector(forceThis,vectors->this_force,0,m_comm.get(),1);

		VecCopy(modal,forceNext);
		VecScale(forceNext,amplitude*sin(2*M_PI*frequency*small_deltat*(index+timestep)+initialPhase));
		carl::write_PETSC_vector(forceNext,vectors->next_force,0,m_comm.get(),1);


		#ifdef PRINT_MATLAB_DEBUG
			libMesh::PetscVector<libMesh::Number> vector_force_this(forceThis,m_comm);
   			vector_force_this.print_matlab(vectors->this_force+".m");
			libMesh::PetscVector<libMesh::Number> vector_force_next(forceNext,m_comm);
   			vector_force_next.print_matlab(vectors->next_force+".m");
		#endif

		VecDestroy(&modal);
		VecDestroy(&forceThis);
		VecDestroy(&forceNext);
	}

	void FETI_Dyn_Operations::prepare_force_vector_by_modal_and_slope(std::string force_path,
    	carl::DynSystemVectorPath* vectors,
    	double slope,
    	double saturation,
    	double offset,
    	double small_deltat,
    	int index,
    	int timestep){

		Vec modal;
		Vec forceThis,forceNext;

		VecCreate(m_comm.get(),&modal);
		carl::read_PETSC_vector(modal,force_path,m_comm.get());
		VecDuplicate(modal,&forceThis);
		VecDuplicate(modal,&forceNext);


		VecCopy(modal,forceThis);
		VecScale(forceThis,std::min(slope*small_deltat*index+offset, saturation));
		carl::write_PETSC_vector(forceThis,vectors->this_force,0,m_comm.get(),1);

		VecCopy(modal,forceNext);
		VecScale(forceNext,std::min(slope*small_deltat*(index+timestep)+offset,saturation));
		carl::write_PETSC_vector(forceNext,vectors->next_force,0,m_comm.get(),1);


		#ifdef PRINT_MATLAB_DEBUG
			libMesh::PetscVector<libMesh::Number> vector_force_this(forceThis,m_comm);
   			vector_force_this.print_matlab(vectors->this_force+".m");
			libMesh::PetscVector<libMesh::Number> vector_force_next(forceNext,m_comm);
   			vector_force_next.print_matlab(vectors->next_force+".m");
		#endif

		VecDestroy(&modal);
		VecDestroy(&forceThis);
		VecDestroy(&forceNext);

	}

	void FETI_Dyn_Operations::prepare_force_vector_by_modal_and_product(std::string force_path,
    	carl::DynSystemVectorPath* vectors,
    	std::string time_series_file,
    	int intepolation_method,
    	double small_deltat,
    	int index,
    	int timestep){

		Vec modal;
		Vec forceThis,forceNext;

		VecCreate(m_comm.get(),&modal);
		carl::read_PETSC_vector(modal,force_path,m_comm.get());
		VecDuplicate(modal,&forceThis);
		VecDuplicate(modal,&forceNext);

		// Read time series file
		std::ifstream inFile(time_series_file);
		std::vector<double> times,values;
		std::string line;
		double time, value;

		if(inFile){
			while(getline(inFile, line, '\n')){
				std::istringstream ss(line);

				ss >> time >> value;
				times.push_back(time);
            	values.push_back(value);
            }
		}else{
			homemade_error_msg("[CArl Error] Time series file can't read!");
		}

		double this_scaling_factor, next_scaling_factor;
		

		if(intepolation_method == carl::ForceInterpolationMethod::OFF){
			double epsilon = 1e-3 * small_deltat;// Tolerance level
			std::vector<double>::iterator it_this, it_next;

			std::cout << "Force preparation: Get scaling factor at "<< small_deltat*index << "without interpolation ... ";

			it_this = std::find_if(times.begin(),times.end(),[small_deltat,index,epsilon](double b) { return abs( small_deltat*index - b) < epsilon; });
			if (it_this != times.end()){
				this_scaling_factor = values.at(std::distance(times.begin(), it_this));
				std::cout << "done!" <<std::endl;
			}else{
				homemade_error_msg("[CArl Error] Time"+std::to_string(small_deltat*index)+"can't be find in the time series file.");
			}

			std::cout << "Force preparation: Get scaling factor at "<< small_deltat*(index+timestep) << "without interpolation ...";
			it_next = std::find_if(times.begin(),times.end(),[small_deltat,index,epsilon,timestep](double b) { return abs( small_deltat*(index+timestep) - b) < epsilon; });
			if (it_next != times.end()){
				next_scaling_factor = values.at(std::distance(times.begin(), it_next));
				std::cout << "done!" <<std::endl;
			}else{
				homemade_error_msg("[CArl Error] Time"+std::to_string(small_deltat*(index+timestep))+"can't be find in the time series file.");
			}
		}else if(intepolation_method == carl::ForceInterpolationMethod::NEAREST){

			std::cout << "Force preparation: Get scaling factor at "<< small_deltat*index << "with nearest interpolation ... ";

			for(auto& element : times){
    			element = abs(element - small_deltat*index);
			}
			std::cout << "done!" <<std::endl;
			
			this_scaling_factor = values.at(std::distance(std::begin(times), std::min_element(std::begin(times), std::end(times))));
			
			std::cout << "Force preparation: Get scaling factor at "<< small_deltat*(index+timestep) << "with nearest interpolation ... ";

			for(auto& element : times){
    			element = abs(element - std::copysign(small_deltat*timestep,element));
			}
			
			next_scaling_factor = values.at(std::distance(std::begin(times), std::min_element(std::begin(times), std::end(times))));
			std::cout << "done!" <<std::endl;

		}else if(intepolation_method == carl::ForceInterpolationMethod::LINEAR){
			double before_time, back_time;
			double before_value, back_value;

			int gotten_index, total_index;

			std::cout << "Force preparation: Get scaling factor at "<< small_deltat*index << "with linear interpolation ... ";

			for(auto& element : times){
    			element = abs(element - small_deltat*index);
			}
			
			gotten_index = std::distance(std::begin(times), std::min_element(std::begin(times), std::end(times)));
			total_index = std::distance(std::begin(times), std::end(times));
			
			before_value = values.at(gotten_index);
			if((before_value <= 0 && gotten_index < total_index) || gotten_index == 0){
				before_time = times.at(gotten_index);
				back_time = times.at(gotten_index + 1 );
				back_value = values.at(gotten_index +1 );
			}else{
				back_value = before_value;
				back_time = times.at(gotten_index);
				before_value = values.at(gotten_index - 1);
				before_time = times.at(gotten_index - 1);
			}

			this_scaling_factor = (back_value-before_value)/(back_time-before_time)*(time-before_time)+before_value;
			std::cout << "done!" <<std::endl;

			std::cout << "Force preparation: Get scaling factor at "<< small_deltat*(index+timestep) << "with linear interpolation ... ";
			for(auto& element : times){
    			element = abs(element - std::copysign(small_deltat*timestep,element));
			}

			gotten_index = std::distance(std::begin(times), std::min_element(std::begin(times), std::end(times)));
			
			before_value = values.at(gotten_index);
			if((before_value <= 0 && gotten_index < total_index) || gotten_index == 0){
				before_time = times.at(gotten_index);
				back_time = times.at(gotten_index + 1 );
				back_value = values.at(gotten_index +1 );
			}else{
				back_value = before_value;
				back_time = times.at(gotten_index);
				before_value = values.at(gotten_index - 1);
				before_time = times.at(gotten_index - 1);
			}

			next_scaling_factor = (back_value-before_value)/(back_time-before_time)*(time-before_time)+before_value;
			std::cout << "done!" <<std::endl;
		}

		// This force at small_deltat*index
		VecCopy(modal,forceThis);
		VecScale(forceThis,this_scaling_factor);
		carl::write_PETSC_vector(forceThis,vectors->this_force,0,m_comm.get(),1);

		// Next force at small_deltat*(index+timestep)
		VecCopy(modal,forceNext);
		VecScale(forceNext,next_scaling_factor);
		carl::write_PETSC_vector(forceNext,vectors->next_force,0,m_comm.get(),1);


		#ifdef PRINT_MATLAB_DEBUG
			libMesh::PetscVector<libMesh::Number> vector_force_this(forceThis,m_comm);
   			vector_force_this.print_matlab(vectors->this_force+".m");
			libMesh::PetscVector<libMesh::Number> vector_force_next(forceNext,m_comm);
   			vector_force_next.print_matlab(vectors->next_force+".m");
		#endif

		VecDestroy(&modal);
		VecDestroy(&forceThis);
		VecDestroy(&forceNext);

	}

	    //TODO: ADD ACCELERATION INITIAL CONDITION!
	void FETI_Dyn_Operations::set_initial_condition(carl::DynInitialVectorPath* initial_A,
      carl::DynInitialVectorPath* initial_B,
      carl::DynSystemMatrixPath* matrix_A,
      carl::DynSystemMatrixPath* matrix_B,
      carl::NewmarkParams* newmark_A,
      carl::NewmarkParams* newmark_B,
      carl::DynSystemVectorPath* vectors_A,
      carl::DynSystemVectorPath* vectors_B,
      std::string result_folder_path)
	{
		Vec dispA,dispB;
	  	Vec speedA,speedB;
	  	Vec forceA,forceB;
	  	Vec rhsA,rhsB;
	  	Mat stiffnessA,stiffnessB;
	  	Mat dampingA,dampingB;
	  	PetscInt MA,NA,MB,NB;

	  	VecCreate(m_comm.get(),&dispA);
		VecCreate(m_comm.get(),&dispB);
		VecCreate(m_comm.get(),&speedA);
		VecCreate(m_comm.get(),&speedB);
		VecCreate(m_comm.get(),&forceA);
		VecCreate(m_comm.get(),&forceB);
		VecCreate(m_comm.get(),&rhsA);
		VecCreate(m_comm.get(),&rhsB);
		MatCreate(m_comm.get(),&stiffnessA);
		MatCreate(m_comm.get(),&stiffnessB);
		MatCreate(m_comm.get(),&dampingA);
		MatCreate(m_comm.get(),&dampingB);

		carl::read_PETSC_matrix(stiffnessA,matrix_A->stiffness,m_comm.get());
		carl::read_PETSC_matrix(stiffnessB,matrix_B->stiffness,m_comm.get());
		carl::read_PETSC_matrix(dampingA,matrix_A->damping,m_comm.get());
		carl::read_PETSC_matrix(dampingB,matrix_B->damping,m_comm.get());

		MatGetSize(stiffnessA,&MA,&NA);
		MatGetSize(stiffnessB,&MB,&NB);

		if(1-(initial_A->disp).compare("")){
			//No imposed initial condition
			VecSetSizes(dispA,PETSC_DECIDE,NA);
  			VecSetFromOptions(dispA);
  			VecZeroEntries(dispA);
		}else{
			carl::read_PETSC_vector(dispA,initial_A->disp,m_comm.get());
		}

		if(1-(initial_B->disp).compare("")){
			//No imposed initial condition
			VecSetSizes(dispB,PETSC_DECIDE,NB);
  			VecSetFromOptions(dispB);
  			VecZeroEntries(dispB);
		}else{
			carl::read_PETSC_vector(dispB,initial_B->disp,m_comm.get());
		}

		if(1-(initial_A->speed).compare("")){
			//No imposed initial condition
			VecSetSizes(speedA,PETSC_DECIDE,NA);
  			VecSetFromOptions(speedA);
  			VecZeroEntries(speedA);
		}else{
			carl::read_PETSC_vector(speedA,initial_A->speed,m_comm.get());
		}

		if(1-(initial_B->speed).compare("")){
			//No imposed initial condition
			VecSetSizes(speedB,PETSC_DECIDE,NB);
  			VecSetFromOptions(speedB);
  			VecZeroEntries(speedB);
		}else{
			carl::read_PETSC_vector(speedB,initial_B->speed,m_comm.get());
		}


		carl::write_PETSC_vector(dispA,result_folder_path+"/disp_A_"+std::to_string(0)+".petscvec",0,m_comm.get(),1);
		carl::write_PETSC_vector(dispB,result_folder_path+"/disp_B_"+std::to_string(0)+".petscvec",0,m_comm.get(),1);
		#ifdef OUTPUT_SPEED_DEBUG
			carl::write_PETSC_vector(speedA,result_folder_path+"/speed_A_"+std::to_string(0)+".petscvec",0,m_comm.get(),1);
			carl::write_PETSC_vector(speedB,result_folder_path+"/speed_B_"+std::to_string(0)+".petscvec",0,m_comm.get(),1);
		#endif 
		carl::write_PETSC_vector(dispA,vector_A.prev_disp,0,m_comm.get(),1);
		carl::write_PETSC_vector(dispB,vector_B.prev_disp,0,m_comm.get(),1);
		carl::write_PETSC_vector(speedA,vector_A.prev_speed,0,m_comm.get(),1);
		carl::write_PETSC_vector(speedB,vector_B.prev_speed,0,m_comm.get(),1);


		VecSetSizes(rhsA,PETSC_DECIDE,MA);
  		VecSetFromOptions(rhsA);

  		VecSetSizes(rhsB,PETSC_DECIDE,MB);
  		VecSetFromOptions(rhsB);

  		//RHS=
		carl::read_PETSC_vector(forceA,vectors_A->next_force,m_comm.get());
		VecAXPY(dispA,newmark_A->deltat,speedA);
		//TODO:ADD IF HAVE ACCELERATION
		MatMult(stiffnessA,dispA,rhsA);
		//TODO:ADD IF HAVE ACCELERATION
		MatMultAdd(dampingA,speedA,rhsA,rhsA);
		VecScale(rhsA,-1);
		VecAXPY(rhsA,1,forceA);

		carl::read_PETSC_vector(forceB,vectors_B->next_force,m_comm.get());
		VecAXPY(dispB,newmark_B->deltat,speedB);
		//TODO:ADD IF HAVE ACCELERATION
		MatMult(stiffnessB,dispB,rhsB);
		//TODO:ADD IF HAVE ACCELERATION
		MatMultAdd(dampingB,speedB,rhsB,rhsB);
		VecScale(rhsB,-1);
		VecAXPY(rhsB,1,forceB);

		carl::write_PETSC_vector(rhsA,vector_A.rhs_free,0,m_comm.get(),1);
		carl::write_PETSC_vector(rhsB,vector_B.rhs_free,0,m_comm.get(),1);

		VecDestroy(&dispA);
		VecDestroy(&dispB);
	  	VecDestroy(&speedA);
	  	VecDestroy(&speedB);
	  	VecDestroy(&forceA);
	  	VecDestroy(&forceB);
	  	VecDestroy(&rhsA);
	  	VecDestroy(&rhsB);
	  	MatDestroy(&stiffnessA);
	  	MatDestroy(&stiffnessB);
	  	MatDestroy(&dampingA);
	  	MatDestroy(&dampingB);
	}
	  



}