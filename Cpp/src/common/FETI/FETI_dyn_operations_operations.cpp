/*
 * \file FETI_dyn_operations_setup.cpp
 *
 *  Created on: Nov 23,2021
 *      Author: Chensheng Luo
 * 
 * \brief **DYN-DI/DYN-CG**   functions responsible for for setup calculations in Solving steps
 */

#include "FETI_dyn_operations.h"

namespace carl
{

	void FETI_Dyn_Operations::rhs_free(carl::NewmarkParams* newmark_params,
      carl::DynSystemVectorPath* vectors,
      carl::DynSystemMatrixPath* matrices){


		homemade_assert_msg(m_bForceAPrepared | m_bForceBPrepared ,"Force not prepared yet!");


		//Used only after move_to_prev

			Vec force;
			Vec displacement;
      Vec speed;
			Vec acceleration;
			Mat stiffness;
			Mat damping;
			Vec rhs;
			PetscInt M,N,Nvec;

			//Initialize vec and mat
			VecCreate(m_comm.get(),&force);
			VecCreate(m_comm.get(),&displacement);
			VecCreate(m_comm.get(),&speed);
			VecCreate(m_comm.get(),&acceleration);
			MatCreate(m_comm.get(),&stiffness);
			MatCreate(m_comm.get(),&damping);
			VecCreate(m_comm.get(),&rhs);

			// Get vec and mat from file
			carl::read_PETSC_vector(displacement,vectors->prev_disp,m_comm.get());
			carl::read_PETSC_vector(speed,vectors->prev_speed,m_comm.get());
			carl::read_PETSC_vector(acceleration,vectors->prev_acc,m_comm.get());
			carl::read_PETSC_matrix(stiffness,matrices->stiffness,m_comm.get());
			carl::read_PETSC_matrix(damping,matrices->damping,m_comm.get());

			//Set the dimension for vec
			MatGetSize(stiffness,&M,&N);
  		VecGetSize(displacement,&Nvec);
  		if(Nvec!=N){
  			homemade_error_msg("Matrix dimension doesn't match!");
  		}
  		VecSetSizes(rhs,PETSC_DECIDE,M);
  		VecSetFromOptions(rhs);

			//disp_predictor=displacement+deltat*speed+(1/2-beta)*deltat^2*acceleration
			VecAXPY(displacement,newmark_params->deltat,speed);
			VecAXPY(displacement,(0.5-newmark_params->beta)*newmark_params->deltat*newmark_params->deltat,acceleration);
			//rhs=(1-alpha)*F+alpha*Fprev-(1-alpha)*K*disp_predictor-alpha*K*disp_prev-D*speed_predictor
			MatMult(stiffness,displacement,rhs);
			VecScale(rhs,(1-newmark_params->alpha));
			//speed_predictor=speed+(1-gamma)*deltat*acceleration
			VecAXPY(speed,(1-newmark_params->gamma)*newmark_params->deltat,acceleration);
			MatMultAdd(damping,speed,rhs,rhs);
			VecScale(rhs,-1);

			//This force
			std::FILE *fh ;
			fh = fopen((vectors->next_force).c_str (),"r");// To examine if there is a force input file

			if(fh!=NULL){
				//There is a force input file
				carl::read_PETSC_vector(force,vectors->next_force,m_comm.get());
				VecAXPY(rhs,1-newmark_params->alpha,force);
			}

			//Prev force
			fh = fopen((vectors->this_force).c_str (),"r");// To examine if there is a force input file

			if(fh!=NULL){
				//There is a force input file
				carl::read_PETSC_vector(force,vectors->this_force,m_comm.get());
				VecAXPY(rhs,newmark_params->alpha,force);
			}

			carl::write_PETSC_vector(rhs,vectors->rhs_free,0,m_comm.get(),1);

			#ifdef PRINT_MATLAB_DEBUG
				libMesh::PetscVector<libMesh::Number> vector_lib(rhs,m_comm);
   			vector_lib.print_matlab(vectors->rhs_free+".m");
			#endif

			VecDestroy(&force);
			VecDestroy(&displacement);
			VecDestroy(&speed);
			VecDestroy(&acceleration);
			MatDestroy(&stiffness);
			MatDestroy(&damping);
			VecDestroy(&rhs);
			
	}
	void FETI_Dyn_Operations::Newmark_speed_free(carl::NewmarkParams* newmark_params,
      carl::DynSystemVectorPath* vectors){

      	// Get acceleration
			Vec speed;
			Vec prev_acceleration;
			Vec this_acceleration_free;
			VecCreate(m_comm.get(),&speed);
			VecCreate(m_comm.get(),&prev_acceleration);
			VecCreate(m_comm.get(),&this_acceleration_free);
			carl::read_PETSC_vector(this_acceleration_free,vectors->this_acc_free,m_comm.get());


		std::FILE *fh ;
		fh = fopen((vectors->prev_acc).c_str (),"r");// To examine if there is a "previous" file
		
		if(fh!=NULL){
			//Not the first moment
			carl::read_PETSC_vector(prev_acceleration,vectors->prev_acc,m_comm.get());
			carl::read_PETSC_vector(speed,vectors->prev_speed,m_comm.get());
			
			//predictor=speed+(1-gamma)*deltat*prev_acceleration
			VecAXPY(speed,(1-newmark_params->gamma)*newmark_params->deltat,prev_acceleration);

			//new_speed=predictor+gamma*deltat*this_acceleration_free
			VecAXPY(speed,newmark_params->gamma*newmark_params->deltat,this_acceleration_free);
		}else{
			//No prev_speed vector exists! We need create one.
			VecDuplicate(this_acceleration_free,&speed);
			VecCopy(this_acceleration_free,speed);
			VecScale(speed,newmark_params->gamma*newmark_params->deltat);
		}

		//Calculate 
		carl::write_PETSC_vector(speed,vectors->this_speed_free,0,m_comm.get(),1);

		VecDestroy(&speed);
		VecDestroy(&prev_acceleration);
		VecDestroy(&this_acceleration_free);
	}


  void FETI_Dyn_Operations::Newmark_displacement_free(carl::NewmarkParams* newmark_params,
      carl::DynSystemVectorPath* vectors){

      Vec displacement;
      Vec prev_speed;
      Vec prev_acceleration;
			Vec this_acceleration_free;
			VecCreate(m_comm.get(),&displacement);
			VecCreate(m_comm.get(),&prev_speed);
			VecCreate(m_comm.get(),&prev_acceleration);
			VecCreate(m_comm.get(),&this_acceleration_free);
			carl::read_PETSC_vector(this_acceleration_free,vectors->this_acc_free,m_comm.get());

		std::FILE *fh;
		fh = fopen((vectors->prev_acc).c_str (),"r");// To examine if there is a "previous" file

		if(fh!=NULL){
			//Not the first moment
			carl::read_PETSC_vector(displacement,vectors->prev_disp,m_comm.get());
			carl::read_PETSC_vector(prev_speed,vectors->prev_speed,m_comm.get());
			carl::read_PETSC_vector(prev_acceleration,vectors->prev_acc,m_comm.get());
			
			//predictor=displacement+deltat*prev_speed+(1/2-beta)*deltat^2*prev_acceleration
			VecAXPY(displacement,newmark_params->deltat,prev_speed);
			VecAXPY(displacement,(0.5-newmark_params->beta)*newmark_params->deltat*newmark_params->deltat,prev_acceleration);

			//new_displacement=predictor+beta*deltat*deltat*this_acceleration_free
			VecAXPY(displacement,newmark_params->beta*newmark_params->deltat*newmark_params->deltat,this_acceleration_free);
		}else{
			//No prev_displacement vector exists! We need create one.
			VecDuplicate(this_acceleration_free,&displacement);
			VecCopy(this_acceleration_free,displacement);
			VecScale(displacement,newmark_params->beta*newmark_params->deltat*newmark_params->deltat);
		}

		//Calculate 
		carl::write_PETSC_vector(displacement,vectors->this_disp_free,0,m_comm.get(),1);

		VecDestroy(&displacement);
		VecDestroy(&prev_speed);
		VecDestroy(&prev_acceleration);
		VecDestroy(&this_acceleration_free);
    }

  void FETI_Dyn_Operations::interpolate_A_disp(int jjj,int m){
  	Vec prev_vector;
  	Vec this_vector;
  	std::string prev_vec_path=vector_A.prev_disp_free;
  	std::string this_vec_path=vector_A.this_disp_free;

  	VecCreate(m_comm.get(),&prev_vector);
		VecCreate(m_comm.get(),&this_vector);
		carl::read_PETSC_vector(this_vector,this_vec_path,m_comm.get());

		std::FILE *fh;
		fh = fopen(prev_vec_path.c_str (),"r");// To examine if this is a "previous" file

		if(fh!=NULL){
			//Not the first moment
			carl::read_PETSC_vector(prev_vector,prev_vec_path,m_comm.get());
			
			//U_free_j=(1-j/m)*prev_acc+j/m*this_acc
			VecAXPBY(this_vector,(1-(jjj)/m),(jjj)/m,prev_vector);
		}else{
			VecScale(this_vector,(jjj)/m);
		}

		carl::write_PETSC_vector(this_vector,vector_A.inter_disp_free,0,m_comm.get(),1);

		VecDestroy(&prev_vector);
		VecDestroy(&this_vector);
  }

  void FETI_Dyn_Operations::rhs_interpolation(std::string& coupling_matrix_A_path,
  	std::string& coupling_matrix_B_path){
  		Mat coupling_matrix_A;
  		Mat coupling_matrix_B;
  		Vec disp_A;
  		Vec disp_B;
  		Vec result1;
  		Vec result2;
  		PetscInt M,N,Nvec;

  		MatCreate(m_comm.get(),&coupling_matrix_A);
  		MatCreate(m_comm.get(),&coupling_matrix_B);
  		VecCreate(m_comm.get(),&disp_A);
  		VecCreate(m_comm.get(),&disp_B);
  		VecCreate(m_comm.get(),&result1);
  		VecCreate(m_comm.get(),&result2);

  		carl::read_PETSC_vector(disp_A,vector_A.inter_disp_free);
  		carl::read_PETSC_vector(disp_B,vector_B.this_disp_free);
  		carl::read_PETSC_matrix(coupling_matrix_A,coupling_matrix_A_path);
  		carl::read_PETSC_matrix(coupling_matrix_B,coupling_matrix_B_path);

  		MatGetSize(coupling_matrix_A,&M,&N);
  		VecGetSize(disp_A,&Nvec);
  		if(Nvec!=N){
  			homemade_error_msg("Matrix dimension doesn't match!");
  		}

  		VecSetSizes(result1,PETSC_DECIDE,M);
  		VecSetFromOptions(result1);
  		VecSetSizes(result2,PETSC_DECIDE,M);
  		VecSetFromOptions(result2);

  		// result1=C_A U_A
  		MatMult(coupling_matrix_A,disp_A,result1);

  		// result2=C_B U_B
  		MatMult(coupling_matrix_B,disp_B,result2);

  		// result2=result1-result2
  		VecAYPX(result2,-1,result1);
  		//VecAXPBY(result2,0,0,result1);//TEST USAGE FOR CASE WITHOUT COUPLING!
  		carl::write_PETSC_vector(result2,m_rhs_interpolation,0,m_comm.get(),1);

  		#ifdef PRINT_MATLAB_DEBUG
			libMesh::PetscVector<libMesh::Number> vector_lib(result2,m_comm);
   			vector_lib.print_matlab(m_rhs_interpolation+".m");
			#endif

  		MatDestroy(&coupling_matrix_A);
  		MatDestroy(&coupling_matrix_B);
  		VecDestroy(&disp_A);
  		VecDestroy(&disp_B);
  		VecDestroy(&result1);
  		VecDestroy(&result2);
    }

  void FETI_Dyn_Operations::rhs_link(carl::DynSystemVectorPath* vectors,
  		std::string coupling_vector_path,
      carl::DynSystemMatrixPath* matrices){
  		Mat coupling_matrix;
  		Vec interpolation_vector;
  		Vec result;
  		PetscInt M,N,Nvec;

  		MatCreate(m_comm.get(),&coupling_matrix);
  		VecCreate(m_comm.get(),&interpolation_vector);
  		VecCreate(m_comm.get(),&result);


  		carl::read_PETSC_matrix(coupling_matrix,matrices->coupling);
  		carl::read_PETSC_vector(interpolation_vector,coupling_vector_path);

  		MatGetSize(coupling_matrix,&M,&N);
  		VecGetSize(interpolation_vector,&Nvec);
  		if(Nvec!=M){
  			homemade_error_msg("Matrix dimension doesn't match between coupling matrix and coupling vector!");
  		}

  		VecSetSizes(result,PETSC_DECIDE,N);
  		VecSetFromOptions(result);

  		MatMultTranspose(coupling_matrix,interpolation_vector,result);

  		VecScale(result,vectors->coupling_sign);

  		carl::write_PETSC_vector(result,vectors->rhs_link,0,m_comm.get(),1);

  		MatDestroy(&coupling_matrix);
  		VecDestroy(&interpolation_vector);
  		VecDestroy(&result);
  }

  void FETI_Dyn_Operations::Newmark_speed_link(carl::NewmarkParams* newmark_params,
      carl::DynSystemVectorPath* vectors){
  		Vec this_acceleration_link;

  		VecCreate(m_comm.get(),&this_acceleration_link);

  		carl::read_PETSC_vector(this_acceleration_link,vectors->this_acc_link);

  		VecScale(this_acceleration_link,newmark_params->gamma*newmark_params->deltat);

  		carl::write_PETSC_vector(this_acceleration_link,vectors->this_speed_link,0,m_comm.get(),1);

  		VecDestroy(&this_acceleration_link);

  }
    void FETI_Dyn_Operations::Newmark_displacement_link(carl::NewmarkParams* newmark_params,
      carl::DynSystemVectorPath* vectors){
    	Vec this_acceleration_link;

  		VecCreate(m_comm.get(),&this_acceleration_link);

  		carl::read_PETSC_vector(this_acceleration_link,vectors->this_acc_link);

  		VecScale(this_acceleration_link,newmark_params->beta*newmark_params->deltat*newmark_params->deltat);

  		carl::write_PETSC_vector(this_acceleration_link,vectors->this_disp_link,0,m_comm.get(),1);

  		VecDestroy(&this_acceleration_link);

    }

    void FETI_Dyn_Operations::add_free_link(carl::DynSystemVectorPath* vectors){
    	this->add_free_link_one(vectors->this_acc_free,vectors->this_acc_link,vectors->this_acc);
    	this->add_free_link_one(vectors->this_disp_free,vectors->this_disp_link,vectors->this_disp);
    	this->add_free_link_one(vectors->this_speed_free,vectors->this_speed_link,vectors->this_speed);
    }

    void FETI_Dyn_Operations::add_free_link_one(std::string& free_term_path,
      std::string& link_term_path,
      std::string& output_path){
    	Vec free_term;
    	Vec link_term;

    	VecCreate(m_comm.get(),&free_term);
    	VecCreate(m_comm.get(),&link_term);

    	carl::read_PETSC_vector(free_term,free_term_path);
    	carl::read_PETSC_vector(link_term,link_term_path);

    	VecAXPY(link_term,1,free_term);

    	carl::write_PETSC_vector(link_term,output_path,0,m_comm.get(),1);

    	VecDestroy(&free_term);
    	VecDestroy(&link_term);
    }

    void FETI_Dyn_Operations::init_test_coupling(){
    	if(m_comm.rank() == 0)
  		{
  			std::ofstream result_file;
				result_file.open(m_scratch_folder_path + "/result_test_inversion.txt",std::ofstream::app);
				result_file.precision(15);
				result_file << "Round" << "   " << "Test norm "<< "   " << "Product norm"<< "   " << "Quotient" << std::endl;
				result_file.close();
			}
    }

    void FETI_Dyn_Operations::test_coupling(std::string& coupling_matrix_A_path,
    	std::string& coupling_matrix_B_path,
    	int index){
    	Mat CA,CB;
  		Vec UA,UB;
  		Vec result1,result2;
  		PetscInt M,N,Nvec;
  		PetscReal norm1;
			PetscReal norm2;

  		MatCreate(m_comm.get(),&CA);
  		VecCreate(m_comm.get(),&UA);
  		MatCreate(m_comm.get(),&CB);
  		VecCreate(m_comm.get(),&UB);
  		VecCreate(m_comm.get(),&result1);
  		VecCreate(m_comm.get(),&result2);

  		carl::read_PETSC_matrix(CA,coupling_matrix_A_path);
  		carl::read_PETSC_matrix(CB,coupling_matrix_B_path);
  		carl::read_PETSC_vector(UA,m_result_folder_path+"/disp_A_"+std::to_string(index)+".petscvec");
  		carl::read_PETSC_vector(UB,m_result_folder_path+"/disp_B_"+std::to_string(index)+".petscvec");

  		MatGetSize(CA,&M,&N);
  		VecGetSize(UA,&Nvec);
  		if(Nvec!=N){
  			homemade_error_msg("Matrix dimension doesn't match!");
  		}

  		VecSetSizes(result1,PETSC_DECIDE,M);
  		VecSetFromOptions(result1);
  		VecSetSizes(result2,PETSC_DECIDE,M);
  		VecSetFromOptions(result2);

  		// result1=C_A U_A
  		MatMult(CA,UA,result1);

  		// result2=C_B U_B
  		MatMult(CB,UB,result2);

  		// result2=result1-result2
  		VecAYPX(result2,-1,result1);

  		VecNorm(result2,NORM_1,&norm1);

  		VecNorm(result1,NORM_1,&norm2);
  		if(m_comm.rank() == 0)
  		{
  		std::ofstream result_file;
			result_file.open(m_scratch_folder_path + "/result_test_inversion.txt",std::ofstream::app);
			result_file.precision(15);
			result_file << "  " << index << "   " << norm1 << "    "<< norm2<< "   "<< norm1/norm2<< std::endl;
			result_file.close();

  		VecDestroy(&UA);
    	VecDestroy(&UB);
    	VecDestroy(&result1);
    	VecDestroy(&result2);
    	MatDestroy(&CA);
    	MatDestroy(&CB);

    	std::FILE *fh ;
    	std::string filename=m_scratch_folder_path+"/CG_solver/FETI_convergence.dat";
			fh = fopen(filename.c_str (),"r");// To examine if there is a cv file

			if(fh!=NULL){
    
    			std::string command_string;
    			command_string = "cp " + filename + " "+ m_scratch_folder_path+"/FETI_convergence_"+std::to_string(index)+".dat";
    			carl::exec_command(command_string.c_str());

  		}

    }

    }

    void FETI_Dyn_Operations::prepare_A_next_force(carl::DynSystemVectorPath* vectors,
    	int index,
    	int step,
    	int prepare_mode,
    	std::string& force_prepare_file,
    	double small_deltat){
    	// File parser
  		GetPot field_parser;

  		carl::dyn_force_params input_params;
  		// If there is an input file, parse it to get the parameters. Else, parse the command line
      field_parser.parse_input_file(force_prepare_file, "#", "\n", " \t\n");
  	  get_input_params(field_parser, prepare_mode,input_params);

		if(prepare_mode == ForcePrepareMethod::MODAL_SINUS){

				this->prepare_force_vector_by_modal_and_sinus(input_params.modal_A,
					vectors,
					input_params.amplitude_A,
					input_params.frequency_A,
					input_params.initialPhase_A,
					small_deltat,
					index,
					step);
			}
		else if(prepare_mode == ForcePrepareMethod::MODAL_CONSTANT){
				this->prepare_force_vector_by_modal_and_constant(input_params.modal_A,
					vectors,
					input_params.amplitude_A);
			}
		else if(prepare_mode == ForcePrepareMethod::MODAL_LINEAR){
				this->prepare_force_vector_by_modal_and_slope(input_params.modal_A,
					vectors,
					input_params.slope_A,
					small_deltat,
					index,
					step);
			}
		else if(prepare_mode == ForcePrepareMethod::MODAL_PRODUCT){
				homemade_error_msg("Not Implemented!");
		}else{
				homemade_error_msg("Please choose a correct force preparation mode!");
		}

		std::cout << "Force A preparation finish!" <<std::endl;

		m_bForceAPrepared = true;
    }


    void FETI_Dyn_Operations::prepare_B_next_force(carl::DynSystemVectorPath* vectors,
      int index,
      int prepare_mode,
      std::string& force_prepare_file,
      double small_deltat){



    	// File parser
  		GetPot field_parser;

  		carl::dyn_force_params input_params;
  		// If there is an input file, parse it to get the parameters. Else, parse the command line
      field_parser.parse_input_file(force_prepare_file, "#", "\n", " \t\n");
  	  get_input_params(field_parser, prepare_mode,input_params);

		if(prepare_mode == ForcePrepareMethod::MODAL_SINUS){

				this->prepare_force_vector_by_modal_and_sinus(input_params.modal_B,
					vectors,
					input_params.amplitude_B,
					input_params.frequency_B,
					input_params.initialPhase_B,
					small_deltat,
					index,
					1);
			}
		else if(prepare_mode == ForcePrepareMethod::MODAL_CONSTANT){
				this->prepare_force_vector_by_modal_and_constant(input_params.modal_B,
					vectors,
					input_params.amplitude_B);
			}
		else if(prepare_mode == ForcePrepareMethod::MODAL_LINEAR){
				this->prepare_force_vector_by_modal_and_slope(input_params.modal_B,
					vectors,
					input_params.slope_B,
					small_deltat,
					index,
					1);
			}
		else if(prepare_mode == ForcePrepareMethod::MODAL_PRODUCT){
				homemade_error_msg("Not Implemented!");
		}else{
				homemade_error_msg("Please choose a correct force preparation mode!");
		}

		std::cout << "Force B preparation finish!" <<std::endl;

		m_bForceBPrepared = true;


    }
}