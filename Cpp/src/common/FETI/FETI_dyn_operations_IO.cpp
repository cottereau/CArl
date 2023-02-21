/*
 * \file FETI_dyn_operations_IO.cpp
 *
 *  Created on: Nov 23,2021
 *      Author: Chensheng Luo
 * 
 * \brief **DYN**   functions responsible for result output, copy/paste and others in Solving steps
 */

#include "FETI_dyn_operations.h"

namespace carl
{

	void FETI_Dyn_Operations::scale_copy_PETSC_matrix(std::string copy_path, std::string paste_path,PetscScalar scale){
		Mat matrix;
		MatCreate(m_comm.get(),&matrix);
		carl::read_PETSC_matrix(matrix,copy_path,m_comm.get());
		MatScale(matrix,scale);
		carl::write_PETSC_matrix(matrix,paste_path,0,m_comm.get(),1);
		
		#ifdef PRINT_MATLAB_DEBUG
			libMesh::PetscMatrix<libMesh::Number> matrix_lib(matrix,m_comm);
   			matrix_lib.print_matlab(paste_path+".m");
		#endif
	}

	void FETI_Dyn_Operations::scale_copy_PETSC_vector(std::string copy_path, std::string paste_path,PetscScalar scale){
		Vec vector;
		VecCreate(m_comm.get(),&vector);
		carl::read_PETSC_vector(vector,copy_path,m_comm.get());
		VecScale(vector,scale);
		carl::write_PETSC_vector(vector,paste_path,0,m_comm.get(),1);
		
		#ifdef PRINT_MATLAB_DEBUG
			libMesh::PetscVector<libMesh::Number> vector_lib(vector,m_comm);
   			vector_lib.print_matlab(paste_path+".m");
		#endif
	}

	void FETI_Dyn_Operations::copy_PETSC_vector(std::string copy_path, std::string paste_path){
		this->scale_copy_PETSC_vector(copy_path,paste_path,1);
	}

	void FETI_Dyn_Operations::set_NAN_vector(std::string& vector_path){
		Vec vector;
		VecCreate(m_comm.get(),&vector);
		carl::read_PETSC_vector(vector,vector_path,m_comm.get());
		VecSet(vector,NAN);
		carl::write_PETSC_vector(vector,vector_path,0,m_comm.get(),1);
	}

	


	void FETI_Dyn_Operations::output_A_result(int index){
		#ifdef OUTPUT_ACC_DEBUG
      		this->copy_PETSC_vector(vector_A.this_acc,m_result_folder_path+"/acc_A_"+std::to_string(index)+".petscvec");
		#endif
      	#ifdef OUTPUT_SPEED_DEBUG
			this->copy_PETSC_vector(vector_A.this_speed,m_result_folder_path+"/speed_A_"+std::to_string(index)+".petscvec");
		#endif
		this->copy_PETSC_vector(vector_A.this_disp,m_result_folder_path+"/disp_A_"+std::to_string(index)+".petscvec");

		#ifdef OUTPUT_FREE_DEBUG
			#ifdef OUTPUT_ACC_DEBUG
	      		this->copy_PETSC_vector(vector_A.this_acc_free,m_result_folder_path+"/acc_free_A_"+std::to_string(index)+".petscvec");
			#endif
      		#ifdef OUTPUT_SPEED_DEBUG
				this->copy_PETSC_vector(vector_A.this_speed_free,m_result_folder_path+"/speed_free_A_"+std::to_string(index)+".petscvec");
			#endif
			this->copy_PETSC_vector(vector_A.this_disp_free,m_result_folder_path+"/disp_free_A_"+std::to_string(index)+".petscvec");
		#endif

		#ifdef OUTPUT_LINK_DEBUG
			#ifdef OUTPUT_ACC_DEBUG
				this->copy_PETSC_vector(vector_A.this_acc_link,m_result_folder_path+"/acc_link_A_"+std::to_string(index)+".petscvec");
			#endif
      		#ifdef OUTPUT_SPEED_DEBUG
				this->copy_PETSC_vector(vector_A.this_speed_link,m_result_folder_path+"/speed_link_A_"+std::to_string(index)+".petscvec");
			#endif
			this->copy_PETSC_vector(vector_A.this_disp_link,m_result_folder_path+"/disp_link_A_"+std::to_string(index)+".petscvec");
		#endif
	}

	void FETI_Dyn_Operations::move_to_prev_A(){
      	this->copy_PETSC_vector(vector_A.this_acc,vector_A.prev_acc);
      	this->copy_PETSC_vector(vector_A.this_disp,vector_A.prev_disp);
      	this->copy_PETSC_vector(vector_A.this_disp_free,vector_A.prev_disp_free);
      	this->copy_PETSC_vector(vector_A.this_speed,vector_A.prev_speed);
      	this->copy_PETSC_vector(vector_A.this_speed_free,vector_A.prev_speed_free);
      	m_bAMovedToPrev = true;
	}

	void FETI_Dyn_Operations::delete_A_this_vector(){
		this->set_NAN_vector(vector_A.rhs_free);
		this->set_NAN_vector(vector_A.rhs_link);

		this->set_NAN_vector(vector_A.this_acc);
		this->set_NAN_vector(vector_A.this_acc_free);
		this->set_NAN_vector(vector_A.this_acc_link);

		this->set_NAN_vector(vector_A.this_disp);
		this->set_NAN_vector(vector_A.this_disp_free);
		this->set_NAN_vector(vector_A.this_disp_link);

		this->set_NAN_vector(vector_A.this_speed);
		this->set_NAN_vector(vector_A.this_speed_free);
		this->set_NAN_vector(vector_A.this_speed_link);
	}

	void FETI_Dyn_Operations::output_B_result(int index){
      	#ifdef OUTPUT_ACC_DEBUG
      		this->copy_PETSC_vector(vector_B.this_acc,m_result_folder_path+"/acc_B_"+std::to_string(index)+".petscvec");
      	#endif
      	#ifdef OUTPUT_SPEED_DEBUG
			this->copy_PETSC_vector(vector_B.this_speed,m_result_folder_path+"/speed_B_"+std::to_string(index)+".petscvec");
		#endif
		this->copy_PETSC_vector(vector_B.this_disp,m_result_folder_path+"/disp_B_"+std::to_string(index)+".petscvec");


		#ifdef OUTPUT_FREE_DEBUG
			#ifdef OUTPUT_ACC_DEBUG
	      		this->copy_PETSC_vector(vector_B.this_acc_free,m_result_folder_path+"/acc_free_B_"+std::to_string(index)+".petscvec");
			#endif
      		#ifdef OUTPUT_SPEED_DEBUG
				this->copy_PETSC_vector(vector_B.this_speed_free,m_result_folder_path+"/speed_free_B_"+std::to_string(index)+".petscvec");
			#endif
			this->copy_PETSC_vector(vector_B.this_disp_free,m_result_folder_path+"/disp_free_B_"+std::to_string(index)+".petscvec");
		#endif

		#ifdef OUTPUT_LINK_DEBUG
			#ifdef OUTPUT_ACC_DEBUG
				this->copy_PETSC_vector(vector_B.this_acc_link,m_result_folder_path+"/acc_link_B_"+std::to_string(index)+".petscvec");
			#endif
      		#ifdef OUTPUT_SPEED_DEBUG
				this->copy_PETSC_vector(vector_B.this_speed_link,m_result_folder_path+"/speed_link_B_"+std::to_string(index)+".petscvec");
			#endif
			this->copy_PETSC_vector(vector_B.this_disp_link,m_result_folder_path+"/disp_link_B_"+std::to_string(index)+".petscvec");
		#endif
	}

	void FETI_Dyn_Operations::move_to_prev_B(){
      	this->copy_PETSC_vector(vector_B.this_acc,vector_B.prev_acc);
      	this->copy_PETSC_vector(vector_B.this_disp,vector_B.prev_disp);
      	this->copy_PETSC_vector(vector_B.this_disp_free,vector_B.prev_disp_free);
      	this->copy_PETSC_vector(vector_B.this_speed,vector_B.prev_speed);
      	this->copy_PETSC_vector(vector_B.this_speed_free,vector_B.prev_speed_free);
      	m_bBMovedToPrev = true;
	}

	void FETI_Dyn_Operations::delete_B_this_vector(){
		this->set_NAN_vector(vector_B.rhs_free);
		this->set_NAN_vector(vector_B.rhs_link);

		this->set_NAN_vector(vector_B.this_acc);
		this->set_NAN_vector(vector_B.this_acc_free);
		this->set_NAN_vector(vector_B.this_acc_link);

		this->set_NAN_vector(vector_B.this_disp);
		this->set_NAN_vector(vector_B.this_disp_free);
		this->set_NAN_vector(vector_B.this_disp_link);

		this->set_NAN_vector(vector_B.this_speed);
		this->set_NAN_vector(vector_B.this_speed_free);
		this->set_NAN_vector(vector_B.this_speed_link);
	}

	void FETI_Dyn_Operations::delete_coupling_vector(){
		this->set_NAN_vector(m_coupling_vector);
		this->set_NAN_vector(m_rhs_interpolation);
	}

	// // [DYN-CG]
	// void FETI_Dyn_Operations::prepare_CG_free_result(std::string& depature_path,
	// 	std::string destination_path){
	// 	this->copy_PETSC_vector(vector_A.inter_disp_free,destination_path+"/ext_solver_u0_A_sys_sol_vec.petscvec");
	// 	this->copy_PETSC_vector(vector_B.this_disp_free,destination_path+"/ext_solver_u0_B_sys_sol_vec.petscvec");

	// }

	// // [DYN-CG]
	// void FETI_Dyn_Operations::prepare_CG_scaled_matrix(std::string& M_path_A,std::string& M_path_B,
    // 	std::string destination_path,carl::NewmarkParams* newmark_A,carl::NewmarkParams* newmark_B){
	// 	this->scale_copy_PETSC_matrix(M_path_A,destination_path+"/Mtilde_A.petscmat", 1.0/(newmark_A->beta*newmark_A->deltat*newmark_A->deltat));
	// 	this->scale_copy_PETSC_matrix(M_path_B,destination_path+"/Mtilde_B.petscmat", 1.0/(newmark_B->beta*newmark_B->deltat*newmark_B->deltat));
	// }

	// // [DYN-CG]
	// void FETI_Dyn_Operations::prepare_CG_scaled_vector(std::string& rb_path,int nb_rb,std::string force_path,
    // 	std::string destination_path,carl::NewmarkParams* newmark){
	// 	this->scale_copy_PETSC_vector(force_path,destination_path+"/current_force.petscvec",-1);

	// 	for(int iii = 0; iii < nb_rb; ++iii){
	// 		std::string input_filename = rb_path + "_" + std::to_string(iii) + "_n_" + std::to_string(nb_rb) + ".petscvec";
	// 		std::string output_filename = destination_path + "/prepared_rb_vector_" + std::to_string(iii) + "_n_" + std::to_string(nb_rb) + ".petscvec";
	// 		this->scale_copy_PETSC_vector(input_filename,output_filename,-(newmark->beta*newmark->deltat*newmark->deltat));
	// 	}
	// }

	void FETI_Dyn_Operations::export_calculation_time(int index,std::string stage)
	{

  		// ONLY write in proc 0!
  		std::ofstream time_data;

  		if(m_comm.rank() == 0){
    		std::time_t now=time(0);
    		char* dt=ctime(&now);

    		// Export the scalar data
    		time_data.open(m_scratch_folder_path + "/Time_data.dat",std::ofstream::app);
    		time_data.precision(5);
    		time_data << std::left << std::setw(7) << index << std::setw(20) << stage << std::setw(25) << dt << std::endl;
    		
    		// std::FILE *fh ;
    		// std::string filename=m_scratch_folder_path+"/CG_solver/FETI_convergence.dat";
			// fh = fopen(filename.c_str (),"r");// To examine if there is a cv file

			// if(fh!=NULL){
			// 	int count;
			// 	char line[256]={0};
			// 	while(std::fgets(line,255,fh)){
			// 		count++;
			// 	}
			// 	//There is a cv file
			// 	time_data << " "<< count;
			// }
    		// time_data << std::endl;
    		// time_data.close();

    	}
    	std::cout << index << " "<< stage << " Time stored!"<< std::endl;
	}
}