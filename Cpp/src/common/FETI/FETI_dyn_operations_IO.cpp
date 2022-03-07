#include "FETI_dyn_operations.h"

namespace carl
{

	void FETI_Dyn_Operations::copy_PETSC_vector(std::string copy_path, std::string paste_path){
		Vec vector;
		VecCreate(m_comm.get(),&vector);
		carl::read_PETSC_vector(vector,copy_path,m_comm.get());
		carl::write_PETSC_vector(vector,paste_path+".petscvec",0,m_comm.get(),1);
		
		#ifdef PRINT_MATLAB_DEBUG
			libMesh::PetscVector<libMesh::Number> vector_lib(vector,m_comm);
   			vector_lib.print_matlab(paste_path+".m");
		#endif
	}

	void FETI_Dyn_Operations::output_B_result(std::string result_folder_path,
      int index){
      	this->copy_PETSC_vector(m_scratch_folder_path+"/this_acc_B.petscvec",result_folder_path+"/acc_B_"+std::to_string(index));
		this->copy_PETSC_vector(m_scratch_folder_path+"/this_speed_B.petscvec",result_folder_path+"/speed_B_"+std::to_string(index));
		this->copy_PETSC_vector(m_scratch_folder_path+"/this_disp_B.petscvec",result_folder_path+"/disp_B_"+std::to_string(index));
	}

	void FETI_Dyn_Operations::move_to_prev_B(){
      	this->copy_PETSC_vector(m_scratch_folder_path+"/this_acc_B.petscvec",m_scratch_folder_path+"/prev_acc_B");
      	this->copy_PETSC_vector(m_scratch_folder_path+"/this_speed_B.petscvec",m_scratch_folder_path+"/prev_speed_B");
      	this->copy_PETSC_vector(m_scratch_folder_path+"/this_disp_B.petscvec",m_scratch_folder_path+"/prev_disp_B");

	}

	void FETI_Dyn_Operations::output_A_result(std::string result_folder_path,
      int index){
      	this->copy_PETSC_vector(m_scratch_folder_path+"/this_acc_A.petscvec",result_folder_path+"/acc_A_"+std::to_string(index));
		this->copy_PETSC_vector(m_scratch_folder_path+"/this_speed_A.petscvec",result_folder_path+"/speed_A_"+std::to_string(index));
		this->copy_PETSC_vector(m_scratch_folder_path+"/this_disp_A.petscvec",result_folder_path+"/disp_A_"+std::to_string(index));
	}

	void FETI_Dyn_Operations::move_to_prev_A(){
      	this->copy_PETSC_vector(m_scratch_folder_path+"/this_acc_A.petscvec",m_scratch_folder_path+"/prev_acc_A");
      	this->copy_PETSC_vector(m_scratch_folder_path+"/this_speed_A.petscvec",m_scratch_folder_path+"/prev_speed_A");
      	this->copy_PETSC_vector(m_scratch_folder_path+"/this_disp_A.petscvec",m_scratch_folder_path+"/prev_disp_A");

	}
}