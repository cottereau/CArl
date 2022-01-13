#include "FETI_dyn_operations.h"

namespace carl
{
	void FETI_Dyn_Operations::init_prepare_rhs_vector(std::string& rhs_vector_A_path, 
      std::string& rhs_vector_B_path)
	{

		this->copy_PETSC_vector(rhs_vector_A_path,m_scratch_folder_path+"/rhs_vec_A_free");
		this->copy_PETSC_vector(rhs_vector_B_path,m_scratch_folder_path+"/rhs_vec_B_free");
	}


}