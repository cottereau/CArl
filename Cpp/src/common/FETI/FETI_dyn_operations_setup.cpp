#include "FETI_dyn_operations.h"

namespace carl
{

	void FETI_Dyn_Operations::init_prepare_force_vector_by_modal_and_sinus(std::string& force_A_path,
		std::string& force_B_path,
		double amplitude,
		double frequency,
		double initialPhase,
		double deltat,
		int innerTimes,
		int outerTimes)
	{
		Vec modalA,modalB;
		Vec forceA,forceB;

		VecCreate(m_comm.get(),&modalA);
		VecCreate(m_comm.get(),&modalB);
		carl::read_PETSC_vector(modalA,force_A_path,m_comm.get());
		carl::read_PETSC_vector(modalB,force_B_path,m_comm.get());
		VecDuplicate(modalA,&forceA);
		VecDuplicate(modalB,&forceB);

		for(int i = 1; i <= innerTimes*outerTimes; i++){
			if(i%innerTimes==0){
				VecCopy(modalA,forceA);
				VecScale(forceA,amplitude*sin(2*M_PI*frequency*deltat*i+initialPhase));
				carl::write_PETSC_vector(forceA,m_scratch_folder_path+"/force_A/force_"+std::to_string(i)+".petscvec",0,m_comm.get(),1);
			}
			VecCopy(modalB,forceB);
			VecScale(forceB,amplitude*sin(2*M_PI*frequency*deltat*i+initialPhase));
			carl::write_PETSC_vector(forceB,m_scratch_folder_path+"/force_B/force_"+std::to_string(i)+".petscvec",0,m_comm.get(),1);
		}

		VecCopy(modalA,forceA);
		VecScale(forceA,amplitude*sin(2*M_PI*frequency*deltat+initialPhase));
		carl::write_PETSC_vector(forceA,m_scratch_folder_path+"/rhs_vec_A_free.petscvec",0,m_comm.get(),1);
		

		VecCopy(modalB,forceB);
		VecScale(forceB,amplitude*sin(2*M_PI*frequency*deltat+initialPhase));
		carl::write_PETSC_vector(forceB,m_scratch_folder_path+"/rhs_vec_B_free.petscvec",0,m_comm.get(),1);

		VecDestroy(&modalA);
		VecDestroy(&modalB);
		VecDestroy(&forceA);
		VecDestroy(&forceB);
	}

}