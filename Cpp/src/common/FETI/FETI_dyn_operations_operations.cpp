#include "FETI_dyn_operations.h"

namespace carl
{

	void FETI_Dyn_Operations::prepare_rhs_vector(double beta,
      double deltat,
      std::string force_vector_folder, 
      int index,
      std::string this_acc_path,
      std::string this_speed_path,
      std::string this_displace_path,
      std::string stiffness_matrix_path,
      std::string output_rhs_path){

			

			Vec force;
			Vec displacement;
      Vec this_speed;
			Vec this_acceleration;
			Mat stiffness;
			Vec rhs;
			VecCreate(m_comm.get(),&force);
			VecCreate(m_comm.get(),&displacement);
			VecCreate(m_comm.get(),&this_speed);
			VecCreate(m_comm.get(),&this_acceleration);
			MatCreate(m_comm.get(),&stiffness);
			VecCreate(m_comm.get(),&rhs);

			carl::read_PETSC_vector(displacement,this_displace_path,m_comm.get());
			carl::read_PETSC_vector(this_speed,this_speed_path,m_comm.get());
			carl::read_PETSC_vector(this_acceleration,this_acc_path,m_comm.get());
			carl::read_PETSC_matrix(stiffness,stiffness_matrix_path,m_comm.get());

			//predictor=displacement+deltat*this_speed+(1/2-beta)*deltat^2*this_acceleration
			VecAXPY(displacement,deltat,this_speed);
			VecAXPY(displacement,(1/2-beta)*deltat*deltat,this_acceleration);

			MatMult(stiffness,displacement,rhs);

			std::string force_vector_path = force_vector_folder+std::to_string(index);

		std::FILE *fh ;
		fh = fopen(force_vector_path.c_str (),"r");// To examine if this is a force input file

		if(fh!=NULL){
			//There is a force input file
			
			carl::read_PETSC_vector(force,force_vector_path,m_comm.get());
			VecAXPY(rhs,1,force);
			}

			carl::write_PETSC_vector(rhs,output_rhs_path,0,m_comm.get(),1);

			VecDestroy(&force);
			VecDestroy(&displacement);
			VecDestroy(&this_speed);
			VecDestroy(&this_acceleration);
			MatDestroy(&stiffness);
			VecDestroy(&rhs);
			
	}
	void FETI_Dyn_Operations::Newmark_speed_free(double gamma,
      double deltat,
      std::string prev_acc_path,
      std::string this_acc_path,
      std::string input_speed_path,
      std::string output_speed_path){

      	// Get acceleration
			Vec speed;
			Vec prev_acceleration;
			Vec this_acceleration;
			VecCreate(m_comm.get(),&speed);
			VecCreate(m_comm.get(),&prev_acceleration);
			VecCreate(m_comm.get(),&this_acceleration);
			carl::read_PETSC_vector(this_acceleration,this_acc_path,m_comm.get());


		std::FILE *fh ;
		fh = fopen(input_speed_path.c_str (),"r");// To examine if this is a "previous" file
		
		if(fh!=NULL){
			//Not the first moment
			carl::read_PETSC_vector(prev_acceleration,prev_acc_path,m_comm.get());
			carl::read_PETSC_vector(speed,input_speed_path,m_comm.get());
			
			//predictor=speed+(1-gamma)*deltat*prev_acceleration
			VecAXPY(speed,(1-gamma)*deltat,prev_acceleration);

			//new_speed=predictor+gamma*deltat*this_acceleration
			VecAXPY(speed,gamma*deltat,this_acceleration);
		}else{
			//No prev_speed vector exists! We need create one.
			VecCopy(this_acceleration,speed);
			VecScale(speed,gamma*deltat);
		}

		//Calculate 
		carl::write_PETSC_vector(speed,output_speed_path,0,m_comm.get(),1);

		VecDestroy(&speed);
		VecDestroy(&prev_acceleration);
		VecDestroy(&this_acceleration);
	}


  void FETI_Dyn_Operations::Newmark_displacement_free(double beta,
      double deltat,
      std::string prev_acc_path,
      std::string this_acc_path,
      std::string prev_speed_path,
      std::string input_displace_path,
      std::string output_displace_path){

      Vec displacement;
      Vec prev_speed;
      Vec prev_acceleration;
			Vec this_acceleration;
			VecCreate(m_comm.get(),&displacement);
			VecCreate(m_comm.get(),&prev_speed);
			VecCreate(m_comm.get(),&prev_acceleration);
			VecCreate(m_comm.get(),&this_acceleration);
			carl::read_PETSC_vector(this_acceleration,this_acc_path,m_comm.get());

		std::FILE *fh;
		fh = fopen(input_displace_path.c_str (),"r");// To examine if this is a "previous" file

		if(fh!=NULL){
			//Not the first moment
			carl::read_PETSC_vector(displacement,input_displace_path,m_comm.get());
			carl::read_PETSC_vector(prev_speed,prev_speed_path,m_comm.get());
			carl::read_PETSC_vector(prev_acceleration,prev_acc_path,m_comm.get());
			
			//predictor=displacement+deltat*prev_speed+(1/2-beta)*deltat^2*prev_acceleration
			VecAXPY(displacement,deltat,prev_speed);
			VecAXPY(displacement,(1/2-beta)*deltat*deltat,prev_acceleration);

			//new_displacement=predictor+beta*deltat*deltat*this_acceleration
			VecAXPY(displacement,beta*deltat*deltat,this_acceleration);
		}else{
			//No prev_displacement vector exists! We need create one.
			VecCopy(this_acceleration,displacement);
			VecScale(displacement,beta*deltat*deltat);
		}

		//Calculate 
		carl::write_PETSC_vector(displacement,output_displace_path,0,m_comm.get(),1);

		VecDestroy(&displacement);
		VecDestroy(&prev_speed);
		VecDestroy(&prev_acceleration);
		VecDestroy(&this_acceleration);
    }

  void FETI_Dyn_Operations::interpolate_A_acceleration(int jjj,int m){
  	Vec prev_acceleration;
  	Vec this_acceleration;
  	std::string prev_acc_path=m_scratch_folder_path+"/prev_acc_A.petscvec";
  	std::string this_acc_path=m_scratch_folder_path+"/this_acc_A_free.petscvec";

  	VecCreate(m_comm.get(),&prev_acceleration);
		VecCreate(m_comm.get(),&this_acceleration);
		carl::read_PETSC_vector(this_acceleration,this_acc_path,m_comm.get());

		std::FILE *fh;
		fh = fopen(prev_acc_path.c_str (),"r");// To examine if this is a "previous" file

		if(fh!=NULL){
			//Not the first moment
			carl::read_PETSC_vector(prev_acceleration,prev_acc_path,m_comm.get());
			
			//U_free_j=(1-j/m)*prev_acc+j/m*this_acc
			VecAXPBY(this_acceleration,(1-jjj/m),jjj/m,prev_acceleration);
		}else{
			VecScale(this_acceleration,jjj/m);
		}

		carl::write_PETSC_vector(this_acceleration,m_scratch_folder_path+"/inter_acc_A_free",0,m_comm.get(),1);

		VecDestroy(&prev_acceleration);
		VecDestroy(&this_acceleration);
  }

  void FETI_Dyn_Operations::rhs_interpolation(std::string coupling_matrix_A_path,std::string coupling_matrix_B_path){
  		Mat coupling_matrix_A;
  		Mat coupling_matrix_B;
  		Vec acceleration_A;
  		Vec acceleration_B;
  		Vec result1;
  		Vec result2;

  		MatCreate(m_comm.get(),&coupling_matrix_A);
  		MatCreate(m_comm.get(),&coupling_matrix_B);
  		VecCreate(m_comm.get(),&acceleration_A);
  		VecCreate(m_comm.get(),&acceleration_B);
  		VecCreate(m_comm.get(),&result1);
  		VecCreate(m_comm.get(),&result2);

  		carl::read_PETSC_vector(acceleration_A,m_scratch_folder_path+"/inter_acc_A_free.petscvec");
  		carl::read_PETSC_vector(acceleration_B,m_scratch_folder_path+"/this_acc_A_free.petscvec");
  		carl::read_PETSC_matrix(coupling_matrix_A,coupling_matrix_A_path);
  		carl::read_PETSC_matrix(coupling_matrix_B,coupling_matrix_B_path);

  		MatMult(coupling_matrix_A,acceleration_A,result1);
  		MatMultAdd(coupling_matrix_B,acceleration_B,result1,result2);
  		VecScale(result2, -1);


  		carl::write_PETSC_vector(result2,m_scratch_folder_path+"/rhs_interpolation_vec",0,m_comm.get(),1);
    

			MatDestroy(&coupling_matrix_A);
  		MatDestroy(&coupling_matrix_B);
  		VecDestroy(&acceleration_A);
  		VecDestroy(&acceleration_B);
  		VecDestroy(&result1);
  		VecDestroy(&result2);
    }
  void FETI_Dyn_Operations::rhs_link(std::string coupling_matrix_path,std::string interpolation_vector_path,std::string output_path){
  		Mat coupling_matrix;
  		Vec interpolation_vector;
  		Vec result;

  		MatCreate(m_comm.get(),&coupling_matrix);
  		VecCreate(m_comm.get(),&interpolation_vector);
  		VecCreate(m_comm.get(),&result);

  		carl::read_PETSC_matrix(coupling_matrix,coupling_matrix_path);
  		carl::read_PETSC_vector(interpolation_vector,interpolation_vector_path);

  		MatMultTranspose(coupling_matrix,interpolation_vector,result);

  		carl::write_PETSC_vector(result,output_path,0,m_comm.get(),1);

  		MatDestroy(&coupling_matrix);
  		VecDestroy(&interpolation_vector);
  		VecDestroy(&result);
  }

  void FETI_Dyn_Operations::Newmark_speed_link(double gamma,
      double deltat,
      std::string this_acc_path,
      std::string output_speed_path){
  		Vec this_acceleration;

  		VecCreate(m_comm.get(),&this_acceleration);

  		carl::read_PETSC_vector(this_acceleration,this_acc_path);

  		VecScale(this_acceleration,gamma*deltat);

  		carl::write_PETSC_vector(this_acceleration,output_speed_path,0,m_comm.get(),1);

  		VecDestroy(&this_acceleration);

  }
    void FETI_Dyn_Operations::Newmark_displacement_link(double beta,
      double deltat,
      std::string this_acc_path,
      std::string output_displace_path){
    	Vec this_acceleration;

  		VecCreate(m_comm.get(),&this_acceleration);

  		carl::read_PETSC_vector(this_acceleration,this_acc_path);

  		VecScale(this_acceleration,beta*deltat*deltat);

  		carl::write_PETSC_vector(this_acceleration,output_displace_path,0,m_comm.get(),1);

  		VecDestroy(&this_acceleration);

    }

    void FETI_Dyn_Operations::add_free_link(std::string free_term_path,
      std::string link_term_path,
      std::string output_path){
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
}