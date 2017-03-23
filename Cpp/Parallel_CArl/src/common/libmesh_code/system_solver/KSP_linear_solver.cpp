/*
 * KSP_linear_solver.cpp
 *
 *  Created on: Nov 10, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#include "KSP_linear_solver.h"
namespace carl
{

libMesh::PetscLinearSolver<libMesh::Number> * KSP_linear_solver::solver_ptr()
{
	return m_KSP_solver.get();
};

libMesh::PetscMatrix<libMesh::Number> * KSP_linear_solver::matrix_ptr()
{
	return m_Matrix;
}
void KSP_linear_solver::set_solver(libMesh::PetscMatrix<libMesh::Number>& matrix, const std::string name)
{
	m_KSP_solver->reuse_preconditioner(true);
	m_solver_name = name;
	this->set_matrix(matrix);
	m_KSP_solver->init(m_Matrix, name.c_str());
	m_perf_log_ptr = std::unique_ptr<libMesh::PerfLog>(new libMesh::PerfLog(m_solver_name));
	m_previous_time = 0;
}

void KSP_linear_solver::set_matrix(libMesh::PetscMatrix<libMesh::Number>& matrix)
{
	m_Matrix = &matrix;
	libMesh::PetscVector<libMesh::Number> * dummy_vec_1 = m_aux_vec_1.get();
	libMesh::PetscVector<libMesh::Number> * dummy_vec_2 = m_aux_vec_2.get();

	PetscInt N, local_N;
	MatGetSize(matrix.mat(),&N,NULL);
	MatGetLocalSize(matrix.mat(),&local_N,NULL);

	dummy_vec_1->init((unsigned int)N,(unsigned int)local_N);
	dummy_vec_2->init((unsigned int)N,(unsigned int)local_N);

	m_bMatrixSet = true;
};

void KSP_linear_solver::set_coupling_matrix(libMesh::PetscMatrix<libMesh::Number>& matrix)
{
	m_Coupling = &matrix;
	m_bCouplingSet = true;
};

void KSP_linear_solver::set_rhs(libMesh::PetscVector<libMesh::Number>& vector)
{
	m_rhs = &vector;
	m_bRhsSet = true;
};

libMesh::PetscMatrix<libMesh::Number>& KSP_linear_solver::get_matrix()
{
	return *m_Matrix;
}

libMesh::PetscMatrix<libMesh::Number>& KSP_linear_solver::get_coupling_matrix()
{
	return *m_Coupling;
}

libMesh::PetscVector<libMesh::Number>& KSP_linear_solver::get_rhs()
{
	return *m_rhs;
}

KSPConvergedReason KSP_linear_solver::get_converged_reason()
{
	return m_conv_reason;
}

void KSP_linear_solver::get_system_dimensions(unsigned int& M_out, unsigned int& M_local_out)
{
	homemade_assert_msg(m_bMatrixSet, "Must be called after setting the system matrix!\n");
	M_out = m_Matrix->m();
	int silly_local;
	MatGetLocalSize(m_Matrix->mat(),&silly_local,NULL);
	M_local_out = silly_local;
}

void KSP_linear_solver::get_coupling_dimensions(unsigned int& M_out, unsigned int& N_out, unsigned int& M_local_out, unsigned int& N_local_out)
{
	homemade_assert_msg(m_bMatrixSet, "Must be called after setting the coupling matrix!\n");
	M_out = m_Coupling->m();
	N_out = m_Coupling->n();
	int silly_local_N;
	int silly_local_M;
	MatGetLocalSize(m_Coupling->mat(),&silly_local_M,&silly_local_N);
	M_local_out = silly_local_M;
	N_local_out = silly_local_N;
}

void KSP_linear_solver::apply_ZMZt(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out)
{
	homemade_assert_msg(m_bMatrixSet, "Solver matrix not set yet!\n");
	homemade_assert_msg(m_bCouplingSet, "Coupling matrix not set yet!\n");

	libMesh::PetscVector<libMesh::Number> * aux_vec_rhs_ptr = m_aux_vec_1.get();
	libMesh::PetscVector<libMesh::Number> * aux_vec_sol_ptr = m_aux_vec_2.get();

	// MatMultTranspose(Mat mat,Vec x,Vec y) : y = A^t * x
	m_perf_log_ptr->push("Transpose coupling mult","Matrix operations");
	MatMultTranspose(m_Coupling->mat(),v_in.vec(),aux_vec_rhs_ptr->vec());
	m_perf_log_ptr->pop("Transpose coupling mult","Matrix operations");

	m_perf_log_ptr->push("Sys. mat mult","Matrix operations");
	m_Matrix->vector_mult(*aux_vec_sol_ptr,*aux_vec_rhs_ptr);
	m_perf_log_ptr->pop("Sys. mat mult","Matrix operations");

	// vector_mult(dest,arg)
	m_perf_log_ptr->push("Coupling mult","Matrix operations");
	m_Coupling->vector_mult(v_out,*aux_vec_sol_ptr);
	m_perf_log_ptr->pop("Coupling mult","Matrix operations");
}

void KSP_linear_solver::apply_ZMinvZt(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out)
{
	homemade_assert_msg(m_bMatrixSet, "Solver matrix not set yet!\n");
	homemade_assert_msg(m_bCouplingSet, "Coupling matrix not set yet!\n");

	libMesh::PetscVector<libMesh::Number> * aux_vec_rhs_ptr = m_aux_vec_1.get();
	libMesh::PetscVector<libMesh::Number> * aux_vec_sol_ptr = m_aux_vec_2.get();

	// MatMultTranspose(Mat mat,Vec x,Vec y) : y = A^t * x
	m_perf_log_ptr->push("Transpose coupling mult","Matrix operations");
	MatMultTranspose(m_Coupling->mat(),v_in.vec(),aux_vec_rhs_ptr->vec());
	m_perf_log_ptr->pop("Transpose coupling mult","Matrix operations");

	// solve(v_in,v_out)
	this->solve(*aux_vec_rhs_ptr,*aux_vec_sol_ptr);

	// vector_mult(dest,arg)
	m_perf_log_ptr->push("Coupling mult","Matrix operations");
	m_Coupling->vector_mult(v_out,*aux_vec_sol_ptr);
	m_perf_log_ptr->pop("Coupling mult","Matrix operations");
}

void KSP_linear_solver::apply_MinvZt(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out)
{
	homemade_assert_msg(m_bMatrixSet, "Solver matrix not set yet!\n");
	homemade_assert_msg(m_bCouplingSet, "Coupling matrix not set yet!\n");

	libMesh::PetscVector<libMesh::Number> * aux_vec_rhs_ptr = m_aux_vec_1.get();

	m_perf_log_ptr->push("Transpose coupling mult","Matrix operations");
	MatMultTranspose(m_Coupling->mat(),v_in.vec(),aux_vec_rhs_ptr->vec());
	m_perf_log_ptr->pop("Transpose coupling mult","Matrix operations");
	this->solve(*aux_vec_rhs_ptr,v_out);
}

void KSP_linear_solver::apply_ZMinv(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out)
{
	homemade_assert_msg(m_bMatrixSet, "Solver matrix not set yet!\n");
	homemade_assert_msg(m_bCouplingSet, "Coupling matrix not set yet!\n");

	libMesh::PetscVector<libMesh::Number> * aux_vec_sol_ptr = m_aux_vec_1.get();

	this->solve(v_in,*aux_vec_sol_ptr);
	m_perf_log_ptr->push("Coupling mult","Matrix operations");
	MatMult(m_Coupling->mat(),aux_vec_sol_ptr->vec(),v_out.vec());
	m_perf_log_ptr->pop("Coupling mult","Matrix operations");
}

// Implementations of virtual methods
void KSP_linear_solver::solve(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out)
{
	homemade_assert_msg(m_bMatrixSet, "Solver matrix not set yet!\n");
	m_perf_log_ptr->push("KSP solver");
	m_KSP_solver->solve(*m_Matrix,v_out,v_in,m_KSP_eps,m_KSP_iter_max);
	m_perf_log_ptr->pop("KSP solver");

	int dummy_int = 0;
	KSPGetIterationNumber(m_KSP_solver->ksp(),&dummy_int);
	m_KSP_tot_iter += dummy_int;

	if(bLogTiming)
	{
		m_solve_data_time = m_perf_log_ptr->get_perf_data("KSP solver");

		if(m_comm->rank() == 0)
		{
			std::ofstream timing_data(m_log_filename, std::ios::app);

			timing_data << m_solve_data_time.count << " " 
		        	    << m_solve_data_time.tot_time - m_previous_time << " " 
				    << dummy_int << " ";
			if(dummy_int != 0)
			{
				timing_data << (m_solve_data_time.tot_time - m_previous_time)/dummy_int << std::endl;
			}
			else
			{
				timing_data << "INF" << std::endl;
			}
			timing_data.close();
			m_previous_time = m_solve_data_time.tot_time;
		}
	}

	KSPGetConvergedReason(m_KSP_solver->ksp(),&m_conv_reason);
}

void KSP_linear_solver::solve(libMesh::PetscVector<libMesh::Number>& v_out)
{
	homemade_assert_msg(m_bRhsSet, "Solver RHS not set yet!\n");
	homemade_assert_msg(m_bMatrixSet, "Solver matrix not set yet!\n");
	m_perf_log_ptr->push("KSP solver");
	m_KSP_solver->solve(*m_Matrix,v_out,*m_rhs,m_KSP_eps,m_KSP_iter_max);
	int dummy_int = 0;
        KSPGetIterationNumber(m_KSP_solver->ksp(),&dummy_int);
        m_KSP_tot_iter += dummy_int;
	m_perf_log_ptr->pop("KSP solver");
	KSPGetConvergedReason(m_KSP_solver->ksp(),&m_conv_reason);
}

void KSP_linear_solver::print_type()
{
	KSPType solver_type_string;
	KSPGetType(m_KSP_solver->ksp(),&solver_type_string);
	std::cout 	<< "|        Solver type : " << solver_type_string << std::endl;
	std::cout 	<< "|        PC     type : " << m_KSP_solver->preconditioner_type() << std::endl;
	std::cout   << "|" << std::endl;
}

void KSP_linear_solver::calculate_pseudo_inverse(const std::string& filename)
{
	Vec dummy_vec;
	Mat matrix_out;

	PetscInt M, N;

	MatCreateVecs(m_Matrix->mat(),&dummy_vec,PETSC_NULL);
	MatGetSize(m_Matrix->mat(),&M,&N);

	Vec *dummy_Id;
	Vec *dummy_sol;

	VecDuplicateVecs(dummy_vec,N,&dummy_Id);
	VecDuplicateVecs(dummy_vec,N,&dummy_sol);

	for(int iii = 0; iii < N; ++iii)
	{
		VecSet(dummy_Id[iii],0);
		VecSetValue(dummy_Id[iii],iii,1,INSERT_VALUES);
		VecSet(dummy_sol[iii],0);

		libMesh::PetscVector<libMesh::Number> dummy_Id_libmesh(dummy_Id[iii],*m_comm);
		libMesh::PetscVector<libMesh::Number> dummy_sol_libmesh(dummy_sol[iii],*m_comm);

		this->solve(dummy_Id_libmesh,dummy_sol_libmesh);
	}

	create_PETSC_dense_matrix_from_vectors(dummy_sol,N,matrix_out);

	PetscViewer    viewer;
	PetscViewerASCIIOpen(m_comm->get(),filename.c_str(),&viewer);
	PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);

	MatView(matrix_out,viewer);

	PetscViewerDestroy(&viewer);

	MatDestroy(&matrix_out);
	VecDestroyVecs(N,&dummy_Id);
	VecDestroyVecs(N,&dummy_sol);
	VecDestroy(&dummy_vec);
}

void KSP_linear_solver::get_perf_log_timing(double& solve_time, int& solve_calls)
{
	libMesh::PerfData performance_data = m_perf_log_ptr->get_perf_data("KSP solver");
	solve_time = performance_data.tot_time_incl_sub;
	solve_calls = performance_data.count;
};

int KSP_linear_solver::get_iter_total()
{
	return m_KSP_tot_iter;
};

};
