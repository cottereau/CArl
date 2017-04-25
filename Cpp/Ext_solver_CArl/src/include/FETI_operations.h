/*
 * FETI_operations.h
 *
 *  Created on: Apr 23, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef FETI_OPERATIONS_H_
#define FETI_OPERATIONS_H_

#include "carl_headers.h"
#include "PETSC_matrix_operations.h"

namespace carl
{
class	FETI_Operations
{
protected:

	enum { maxVecLimit = 6 };
	libMesh::Parallel::Communicator& m_comm;

	PetscInt	m_C_R_micro_M, m_C_R_micro_N, m_C_R_micro_M_local, m_C_R_micro_N_local;
	PetscInt	m_C_R_BIG_M, m_C_R_BIG_N, m_C_R_BIG_M_local, m_C_R_BIG_N_local;
	PetscInt	m_C_RR_M, m_C_RR_M_local;

	bool 	    m_bNullVecsSet;
	PetscInt    m_null_nb_vecs;
	Vec        	m_null_vecs[maxVecLimit];
	Vec		    m_null_coupled_vecs[maxVecLimit];

	bool 	    m_bC_R_BIG_MatrixSet;
	bool 	    m_bC_R_micro_MatrixSet;
	bool 	    m_bC_RR_MatrixSet;
	bool 	    m_bCouplingMatricesSet;
	
	Mat m_C_R_BIG;
	Mat m_C_R_micro;
	Mat m_C_RR;

	Mat m_RITRI_mat;
	Mat m_inv_RITRI_mat;

	RBModesSystem m_RB_modes_system;

	FETI_Operations();


public:
	FETI_Operations(libMesh::Parallel::Communicator& comm) :
		m_comm { comm },
		m_C_R_micro_M { -1},
		m_C_R_micro_N { -1},
		m_C_R_micro_M_local { -1},
		m_C_R_micro_N_local { -1},
		m_C_R_BIG_M { -1},
		m_C_R_BIG_N { -1},
		m_C_R_BIG_M_local { -1},
		m_C_R_BIG_N_local { -1},
		m_C_RR_M { -1},
		m_C_RR_M_local { -1},
		m_bNullVecsSet { false },
		m_null_nb_vecs { -1 },
		m_bC_R_BIG_MatrixSet  { false },
		m_bC_R_micro_MatrixSet  { false },
		m_bC_RR_MatrixSet  { false },
		m_bCouplingMatricesSet  { false },
		m_RB_modes_system { RBModesSystem::MICRO }
	{
	};

	~FETI_Operations()
	{
		if(m_bC_R_BIG_MatrixSet)
			MatDestroy(&m_C_R_BIG);
		if(m_bC_R_micro_MatrixSet)
			MatDestroy(&m_C_R_micro);
		if(m_bC_RR_MatrixSet)
			MatDestroy(&m_C_RR);
		if(m_bNullVecsSet)
		{
			for(int iii = 0; iii < m_null_nb_vecs; ++iii)
				VecDestroy(&m_null_vecs[iii]);
		}
	};
	
	void set_null_space_vecs_micro(const std::string& input_filename_base, const std::string& output_filename_base, int nb_of_vecs);

	void set_coupling_matrix_R_micro(const std::string& filename);

	void set_coupling_matrix_R_BIG(const std::string& filename);

	void set_coupling_matrix_RR(const std::string& filename);

	void set_coupling_matrices(const std::string& filename_base);

	void calculate_phi_0(const std::string& force_micro_path, const std::string& scratch_folder_path);
};
}

#endif /* FETI_OPERATIONS_H_ */