/*
 * generic_solver_interface.cpp
 *
 *  Created on: Nov 10, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#include "generic_solver_interface.h"

namespace carl
{

void generic_solver_interface::get_coupling_dimensions(unsigned int& M_out, unsigned int& N_out, unsigned int& M_local_out, unsigned int& N_local_out)
{
	homemade_assert_msg(m_bCouplingSet, "Must be called after setting the coupling matrix!\n");
	M_out = m_Coupling->m();
	N_out = m_Coupling->n();
	int silly_local_N;
	int silly_local_M;
	MatGetLocalSize(m_Coupling->mat(),&silly_local_M,&silly_local_N);
	M_local_out = silly_local_M;
	N_local_out = silly_local_N;
};

void generic_solver_interface::set_coupling_matrix(libMesh::PetscMatrix<libMesh::Number>& matrix)
{
	m_Coupling = &matrix;
	m_bCouplingSet = true;
};

libMesh::PetscMatrix<libMesh::Number>& generic_solver_interface::get_coupling_matrix()
{
	return *m_Coupling;
};

}