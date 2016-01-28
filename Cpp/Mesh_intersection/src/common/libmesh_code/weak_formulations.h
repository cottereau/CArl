/*
 * weak_formulations.h
 *
 *  Created on: Jan 7, 2016
 *      Author: Thiago Milanetto Schlittler
 *
 *		This file contains the headers of the functions used to add terms to the
 *		FEM matrices, including the coupling matrices used for the Arlequin
 *		method
 */

#ifndef COMMON_LIBMESH_CODE_WEAK_FORMULATIONS_H_
#define COMMON_LIBMESH_CODE_WEAK_FORMULATIONS_H_

#include "common_header.h"
#include "common_header_libmesh.h"

// --- Matrices
void Mass(	libMesh::DenseMatrix<libMesh::Number>& Mass,
			unsigned int qp,
			const std::vector<std::vector<libMesh::Real> >& phi,
			const unsigned int n_dofs,
			const std::vector<libMesh::Real>& JxW,
			const libMesh::Number cte = 1
			);

// --- Vectors

// --- Arlequin coupling matrices
void L2_Coupling(	libMesh::DenseMatrix<libMesh::Number>& Coupl,
					unsigned int qp,
					const std::vector<std::vector<libMesh::Real> >& phi_sysA,
					const std::vector<std::vector<libMesh::Real> >& phi_sysB,
					const unsigned int n_dofs_sysA,
					const unsigned int n_dofs_sysB,
					const std::vector<libMesh::Real>& JxW,
					const libMesh::Number cte
					);

void L2_Coupling(	libMesh::DenseSubMatrix<libMesh::Number>& Coupl,
					unsigned int qp,
					const std::vector<std::vector<libMesh::Real> >& phi_sysA,
					const std::vector<std::vector<libMesh::Real> >& phi_sysB,
					const unsigned int n_dofs_sysA,
					const unsigned int n_dofs_sysB,
					const std::vector<libMesh::Real>& JxW,
					const libMesh::Number cte
					);

#endif /* COMMON_LIBMESH_CODE_WEAK_FORMULATIONS_H_ */
