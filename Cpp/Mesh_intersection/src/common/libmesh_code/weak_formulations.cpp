#include "weak_formulations.h"

// --- Matrices

void Mass(	libMesh::DenseMatrix<libMesh::Number>& Mass,
			unsigned int qp,
			const std::vector<std::vector<libMesh::Real> >& phi,
			const unsigned int n_dofs,
			const std::vector<libMesh::Real>& JxW,
			const libMesh::Number cte
			)
{
	for (unsigned int iii=0; iii< n_dofs; iii++)
	{
		for (unsigned int jjj=0; jjj< n_dofs; jjj++)
		{
			Mass(iii,jjj) += cte*JxW[qp]*phi[iii][qp]*phi[jjj][qp];
		}
	}
};

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
					)
{
	for (unsigned int iii=0; iii< n_dofs_sysA; iii++)
	{
		for (unsigned int jjj=0; jjj< n_dofs_sysB; jjj++)
		{
			Coupl(iii,jjj) += cte*JxW[qp]*phi_sysA[iii][qp]*phi_sysB[jjj][qp];
		}
	}
};
