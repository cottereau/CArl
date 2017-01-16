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

void L2_Coupling(	libMesh::DenseSubMatrix<libMesh::Number>& Coupl,
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

void H1_Coupling_Extra_Term(	libMesh::DenseSubMatrix<libMesh::Number>& Coupl,
								unsigned int qp,
								unsigned int C_i,
								unsigned int C_k,
								unsigned int n_components_A,
								unsigned int n_components_B,
								const std::vector<std::vector<libMesh::RealGradient> >& dphi_sysA,
								const std::vector<std::vector<libMesh::RealGradient> >& dphi_sysB,
								const unsigned int n_dofs_sysA,
								const unsigned int n_dofs_sysB,
								const std::vector<libMesh::Real>& JxW,
								const libMesh::Number cte
								)
{
	for (unsigned int iii=0; iii < n_dofs_sysA; iii++)
	{
		for (unsigned int jjj=0; jjj < n_dofs_sysB; jjj++)
		{
			for(unsigned int C_j=0; C_j<n_components_A; C_j++)
			{
				for(unsigned int C_l=0; C_l<n_components_B; C_l++)
				{
					Coupl(iii,jjj) += cte*JxW[qp]* dphi_sysA[iii][qp](C_j)*dphi_sysB[jjj][qp](C_l) *
											( kronecker_delta(C_i,C_k) * kronecker_delta(C_j,C_l) +
											  kronecker_delta(C_i,C_l) * kronecker_delta(C_j,C_k));
				}
			}
		}
	}
};
