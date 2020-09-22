#include "weak_formulations.h"

// --- Matrices

void Mass(  libMesh::DenseMatrix<libMesh::Number>& Mass,
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
void L2_Coupling( libMesh::DenseMatrix<libMesh::Number>& Coupl,
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

void L2_Coupling( libMesh::DenseSubMatrix<libMesh::Number>& Coupl,
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

void H1_Coupling_Extra_Term(  libMesh::DenseSubMatrix<libMesh::Number>& Coupl,
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
  int d_ik = kronecker_delta(C_i,C_k);
  int d_jl = 0;
  int d_il = 0;
  int d_jk = 0;
  double JW = JxW[qp];

  for (unsigned int iii=0; iii < n_dofs_sysA; iii++)
  {
    const libMesh::RealGradient& dphy_sysA_grad = dphi_sysA[iii][qp];
    for (unsigned int jjj=0; jjj < n_dofs_sysB; jjj++)
    {
      const libMesh::RealGradient& dphy_sysB_grad = dphi_sysB[jjj][qp];
      for(unsigned int C_j=0; C_j<n_components_A; C_j++)
      {
        d_jk = kronecker_delta(C_j,C_k);
        for(unsigned int C_l=0; C_l<n_components_B; C_l++)
        {
          d_il = kronecker_delta(C_i,C_l);
          d_jl = kronecker_delta(C_j,C_l);
          Coupl(iii,jjj) += cte*JW* dphy_sysA_grad(C_j)*dphy_sysB_grad(C_l) * ( d_ik * d_jl + d_il * d_jk ) ;
        }
      }
    }
  }
};
